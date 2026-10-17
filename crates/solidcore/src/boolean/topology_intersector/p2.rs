//! Topology intersector, part 2: construction, the registered surface/surface
//! intersection table, mate removal, vertex squeezing and gap fixing, face deletion,
//! the top-level intersect/insert/relate driver, the face/face intersection pass,
//! classification mate fixups and the edge/edge and edge/face coincidence passes.
//!
//! Upstream: `source/SMLib/src/SmTopologyIntersector.cpp`,
//! `source/SMLib/inc/SmTopologyIntersector.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`TopologyIntersector`] | `SmTopologyIntersector` (CO-PORTED data, see below) |
//! | [`TopologyIntersector::with_context`] | `SmTopologyIntersector(const SmContext&)` |
//! | [`TopologyIntersector::new`] | `SmTopologyIntersector(ctx, pBrep, pOther, tol, angTol)` |
//! | [`IntersectionTable`], [`IntersectionResults`] | `SmIntersectionTable`, `SmIntersectionResults` |
//! | [`TopologyIntersector::register_existing_ssi`] | `RegisterExistingSSI` |
//! | [`TopologyIntersector::check_for_squeeze_vertices`] | `CheckForSqueezeVertices` |
//! | [`TopologyIntersector::delete_faces`] | `DeleteFaces` |
//! | [`TopologyIntersector::fix_gaps`] | `FixGaps` |
//! | [`TopologyIntersector::remove_relationship`] | `RemoveRelationship` |
//! | [`TopologyIntersector::intersect_insert_relate`] | `IntersectInsertRelate` |
//! | [`TopologyIntersector::intersect_topology`] | `IntersectTopology` |
//! | [`TopologyIntersector::set_relationships`] | `SetRelationships` |
//! | [`TopologyIntersector::fixup_class_mates`] | `FixupClassMates` |
//! | [`TopologyIntersector::fixup_split_vertex_pair`] | `FixupSplitVertexPair` |
//! | [`TopologyIntersector::edge_edge_intersect`] | `EdgeEdgeIntersect` |
//! | [`TopologyIntersector::edge_face_coincidence`] | `EdgeFaceCoincidence` |
//!
//! Deviations:
//! - **CO-PORTED data.** The class and its header belong to
//!   `boolean.topology_intersector_p1` (open PR #146, not merged). This unit owns the
//!   constructors, so [`TopologyIntersector`] carries the minimal members these
//!   functions read or write: the option flags, the operand handles, the two mate
//!   maps (`m_pBToO`/`m_pOToB`, `None` until the full constructor allocates them),
//!   the tolerances, the intersection table and the proto topology manager. The
//!   one-line mate queries `GetOtherMate`, `GetThisMate`/`GetBrepMate` and
//!   `IsMatedPair` are co-ported as [`TopologyIntersector::other_mate`],
//!   [`TopologyIntersector::brep_mate`] and [`TopologyIntersector::is_mated_pair`];
//!   [`TopologyIntersector::insert_mate_pair`] is raw map storage (what
//!   `SmMapPtrToPtr::Insert` does), not `Relate`. `m_pSubset` is not modelled (no
//!   function here reads it).
//! - **B-rep access.** Topology objects are the caller's handles `O`
//!   ([`ClassObject`]); what these functions read from vertices, edges and faces is
//!   the adapter [`IntersectorTopology`] (the `ClassGeometry` hook pattern), which the
//!   owner of the B-reps implements with the ported topology functions.
//! - **Callees not ported yet** are private stubs returning [`SolidError::Unported`],
//!   tagged `UNPORTED(<unit>)` with the owning unit. The code after a stub is ported in
//!   full.
//! - Statuses upstream ignores (`SE`, unchecked returns) go through `ignored`; only
//!   [`SolidError::Unported`] propagates from them.
//! - Curves of the intersection table are owned values: the B-spline copies upstream
//!   makes (`new SmBSplineCurve(*SM_CAST_PTR(..))`, a slicing copy of the B-spline
//!   part) and clones of the caller's non-B-spline UV curves (upstream stores the
//!   caller's pointer).
//! - `FixupClassMates` restarts its scan after each successful
//!   `FixupSplitVertexPair`; each restart removes one interval, so the scan restarts at
//!   most `size` times. `EdgeEdgeIntersect` appends split pieces to the lists it walks;
//!   past [`MAX_EDGE_LIST`] entries it returns an error (upstream would keep growing).
//! - Uninitialised C++ reads (the proto topology manager before
//!   `IntersectInsertRelate` creates it) are [`SolidError::InvalidInput`].
//! - The `SM_DEBUG_CODE` and `SM_VALIDATE_TOPOLOGY` blocks are not ported; `SM_USE_TBB`
//!   is not defined upstream's reference build (its cache pre-warm is skipped).
//! - **Recursion.** None.

use std::collections::BTreeMap;

use crate::curve::curve::p2::{Curve, CurveKind};
use crate::curve::curve::p6::DroppedPoint;
use crate::error::{ignore_status, ignored, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::geom::global_solver::SolverOperation;
use crate::geom::solution_array::{ObjectKind, SolutionObject};
use crate::intersect::curve_class::p1::src_zone_tol_3d_of;
use crate::intersect::curve_class::CurveClassification;
use crate::intersect::point_class::{
    ClassGeometry, ClassObject, PointClassType, PointClassification,
};
use crate::kernel::Kernel;
use crate::math::{Point3, BIG_DOUBLE, EFF_ZERO_SQRT};
use crate::surface::surface::p2::{Surface, TsectCurveType};
use crate::topology::proto_topology::p1::proto::can_imprint;
use crate::topology::proto_topology::{ProtoBrep, ProtoEdgeId, ProtoTopologyManager};

const OP: &str = "boolean::topology_intersector::p2";

/// Upper bound on the working edge lists of
/// [`TopologyIntersector::edge_edge_intersect`].
pub const MAX_EDGE_LIST: usize = 1 << 20;

/// A null dereference or an uninitialised read upstream.
const fn undefined(reason: &'static str) -> SolidError {
    SolidError::InvalidInput {
        operation: OP,
        reason,
    }
}

/// The B-rep reads of these functions, implemented by the owner of the handles with
/// the ported topology functions (`topology.vertex`, `topology.edge`, `query.edge`).
pub trait IntersectorTopology<O> {
    /// `SmVertex::GetTolerance` / `SmEdge::GetTolerance`.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn topology_tolerance(&self, topo: O) -> Result<f64>;
    /// `SmVertex::GetPoint`.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn vertex_point(&self, vertex: O) -> Result<Point3>;
    /// `SmVertex::GetEdges`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn vertex_edges(&self, vertex: O) -> Result<Vec<O>>;
    /// `SmVertex::GetFaces`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn vertex_faces(&self, vertex: O) -> Result<Vec<O>>;
    /// `SmEdge::GetCurve`.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn edge_curve(&self, edge: O) -> Result<Option<&Curve>>;
    /// `SmEdge::GetInterval`.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn edge_interval(&self, edge: O) -> Result<Extent1d>;
    /// `SmEdge::CalculateBoundingBox(&box)` (the other arguments at their defaults).
    ///
    /// # Errors
    ///
    /// Those of the box computation.
    fn edge_bounding_box(&self, k: &Kernel, edge: O) -> Result<Extent3d>;
    /// `SmEdge::GetVertex`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn edge_vertex(&self, edge: O) -> Result<Option<O>>;
    /// `SmEdge::GetOtherVertex`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn edge_other_vertex(&self, edge: O, vertex: Option<O>) -> Result<Option<O>>;
    /// `SmEdge::GetFaces`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn edge_faces(&self, edge: O) -> Result<Vec<O>>;
    /// `SmEdge::IsWire`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn edge_is_wire(&self, edge: O) -> Result<bool>;
    /// `SmVertex::IsShellVertex`.
    ///
    /// # Errors
    ///
    /// Those of the topology walk.
    fn vertex_is_shell_vertex(&self, vertex: O) -> Result<bool>;
    /// `SmFace::GetSurface`.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn face_surface_of(&self, face: O) -> Result<Option<&Surface>>;
    /// `SmTol::GetZoneTol3d(pE)` of an edge.
    ///
    /// # Errors
    ///
    /// Those of the lookup.
    fn edge_zone_tol_3d(&self, k: &Kernel, edge: O) -> Result<f64>;
}

/// Which operand B-rep a call addresses (`m_pBrep` or `m_pOther`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operand {
    /// `m_pBrep`.
    Brep,
    /// `m_pOther`.
    Other,
}

/// The curves of one registered surface pair (`SmIntersectionResults` with curves).
#[derive(Debug, Clone, Default)]
pub struct SsiCurves {
    /// `m_p3DCurves`.
    pub curves_3d: Vec<Option<Curve>>,
    /// `m_pUVCurves1`.
    pub uv_curves1: Vec<Option<Curve>>,
    /// `m_pUVCurves2`.
    pub uv_curves2: Vec<Option<Curve>>,
    /// `m_pCurveTypes`.
    pub curve_types: Vec<TsectCurveType>,
    /// `m_pDeviations`.
    pub deviations: Vec<f64>,
}

/// `SmIntersectionResults`: `None` when no curves were registered.
#[derive(Debug, Clone, Default)]
pub struct IntersectionResults {
    /// The registered curves (the arrays are allocated only for a nonempty set).
    pub curves: Option<SsiCurves>,
}

/// `SmIntersectionTable`: registered surface pairs and their results. The surface
/// lists can be longer than the results when a registration failed part way.
#[derive(Debug, Clone)]
pub struct IntersectionTable<O> {
    /// `m_vBrepSurfaces`.
    pub brep_surfaces: Vec<Option<O>>,
    /// `m_vOtherSurfaces`.
    pub other_surfaces: Vec<Option<O>>,
    /// `m_vResults`.
    pub results: Vec<IntersectionResults>,
}

impl<O> Default for IntersectionTable<O> {
    fn default() -> Self {
        Self {
            brep_surfaces: Vec::new(),
            other_surfaces: Vec::new(),
            results: Vec::new(),
        }
    }
}

/// The two mate maps (`m_pBToO`, `m_pOToB`).
#[derive(Debug, Clone)]
pub struct MateMaps<O> {
    /// `m_pBToO`: primary entity to other entity.
    pub b_to_o: BTreeMap<O, O>,
    /// `m_pOToB`: other entity to primary entity.
    pub o_to_b: BTreeMap<O, O>,
}

impl<O> Default for MateMaps<O> {
    fn default() -> Self {
        Self {
            b_to_o: BTreeMap::new(),
            o_to_b: BTreeMap::new(),
        }
    }
}

/// The outputs of [`TopologyIntersector::check_for_squeeze_vertices`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SqueezeOutputs<O> {
    /// `rDeletedV`.
    pub deleted_vertices: Vec<O>,
    /// `rSurvivingV`.
    pub surviving_vertices: Vec<O>,
    /// `rDeletedE`.
    pub deleted_edges: Vec<O>,
}

impl<O> Default for SqueezeOutputs<O> {
    fn default() -> Self {
        Self {
            deleted_vertices: Vec::new(),
            surviving_vertices: Vec::new(),
            deleted_edges: Vec::new(),
        }
    }
}

/// The outputs of [`TopologyIntersector::fix_gaps`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FixGapsOutputs<O> {
    /// `rbModifiedTopology`.
    pub modified_topology: bool,
    /// `rDeletedVertices`.
    pub deleted_vertices: Vec<O>,
    /// `rDeletedEdges` (a `NULL` entry is `None`).
    pub deleted_edges: Vec<Option<O>>,
}

/// The outputs of [`TopologyIntersector::edge_edge_intersect`].
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeEdgeOutputs<O> {
    /// `rSplitEdges1`.
    pub split_edges1: Vec<O>,
    /// `rNewEdges1`.
    pub new_edges1: Vec<O>,
    /// `rSplitEdges2`.
    pub split_edges2: Vec<O>,
    /// `rNewEdges2`.
    pub new_edges2: Vec<O>,
    /// `rNewVerts`.
    pub new_verts: Vec<O>,
    /// `rdMaxEdgeGap`.
    pub max_edge_gap: f64,
    /// `rdMinUnstitchedEdgeGap`.
    pub min_unstitched_edge_gap: f64,
}

/// The outputs of [`TopologyIntersector::edge_face_coincidence`].
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeFaceOutputs<O> {
    /// `rImprintedEdges`.
    pub imprinted_edges: Vec<O>,
    /// `rProcessedFaces`.
    pub processed_faces: Vec<O>,
    /// `rNewFaces`.
    pub new_faces: Vec<O>,
    /// `rSplitFaces`.
    pub split_faces: Vec<O>,
    /// `dMaxEdgeGap`.
    pub max_edge_gap: f64,
    /// `dThisMinUnstitched`.
    pub this_min_unstitched: f64,
}

/// `SmTopologyIntersector` (the members used by this unit; see the module header).
#[derive(Debug, Clone)]
pub struct TopologyIntersector<O> {
    /// `m_bSkipSSI`.
    pub skip_ssi: bool,
    /// `m_bSkipManifoldEdges`.
    pub skip_manifold_edges: bool,
    /// `m_bIntersectTangentEdges`.
    pub intersect_tangent_edges: bool,
    /// `m_lIntersectLaminaEdges`: 0 none, 1 both B-reps, 2 `m_pBrep`, 3 `m_pOther`.
    pub intersect_lamina_edges: u32,
    /// `m_bIntersectWireEdges`.
    pub intersect_wire_edges: bool,
    /// `m_bIntersectShellVertices`.
    pub intersect_shell_vertices: bool,
    /// `m_bIgnoreTangency`.
    pub ignore_tangency: bool,
    /// `m_bDoGluing`.
    pub do_gluing: bool,
    /// `m_bDoingFillets`.
    pub doing_fillets: bool,
    /// `m_bTopologyDeleted`.
    pub topology_deleted: bool,
    /// `m_pBrep`.
    pub brep: Option<O>,
    /// `m_pOther`.
    pub other: Option<O>,
    /// `m_pBToO` / `m_pOToB` (`None` when not allocated).
    pub mates: Option<MateMaps<O>>,
    /// `m_bSwappedOrder`.
    pub swapped_order: bool,
    /// `m_dThisApproxTol3d`.
    pub this_approx_tol_3d: f64,
    /// `m_dThisAngTolRad`.
    pub this_ang_tol_rad: f64,
    /// `m_pIntersectionTable`.
    pub intersection_table: Option<IntersectionTable<O>>,
    /// `m_pProtoTopoMgr` (`None` is the uninitialised or deleted pointer).
    pub proto_topo_mgr: Option<ProtoTopologyManager<O>>,
}

/// `new SmBSplineCurve(*SM_CAST_PTR(SmBSplineCurve, pCurve))`: the slicing copy of the
/// B-spline part of a B-spline-derived curve, `None` when the cast is null.
fn bspline_slice(curve: Option<&Curve>) -> Option<Curve> {
    let base = match curve?.kind() {
        CurveKind::BSpline(c) => c.clone(),
        CurveKind::Line(c) => c.base().clone(),
        CurveKind::Ellipse(c) => c.base().clone(),
        CurveKind::Hyperbola(c) => c.base().clone(),
        CurveKind::Parabola(c) => c.base().clone(),
        _ => return None,
    };
    Some(Curve::from(base))
}

/// `SM_CAST_PTR(SmBSplineCurve, c)` is not null.
fn is_bspline(curve: Option<&Curve>) -> bool {
    curve.is_some_and(|c| {
        matches!(
            c.kind(),
            CurveKind::BSpline(_)
                | CurveKind::Line(_)
                | CurveKind::Ellipse(_)
                | CurveKind::Hyperbola(_)
                | CurveKind::Parabola(_)
        )
    })
}

// ---- stubs for callees in units not ported yet ----------------------------------------

/// `SmBrep::CombineCoincidentVertices` of operand `which`: the deleted edge.
fn brep_combine_coincident_vertices<O>(_which: Operand, _keep: O, _gone: O) -> Result<Option<O>> {
    // UNPORTED(topology.brep_p5): SmBrep::CombineCoincidentVertices
    Err(SolidError::unported(
        "SmBrep::CombineCoincidentVertices",
        "topology.brep_p5",
    ))
}

/// `SmFace::GetEdges`.
fn face_get_edges<O>(_face: O) -> Result<Vec<O>> {
    // UNPORTED(topology.face_p4): SmFace::GetEdges
    Err(SolidError::unported("SmFace::GetEdges", "topology.face_p4"))
}

/// `SmFace::GetVertices`.
fn face_get_vertices<O>(_face: O) -> Result<Vec<O>> {
    // UNPORTED(topology.face_p4): SmFace::GetVertices
    Err(SolidError::unported(
        "SmFace::GetVertices",
        "topology.face_p4",
    ))
}

/// `SmFace::GetBrep`.
fn face_get_brep<O>(_face: O) -> Result<Option<O>> {
    // UNPORTED(topology.face_p4): SmFace::GetBrep
    Err(SolidError::unported("SmFace::GetBrep", "topology.face_p4"))
}

/// `SmBrep::DeleteFace(pF, FALSE)`.
fn brep_delete_face<O>(_brep: O, _face: O) -> Result<()> {
    // UNPORTED(topology.brep_p3): SmBrep::DeleteFace
    Err(SolidError::unported(
        "SmBrep::DeleteFace",
        "topology.brep_p3",
    ))
}

/// `pE->IsWire()`, then `pE->GetBrep()->DeleteEdge(pE)` live in the deleted-face
/// clean-up; `SmBrep::DeleteEdge`.
fn brep_delete_edge<O>(_edge: O) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::DeleteEdge
    Err(SolidError::unported(
        "SmBrep::DeleteEdge",
        "topology.brep_p4",
    ))
}

/// `SmBrep::DeleteVertex`.
fn brep_delete_vertex<O>(_vertex: O) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::DeleteVertex
    Err(SolidError::unported(
        "SmBrep::DeleteVertex",
        "topology.brep_p4",
    ))
}

/// `SmBrep::GetFaces`.
fn brep_get_faces<O>(_brep: Option<O>) -> Result<Vec<O>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetFaces
    Err(SolidError::unported("SmBrep::GetFaces", "topology.brep_p1"))
}

/// `sm_IntersectShellVertices(pBrep, pOther, bSwapped, pMgr)` (status ignored).
fn intersect_shell_vertices<O>(_brep: Option<O>, _other: Option<O>, _swapped: bool) -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): sm_IntersectShellVertices
    Err(SolidError::unported(
        "sm_IntersectShellVertices",
        "boolean.topology_intersector_p1",
    ))
}

/// `sm_IntersectWireLaminaEdges(pBrep, pOther, bDoLamina, bSwapped, pMgr)`.
fn intersect_wire_lamina_edges<O>(
    _brep: Option<O>,
    _other: Option<O>,
    _do_lamina: bool,
    _swapped: bool,
) -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): sm_IntersectWireLaminaEdges
    Err(SolidError::unported(
        "sm_IntersectWireLaminaEdges",
        "boolean.topology_intersector_p1",
    ))
}

/// The registered intersection of a face pair (`CheckIfHaveExistingSSI`): whether one
/// exists.
fn check_if_have_existing_ssi<O>(
    _brep_surface: Option<O>,
    _other_surface: Option<O>,
) -> Result<bool> {
    // UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::CheckIfHaveExistingSSI
    Err(SolidError::unported(
        "SmTopologyIntersector::CheckIfHaveExistingSSI",
        "boolean.topology_intersector_p1",
    ))
}

/// The outputs of `SmFace::FaceIntersect` read here.
struct FaceIntersection<O> {
    curve_count: usize,
    int_pts: Vec<Point3>,
    ccs1: Vec<CurveClassification<O>>,
    ccs2: Vec<CurveClassification<O>>,
}

/// `SmFace::FaceIntersect(ctx, pF2, dTol, 20.0, ..., &bAreCoincident)`.
fn face_intersect<O>(
    _f1: O,
    _f2: O,
    _tol: f64,
    _ang_tol_deg: f64,
    _are_coincident: bool,
) -> Result<FaceIntersection<O>> {
    // UNPORTED(topology.face_p3): SmFace::FaceIntersect
    Err(SolidError::unported(
        "SmFace::FaceIntersect",
        "topology.face_p3",
    ))
}

/// `SmFace::Point3DClassify(pt, tol, bDoBoundary, rPC)`.
fn face_point_3d_classify<O>(
    _face: O,
    _point: Point3,
    _tol: f64,
    _do_boundary: bool,
    _pc: &mut PointClassification<O>,
) -> Result<()> {
    // UNPORTED(query.face_p1): SmFace::Point3DClassify
    Err(SolidError::unported(
        "SmFace::Point3DClassify",
        "query.face_p1",
    ))
}

/// `SmTopologyIntersector::ResolveIntersections`.
fn resolve_intersections() -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::ResolveIntersections
    Err(SolidError::unported(
        "SmTopologyIntersector::ResolveIntersections",
        "boolean.topology_intersector_p1",
    ))
}

/// `SmTopologyIntersector::ImprintIntersections`.
fn imprint_intersections() -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::ImprintIntersections
    Err(SolidError::unported(
        "SmTopologyIntersector::ImprintIntersections",
        "boolean.topology_intersector_p1",
    ))
}

/// `SmTopologyIntersector::Relate(pTopo1, pTopo2)` (status ignored by the caller).
fn relate<O>(_topo1: O, _topo2: O) -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::Relate
    Err(SolidError::unported(
        "SmTopologyIntersector::Relate",
        "boolean.topology_intersector_p1",
    ))
}

/// `FixupEndPointClassMates` / `FixupMidPointClassMates`.
fn fixup_point_class_mates(symbol: &'static str) -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p1): SmTopologyIntersector::FixupEndPointClassMates
    Err(SolidError::unported(
        symbol,
        "boolean.topology_intersector_p1",
    ))
}

/// `sm_TestAndSplitTwoEdges`: `(dDistToEdge, dMinUnStitch, sNewV, sNewE1, sNewE2)`.
#[allow(clippy::type_complexity)]
fn test_and_split_two_edges<O>(
    _e1: O,
    _e2: O,
    _tol: f64,
) -> Result<(f64, f64, Vec<O>, Vec<O>, Vec<O>)> {
    // UNPORTED(boolean.topology_intersector_p1): sm_TestAndSplitTwoEdges
    Err(SolidError::unported(
        "sm_TestAndSplitTwoEdges",
        "boolean.topology_intersector_p1",
    ))
}

/// `SmTol::UpdateObjectTolerance(SmTopology*, ..)`: its topology overload is not in
/// the port of `geom.tol` (see the note in `geom/tol.rs`).
fn update_object_tolerance<O>(_object: O, _dist: f64) -> Result<()> {
    // UNPORTED(geom.tol): SmTol::UpdateObjectTolerance -- justified: plan/findings.md "boolean.topology_intersector_p2: SmTol::UpdateObjectTolerance topology overload"
    Err(SolidError::unported(
        "SmTol::UpdateObjectTolerance",
        "geom.tol",
    ))
}

/// `SmBrep::MakeEdgeInFace`: `(pNewEdge, pNewLoop, pNewFace)`.
#[allow(clippy::type_complexity)]
fn brep_make_edge_in_face<O>(
    _face: O,
    _v: O,
    _other_v: O,
) -> Result<(Option<O>, Option<O>, Option<O>)> {
    // UNPORTED(topology.brep_p5): SmBrep::MakeEdgeInFace
    Err(SolidError::unported(
        "SmBrep::MakeEdgeInFace",
        "topology.brep_p5",
    ))
}

/// `SmBrep::GlueEdgesGeneral`.
fn brep_glue_edges_general<O>(_keep: O, _gone: O) -> Result<()> {
    // UNPORTED(topology.brep_p7): SmBrep::GlueEdgesGeneral
    Err(SolidError::unported(
        "SmBrep::GlueEdgesGeneral",
        "topology.brep_p7",
    ))
}

/// The objects of a coincidence solution (unused here).
#[derive(Debug, Clone, Copy, PartialEq)]
struct NoObject;

impl SolutionObject for NoObject {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Curve
    }
}

/// `SmFace::CalculateUVDomainFromUVTrimCurves`.
fn face_uv_domain_from_trim_curves<O>(_face: O) -> Result<Extent2d> {
    // UNPORTED(topology.face_p1): SmFace::CalculateUVDomainFromUVTrimCurves
    Err(SolidError::unported(
        "SmFace::CalculateUVDomainFromUVTrimCurves",
        "topology.face_p1",
    ))
}

impl<O: ClassObject + Ord> TopologyIntersector<O> {
    /// `SmTopologyIntersector(const SmContext&)`: no operands, no mate maps, zero
    /// tolerances. The proto topology manager is left uninitialised upstream.
    #[must_use]
    pub const fn with_context() -> Self {
        Self {
            skip_ssi: false,
            skip_manifold_edges: true,
            intersect_tangent_edges: true,
            intersect_lamina_edges: 0,
            intersect_wire_edges: true,
            intersect_shell_vertices: true,
            ignore_tangency: false,
            do_gluing: false,
            doing_fillets: false,
            topology_deleted: false,
            brep: None,
            other: None,
            mates: None,
            swapped_order: false,
            this_approx_tol_3d: 0.0,
            this_ang_tol_rad: 0.0,
            intersection_table: None,
            proto_topo_mgr: None,
        }
    }

    /// `SmTopologyIntersector(ctx, pBrep, pOther, dThisApproxTol3d, dThisAngTolRad)`:
    /// the operands, the tolerances and two empty mate maps.
    #[must_use]
    pub fn new(
        brep: Option<O>,
        other: Option<O>,
        this_approx_tol_3d: f64,
        this_ang_tol_rad: f64,
    ) -> Self {
        let mut ti = Self::with_context();
        ti.brep = brep;
        ti.other = other;
        ti.this_approx_tol_3d = this_approx_tol_3d;
        ti.this_ang_tol_rad = this_ang_tol_rad;
        ti.mates = Some(MateMaps::default());
        ti
    }

    // ---- co-ported mate queries (boolean.topology_intersector_p1) -----------------------

    /// `GetOtherMate` (CO-PORTED): the other entity mated to a primary entity.
    #[must_use]
    pub fn other_mate(&self, brep_entity: Option<O>) -> Option<O> {
        let maps = self.mates.as_ref()?;
        maps.b_to_o.get(&brep_entity?).copied()
    }

    /// `GetBrepMate` / `GetThisMate` (CO-PORTED): the primary entity mated to an other
    /// entity.
    #[must_use]
    pub fn brep_mate(&self, other_entity: Option<O>) -> Option<O> {
        let maps = self.mates.as_ref()?;
        maps.o_to_b.get(&other_entity?).copied()
    }

    /// `IsMatedPair` (CO-PORTED).
    #[must_use]
    pub fn is_mated_pair(&self, brep_entity: Option<O>, other_entity: Option<O>) -> bool {
        if brep_entity.is_none() || other_entity.is_none() {
            return false;
        }
        other_entity == self.other_mate(brep_entity)
    }

    /// Raw storage of a mate pair in both maps (`m_pBToO->Insert`, `m_pOToB->Insert`),
    /// not `Relate`: no checks and no tolerance updates. A no-op without maps.
    pub fn insert_mate_pair(&mut self, brep_entity: O, other_entity: O) {
        if let Some(maps) = self.mates.as_mut() {
            maps.b_to_o.insert(brep_entity, other_entity);
            maps.o_to_b.insert(other_entity, brep_entity);
        }
    }

    // ---- owned functions ---------------------------------------------------------------

    /// `RemoveRelationship`: removes a mated pair from both maps.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when the two are not mated.
    pub fn remove_relationship(&mut self, brep_entity: O, other_entity: O) -> Result<()> {
        if !self.is_mated_pair(Some(brep_entity), Some(other_entity)) {
            return Err(SolidError::sm_err(OP));
        }
        // a mated pair implies allocated maps
        if let Some(maps) = self.mates.as_mut() {
            maps.b_to_o.remove(&brep_entity);
            maps.o_to_b.remove(&other_entity);
        }
        Ok(())
    }

    /// `RegisterExistingSSI`: registers the intersection curves of a surface pair,
    /// copying the B-spline curves. A `None` surface turns off the non-registered
    /// intersections of the other.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when a 3D curve is not a B-spline (the surfaces stay registered without
    /// results, as upstream); [`SolidError::InvalidInput`] when a UV curve list is
    /// shorter than the 3D curve list (an out-of-range read upstream).
    #[allow(clippy::too_many_arguments)]
    pub fn register_existing_ssi(
        &mut self,
        brep_surface: Option<O>,
        other_surface: Option<O>,
        curves_3d: &[Option<&Curve>],
        uv_curves1: &[Option<&Curve>],
        uv_curves2: &[Option<&Curve>],
        curve_types: &[TsectCurveType],
        deviations: &[f64],
    ) -> Result<()> {
        if !curves_3d.is_empty()
            && (uv_curves1.len() < curves_3d.len() || uv_curves2.len() < curves_3d.len())
        {
            return Err(undefined("UV curve list shorter than the 3D curve list"));
        }
        let table = self
            .intersection_table
            .get_or_insert_with(IntersectionTable::default);
        table.brep_surfaces.push(brep_surface);
        table.other_surfaces.push(other_surface);
        let mut result = IntersectionResults::default();
        if !curves_3d.is_empty() {
            let mut sets = SsiCurves {
                curves_3d: curves_3d.iter().map(|c| c.cloned()).collect(),
                uv_curves1: uv_curves1.iter().map(|c| c.cloned()).collect(),
                uv_curves2: uv_curves2.iter().map(|c| c.cloned()).collect(),
                curve_types: curve_types.to_vec(),
                deviations: deviations.to_vec(),
            };
            for ii in 0..sets.curves_3d.len() {
                let c3 = curves_3d.get(ii).copied().flatten();
                if is_bspline(c3) {
                    if let Some(slot) = sets.curves_3d.get_mut(ii) {
                        *slot = bspline_slice(c3);
                    }
                } else {
                    return Err(SolidError::sm_err(OP));
                }
                let u1 = uv_curves1.get(ii).copied().flatten();
                if is_bspline(u1) {
                    if let Some(slot) = sets.uv_curves1.get_mut(ii) {
                        *slot = bspline_slice(u1);
                    }
                }
                let u2 = uv_curves2.get(ii).copied().flatten();
                if is_bspline(u2) {
                    if let Some(slot) = sets.uv_curves2.get_mut(ii) {
                        *slot = bspline_slice(u2);
                    }
                }
            }
            result.curves = Some(sets);
        }
        table.results.push(result);
        Ok(())
    }

    /// `CheckForSqueezeVertices`: resolves a conflicting earlier relationship of either
    /// vertex by combining the vertices it mated, unless that relationship was better.
    ///
    /// # Errors
    ///
    /// `SM_ERR` (`Conflicting relationship, previous one is better`); those of
    /// `RemoveRelationship`; [`SolidError::Unported`] for
    /// `SmBrep::CombineCoincidentVertices`.
    pub fn check_for_squeeze_vertices(
        &mut self,
        geometry: &impl IntersectorTopology<O>,
        brep_vertex: O,
        other_vertex: O,
        out: &mut SqueezeOutputs<O>,
    ) -> Result<()> {
        // cbi breakpoint checks (`SER(SM_SUCCESS)` is a no-op), tolerances still read
        let _ = geometry.topology_tolerance(brep_vertex)?;
        let _ = geometry.topology_tolerance(other_vertex)?;

        // Get any already-existing vertex relationships.
        let other_v = self.other_mate(Some(brep_vertex));
        let brep_v = self.brep_mate(Some(other_vertex));

        // Check for a pre-existing, conflicting relationship.
        if let Some(other_v) = other_v.filter(|&v| v != other_vertex) {
            let _ = geometry.topology_tolerance(other_v)?;
            let p = geometry.vertex_point(brep_vertex)?;
            let dist_old = p.distance(geometry.vertex_point(other_v)?);
            let dist_new = p.distance(geometry.vertex_point(other_vertex)?);
            if dist_new > dist_old + EFF_ZERO_SQRT {
                return Err(SolidError::sm_err(OP));
            }
            self.remove_relationship(brep_vertex, other_v)?;
            let deleted_other_edge =
                brep_combine_coincident_vertices(Operand::Other, other_vertex, other_v)?;
            if let Some(e) = deleted_other_edge {
                out.deleted_edges.push(e);
                if let Some(deleted_brep_edge) = self.brep_mate(Some(e)) {
                    self.remove_relationship(deleted_brep_edge, e)?;
                }
            }
            out.surviving_vertices.push(other_vertex);
            out.deleted_vertices.push(other_v);
        }

        // Check for a pre-existing, conflicting relationship.
        if let Some(brep_v) = brep_v.filter(|&v| v != brep_vertex) {
            let _ = geometry.topology_tolerance(brep_v)?;
            let p = geometry.vertex_point(other_vertex)?;
            let dist_old = p.distance(geometry.vertex_point(brep_v)?);
            let dist_new = p.distance(geometry.vertex_point(brep_vertex)?);
            if dist_new > dist_old + EFF_ZERO_SQRT {
                return Err(SolidError::sm_err(OP));
            }
            self.remove_relationship(brep_v, other_vertex)?;
            let deleted_brep_edge =
                brep_combine_coincident_vertices(Operand::Brep, brep_vertex, brep_v)?;
            if let Some(e) = deleted_brep_edge {
                out.deleted_edges.push(e);
                if let Some(deleted_other_edge) = self.other_mate(Some(e)) {
                    self.remove_relationship(e, deleted_other_edge)?;
                }
            }
            out.surviving_vertices.push(brep_vertex);
            out.deleted_vertices.push(brep_v);
        }
        Ok(())
    }

    /// `DeleteFaces`: removes the faces from their B-rep, then the edges left as wires
    /// (unless kept) and the vertices left as shell vertices, unrelating each.
    ///
    /// # Errors
    ///
    /// [`SolidError::Unported`] for the face and B-rep edits (topology units).
    pub fn delete_faces(
        &mut self,
        geometry: &impl IntersectorTopology<O>,
        faces: &[O],
        keep_as_wires: Option<&[O]>,
    ) -> Result<()> {
        let mut faces_vertices: Vec<O> = Vec::new();
        let mut faces_edges: Vec<O> = Vec::new();
        let mut brep: Option<O> = None;

        for &f in faces {
            for e in face_get_edges(f)? {
                if !faces_edges.contains(&e) {
                    faces_edges.push(e);
                }
            }
            for v in face_get_vertices(f)? {
                if !faces_vertices.contains(&v) {
                    faces_vertices.push(v);
                }
            }
            if let Some(of) = self.other_mate(Some(f)) {
                ignore_status(self.remove_relationship(f, of))?;
            }
            brep = face_get_brep(f)?;
            let b = brep.ok_or(undefined("face without B-rep"))?;
            brep_delete_face(b, f)?;
        }

        if brep.is_none() {
            return Ok(());
        }

        // iter every edge of every deleted face looking for wire edges to delete
        for &e in &faces_edges {
            if geometry.edge_is_wire(e)? {
                if keep_as_wires.is_some_and(|k| k.contains(&e)) {
                    continue;
                }
                if let Some(oe) = self.other_mate(Some(e)) {
                    ignore_status(self.remove_relationship(e, oe))?;
                }
                brep_delete_edge(e)?;
            }
        }
        // iter every vertex of every deleted face looking for ShellVertices to delete
        for &v in &faces_vertices {
            if geometry.vertex_is_shell_vertex(v)? {
                if let Some(ov) = self.other_mate(Some(v)) {
                    ignore_status(self.remove_relationship(v, ov))?;
                }
                brep_delete_vertex(v)?;
            }
        }
        Ok(())
    }

    /// `FixGaps`: squeezes together pairs of single primary vertices closer than their
    /// summed tolerances (times 1.01) that share no edge longer than that.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` for a candidate without other mate; those of the topology
    /// reads; [`SolidError::Unported`] for `SmBrep::CombineCoincidentVertices`.
    pub fn fix_gaps(
        &mut self,
        geometry: &impl IntersectorTopology<O>,
        single_vertices: &[O],
        out: &mut FixGapsOutputs<O>,
    ) -> Result<()> {
        out.modified_topology = false;
        out.deleted_vertices.clear();
        out.deleted_edges.clear();
        let n = single_vertices.len();
        for ii in 0..n {
            let Some(&v1) = single_vertices.get(ii) else {
                break;
            };
            for jj in ii + 1..n {
                let Some(&v2) = single_vertices.get(jj) else {
                    break;
                };
                // Increase tol just a bit.  Was 10.0.  [B498]
                let tol_factor = 1.01;
                let tol = tol_factor
                    * (geometry.topology_tolerance(v1)? + geometry.topology_tolerance(v2)?);
                if geometry
                    .vertex_point(v1)?
                    .distance(geometry.vertex_point(v2)?)
                    < tol
                {
                    let edges1 = geometry.vertex_edges(v1)?;
                    let edges2 = geometry.vertex_edges(v2)?;
                    let mut found_common_edge = false;
                    for &test in &edges1 {
                        if edges2.contains(&test) {
                            let curve = geometry
                                .edge_curve(test)?
                                .ok_or(undefined("edge without curve"))?;
                            let interval = geometry.edge_interval(test)?;
                            let (edge_len, _) = curve.approximate_length(&interval, 10, false)?;
                            if edge_len > tol {
                                found_common_edge = true;
                            }
                        }
                    }
                    if !found_common_edge {
                        let other_v1 = self
                            .other_mate(Some(v1))
                            .ok_or(SolidError::null_pointer(OP))?;
                        let other_v2 = self
                            .other_mate(Some(v2))
                            .ok_or(SolidError::null_pointer(OP))?;
                        let deleted_brep_edge =
                            brep_combine_coincident_vertices(Operand::Brep, v2, v1)?;
                        out.deleted_vertices.push(v1);
                        out.deleted_edges.push(deleted_brep_edge);
                        let deleted_other_edge =
                            brep_combine_coincident_vertices(Operand::Other, other_v2, other_v1)?;
                        ignore_status(self.remove_relationship(v1, other_v1))?;
                        if let (Some(db), Some(dother)) = (deleted_brep_edge, deleted_other_edge) {
                            if self.is_mated_pair(Some(db), Some(dother)) {
                                ignore_status(self.remove_relationship(db, dother))?;
                            } else {
                                if let Some(mate) = self.other_mate(Some(db)) {
                                    ignore_status(self.remove_relationship(db, mate))?;
                                }
                                if let Some(mate) = self.brep_mate(Some(dother)) {
                                    ignore_status(self.remove_relationship(mate, dother))?;
                                }
                            }
                        }
                        out.modified_topology = true;
                        break;
                    }
                }
            }
        }
        Ok(())
    }

    /// `IntersectInsertRelate`: creates the proto topology manager, intersects the
    /// topology, resolves and imprints the intersections and sets the relationships;
    /// the manager is deleted afterwards.
    ///
    /// # Errors
    ///
    /// Those of `IntersectTopology` and `ResolveIntersections`; those of the manager
    /// constructor; [`SolidError::Unported`] for the unported callees.
    pub fn intersect_insert_relate<G>(&mut self, k: &Kernel, geometry: &mut G, ti: O) -> Result<()>
    where
        G: ProtoBrep<O> + IntersectorTopology<O>,
    {
        let mut capture_distance = 0.01;
        let topo_int_tol = self.this_approx_tol_3d;
        if capture_distance < topo_int_tol * 2000.0 {
            capture_distance = topo_int_tol * 2000.0;
        }
        let brep1 = self.brep.ok_or(undefined("no primary B-rep"))?;
        let brep2 = self.other.ok_or(undefined("no other B-rep"))?;
        self.proto_topo_mgr = Some(ProtoTopologyManager::new(
            brep1,
            brep2,
            ti,
            self.this_approx_tol_3d,
            capture_distance,
            0.0,
        )?);

        // Intersect Topology; if that didn't work, quit before we mess up the Breps.
        self.intersect_topology(k, geometry)?;
        // Resolve Intersections
        resolve_intersections()?;
        // Imprint Intersections (SE_MSG: status ignored)
        ignore_status(imprint_intersections())?;
        // Set relationship maps.
        let stat = self.set_relationships();
        let _ = self.proto_topo_mgr.take();
        stat
    }

    /// `IntersectTopology`: shell vertex, wire/lamina edge and face/face intersections
    /// into the proto topology, then the proto edge/edge splits and the removal of
    /// coincident edge-less proto vertices.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] without a proto topology manager (an
    /// uninitialised pointer upstream); those of the checked callees;
    /// [`SolidError::Unported`] for the unported callees.
    pub fn intersect_topology<G>(&mut self, k: &Kernel, geometry: &mut G) -> Result<()>
    where
        G: ProtoBrep<O> + IntersectorTopology<O>,
    {
        self.proto_topo_mgr
            .as_mut()
            .ok_or(undefined("no proto topology manager"))?
            .reset();
        let faces1 = brep_get_faces(self.brep)?;
        let faces2 = brep_get_faces(self.other)?;

        // Process ShellVertex Intersections (statuses ignored).
        ignore_status(intersect_shell_vertices(self.brep, self.other, false))?;
        ignore_status(intersect_shell_vertices(self.other, self.brep, true))?;

        // Wire/Lamina Edge intersections (statuses ignored).
        let do_lamina = self.intersect_lamina_edges == 1 || self.intersect_lamina_edges == 2;
        ignore_status(intersect_wire_lamina_edges(
            self.brep, self.other, do_lamina, false,
        ))?;
        let do_lamina = self.intersect_lamina_edges == 1 || self.intersect_lamina_edges == 3;
        ignore_status(intersect_wire_lamina_edges(
            self.other, self.brep, do_lamina, true,
        ))?;

        let mut coinc_pes: Vec<ProtoEdgeId> = Vec::new();
        for &f1 in &faces1 {
            for &f2 in &faces2 {
                let s1 = geometry
                    .face_surface(f1)
                    .ok_or(undefined("bad face handle"))?;
                let s2 = geometry
                    .face_surface(f2)
                    .ok_or(undefined("bad face handle"))?;
                if check_if_have_existing_ssi(s1, s2)? {
                    continue;
                }
                let tol = self.this_approx_tol_3d;
                coinc_pes.clear();
                let surf1 = s1
                    .and_then(|s| geometry.surface(s))
                    .ok_or(undefined("face without surface"))?;
                let surf2 = s2
                    .and_then(|s| geometry.surface(s))
                    .ok_or(undefined("face without surface"))?;
                let are_coincident = surf1
                    .coincidence_check(k, surf2, tol, false)?
                    .are_coincident;
                let ang_tol_deg = 20.0;
                let mut fi = face_intersect(f1, f2, tol, ang_tol_deg, are_coincident)?;

                if fi.curve_count == 0 && fi.int_pts.is_empty() {
                    continue;
                }
                let ti_tol = self.this_approx_tol_3d;
                let mgr = self
                    .proto_topo_mgr
                    .as_mut()
                    .ok_or(undefined("no proto topology manager"))?;
                for (cc1, cc2) in fi.ccs1.iter_mut().zip(fi.ccs2.iter_mut()) {
                    ignore_status(cc1.set_point_class_parameters(k, &*geometry, false))?;
                    ignore_status(cc2.set_point_class_parameters(k, &*geometry, false))?;
                    for ll in 0..cc1.size() {
                        let mut ci1 = *cc1.at(ll)?;
                        let mut ci2 = *cc2.at(ll)?;
                        let found = ignored(mgr.find_or_create_proto_edge(
                            k, geometry, ti_tol, cc1, &mut ci1, cc2, &mut ci2,
                        ))?;
                        *cc1.at_mut(ll)? = ci1;
                        *cc2.at_mut(ll)? = ci2;
                        if let Some((Some(pe), _)) = found {
                            if are_coincident && !coinc_pes.contains(&pe) {
                                coinc_pes.push(pe);
                            }
                        }
                    }
                }

                let f1_zone_tol = geometry
                    .topology_zone_tol_3d(f1)
                    .ok_or(undefined("face tolerance"))?;
                let f2_zone_tol = geometry
                    .topology_zone_tol_3d(f2)
                    .ok_or(undefined("face tolerance"))?;
                for &pt in &fi.int_pts {
                    let mut pc1 = PointClassification::new(f1_zone_tol);
                    let mut pc2 = PointClassification::new(f2_zone_tol);
                    let stat1 =
                        ignored(face_point_3d_classify(f1, pt, f1_zone_tol, true, &mut pc1))?;
                    let stat2 =
                        ignored(face_point_3d_classify(f2, pt, f2_zone_tol, true, &mut pc2))?;
                    if stat1.is_none() || stat2.is_none() {
                        continue;
                    }
                    if !can_imprint(pc1.point_class()) || !can_imprint(pc2.point_class()) {
                        continue;
                    }
                    ignore_status(
                        mgr.find_or_create_proto_vertex_pcs(k, &*geometry, &mut pc1, &mut pc2),
                    )?;
                }

                if are_coincident && !coinc_pes.is_empty() {
                    let pf = ignored(mgr.create_proto_face(Some(f1), Some(f2), &coinc_pes))?;
                    let pf = pf.ok_or(undefined("NULL proto face"))?;
                    for &pe in &coinc_pes {
                        ignore_status(mgr.relate_coincident_pf(pe, pf))?;
                    }
                }
            }
        }

        let mgr = self
            .proto_topo_mgr
            .as_mut()
            .ok_or(undefined("no proto topology manager"))?;
        // 1. split ProtoEdges at ProtoEdge/ProtoEdge XSects
        let mut pes: Vec<ProtoEdgeId> = mgr.proto_edges().to_vec();
        let mut ii = 0;
        while ii < pes.len() {
            let Some(&pe1) = pes.get(ii) else { break };
            let mut jj = ii + 1;
            while jj < pes.len() {
                let Some(&pe2) = pes.get(jj) else { break };
                ignore_status(mgr.split_intersecting_proto_edges(k, geometry, pe1, pe2))?;
                pes = mgr.proto_edges().to_vec();
                jj += 1;
            }
            ii += 1;
        }
        // 2. remove coincident edge-less proto vertices, counting backward
        let pvs = mgr.proto_vertices().to_vec();
        for &pv in pvs.iter().rev() {
            for &pe in &pes {
                let removed = mgr.remove_coincident_no_edge_proto_vertex(k, &*geometry, pe, pv)?;
                if removed {
                    break;
                }
            }
        }
        Ok(())
    }

    /// `SetRelationships`: relates the final vertices of every vertex definition and
    /// the final edges of every edge definition (faces are not related).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] without a proto topology manager;
    /// [`SolidError::Unported`] for `Relate`.
    pub fn set_relationships(&mut self) -> Result<()> {
        let mgr = self
            .proto_topo_mgr
            .as_ref()
            .ok_or(undefined("no proto topology manager"))?;
        for vd in mgr.all_vertex_definitions()? {
            let d = mgr
                .vertex_def(vd)
                .ok_or(undefined("deleted vertex definition"))?;
            match (d.final_vertex(1), d.final_vertex(2)) {
                (Some(t1), Some(t2)) => ignore_status(relate(t1, t2))?,
                _ => continue,
            }
        }
        for ed in mgr.all_edge_definitions()? {
            let Some(d) = mgr.edge_def(ed) else {
                continue;
            };
            match (d.final_edge(1), d.final_edge(2)) {
                (Some(t1), Some(t2)) => ignore_status(relate(t1, t2))?,
                _ => continue,
            }
        }
        Ok(())
    }

    /// `FixupClassMates`: removes short intervals that split mated vertex pairs, then
    /// fixes the end and mid point classifications of a homogenised pair.
    ///
    /// # Errors
    ///
    /// `SM_ERR` for classifications of different sizes; those of
    /// `FixupSplitVertexPair`; [`SolidError::Unported`] for the point fixups.
    pub fn fixup_class_mates<G>(
        &mut self,
        k: &Kernel,
        geometry: &G,
        cc_brep: &mut CurveClassification<O>,
        cc_other: &mut CurveClassification<O>,
    ) -> Result<()>
    where
        G: ClassGeometry<O> + IntersectorTopology<O>,
    {
        let curve_h = cc_brep.curve();
        if curve_h != cc_other.curve() {
            return Ok(());
        }
        let curve = curve_h
            .and_then(|c| geometry.curve(c))
            .ok_or(undefined("classification without curve"))?;
        if cc_brep.size() != cc_other.size() {
            return Err(SolidError::sm_err(OP));
        }
        let mut ii = 0usize;
        let mut restarts = 0usize;
        while ii < cc_brep.size() {
            let ivl_b = *cc_brep.at(ii)?;
            let ivl_o = *cc_other.at(ii)?;
            if !curve.is_closed(&ivl_b.interval, 0.0)? {
                let start_end = ivl_b.start.point_class() == PointClassType::Vertex
                    && ivl_o.end.point_class() == PointClassType::Vertex
                    && self.is_mated_pair(ivl_b.start.object(), ivl_o.end.object());
                let end_start = ivl_b.end.point_class() == PointClassType::Vertex
                    && ivl_o.start.point_class() == PointClassType::Vertex
                    && self.is_mated_pair(ivl_b.end.object(), ivl_o.start.object());
                if start_end || end_start {
                    let changed =
                        self.fixup_split_vertex_pair(k, geometry, ii, cc_brep, cc_other)?;
                    if changed {
                        restarts += 1;
                        if restarts > cc_brep.size() + 1 {
                            return Err(undefined("classification fixup does not converge"));
                        }
                        ii = 0;
                        continue;
                    }
                }
            }
            if ivl_b.start.point_class() != PointClassType::Unknown
                && ivl_o.start.point_class() != PointClassType::Unknown
            {
                fixup_point_class_mates("SmTopologyIntersector::FixupEndPointClassMates")?;
            }
            if ivl_b.mid.point_class() != PointClassType::Unknown
                && ivl_o.mid.point_class() != PointClassType::Unknown
            {
                fixup_point_class_mates("SmTopologyIntersector::FixupMidPointClassMates")?;
            }
            if ii + 1 == cc_brep.size()
                && ivl_b.end.point_class() != PointClassType::Unknown
                && ivl_o.end.point_class() != PointClassType::Unknown
            {
                fixup_point_class_mates("SmTopologyIntersector::FixupEndPointClassMates")?;
            }
            ii += 1;
        }
        Ok(())
    }

    /// `FixupSplitVertexPair`: removes interval `index` when its ends are a mated
    /// vertex pair split by tolerance, moving both vertex classifications to a common
    /// drop parameter on the neighbouring interval ends. Returns whether the
    /// classifications changed.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when neither vertex drops to the curve; [`SolidError::InvalidInput`]
    /// for an index out of range or unresolved curve/vertex handles.
    pub fn fixup_split_vertex_pair<G>(
        &self,
        k: &Kernel,
        geometry: &G,
        index: usize,
        cc_brep: &mut CurveClassification<O>,
        cc_other: &mut CurveClassification<O>,
    ) -> Result<bool>
    where
        G: ClassGeometry<O> + IntersectorTopology<O>,
    {
        let size = cc_brep.size();
        let ivl_b = *cc_brep.at(index)?;
        let ivl_o = *cc_other.at(index)?;
        let brep_src_zone_tol = src_zone_tol_3d_of(k, Some(cc_brep));
        let other_src_zone_tol = src_zone_tol_3d_of(k, Some(cc_other));
        let curve = cc_brep
            .curve()
            .and_then(|c| geometry.curve(c))
            .ok_or(undefined("classification without curve"))?;
        let ivl = cc_brep.interval().union(&cc_other.interval());

        if curve.is_closed(&ivl_b.interval, 0.0)? {
            return Ok(false);
        }

        // sPointClassB/O are constructed with the source zone tolerances, then
        // overwritten by the matched classifications (operator= copies every member).
        let _ = (brep_src_zone_tol, other_src_zone_tol);
        let (mut pc_b, mut pc_o) = if ivl_b.start.point_class() == PointClassType::Vertex
            && ivl_o.end.point_class() == PointClassType::Vertex
            && self.is_mated_pair(ivl_b.start.object(), ivl_o.end.object())
        {
            (ivl_b.start, ivl_o.end)
        } else if ivl_b.end.point_class() == PointClassType::Vertex
            && ivl_o.start.point_class() == PointClassType::Vertex
            && self.is_mated_pair(ivl_b.end.object(), ivl_o.start.object())
        {
            (ivl_b.end, ivl_o.start)
        } else {
            return Ok(false);
        };

        let guess_b = pc_b.t_param();
        let guess_o = pc_o.t_param();
        let vb = pc_b.vertex_object().ok_or(undefined("no vertex object"))?;
        let vo = pc_o.vertex_object().ok_or(undefined("no vertex object"))?;
        // SM_USE_OLDTOL: GetDeviation for the primary, GetGap3d for the other (both the
        // stored gap); statuses ignored, the outputs keep their initial values.
        let drop_b = ignored(curve.drop_point(
            k,
            &ivl,
            geometry.vertex_point(vb)?,
            None,
            10.0 * pc_b.deviation(),
            Some(guess_b),
            SolverOperation::Minimize,
        ))?
        .unwrap_or_else(DroppedPoint::default);
        let drop_o = ignored(curve.drop_point(
            k,
            &ivl,
            geometry.vertex_point(vo)?,
            None,
            10.0 * pc_o.deviation(),
            Some(guess_o),
            SolverOperation::Minimize,
        ))?
        .unwrap_or_else(DroppedPoint::default);

        if !drop_b.success && !drop_o.success {
            return Err(SolidError::sm_err(OP));
        }
        let param = if drop_b.success && drop_o.success && drop_b.distance < drop_o.distance / 5.0 {
            drop_b.parameter
        } else if drop_b.success && drop_o.success {
            (drop_b.parameter + drop_o.parameter) / 2.0
        } else if drop_b.success {
            drop_b.parameter
        } else {
            drop_o.parameter
        };

        pc_b.set_t_param(param, false);
        pc_b.set_deviation(if drop_b.success {
            drop_b.distance
        } else {
            pc_b.deviation()
        });
        pc_o.set_t_param(param, false);
        pc_o.set_deviation(if drop_o.success {
            drop_o.distance
        } else {
            pc_o.deviation()
        });

        if index >= 1 {
            cc_brep.at_mut(index - 1)?.end = pc_b;
            cc_other.at_mut(index - 1)?.end = pc_o;
        }
        if size > 1 && index + 2 <= size {
            cc_brep.at_mut(index + 1)?.start = pc_b;
            cc_other.at_mut(index + 1)?.start = pc_o;
        }
        cc_brep.remove_at(index)?;
        cc_other.remove_at(index)?;
        Ok(true)
    }

    /// `EdgeEdgeIntersect`: intersects the edges of two lists in their interiors,
    /// splitting them at new vertices; split pieces are processed in turn.
    ///
    /// # Errors
    ///
    /// Those of `sm_TestAndSplitTwoEdges` unless `ignore_problems`;
    /// [`SolidError::Unported`] for it.
    pub fn edge_edge_intersect(
        &mut self,
        k: &Kernel,
        geometry: &impl IntersectorTopology<O>,
        input_edges1: &[Option<O>],
        input_edges2: &[Option<O>],
        ignore_problems: bool,
    ) -> Result<EdgeEdgeOutputs<O>> {
        let mut out = EdgeEdgeOutputs {
            split_edges1: Vec::new(),
            new_edges1: Vec::new(),
            split_edges2: Vec::new(),
            new_edges2: Vec::new(),
            new_verts: Vec::new(),
            max_edge_gap: 0.0,
            min_unstitched_edge_gap: BIG_DOUBLE,
        };
        let mut num_edges1 = input_edges1.len();
        if num_edges1 < 1 {
            return Ok(out);
        }
        let mut num_edges2 = input_edges2.len();
        if num_edges2 < 1 {
            return Ok(out);
        }
        let mut edges1: Vec<Option<O>> = input_edges1.to_vec();
        let mut edges2: Vec<Option<O>> = input_edges2.to_vec();

        let mut ii = 0;
        while ii < num_edges1 {
            let Some(Some(e1)) = edges1.get(ii).copied() else {
                ii += 1;
                continue;
            };
            // status ignored: a failed box keeps its initial value
            let e1_box = ignored(geometry.edge_bounding_box(k, e1))?.unwrap_or_default();
            let mut jj = 0;
            while jj < num_edges2 {
                let Some(Some(e2)) = edges2.get(jj).copied() else {
                    jj += 1;
                    continue;
                };
                let e2_box = ignored(geometry.edge_bounding_box(k, e2))?.unwrap_or_default();
                let this_tol = three_max(
                    self.this_approx_tol_3d,
                    geometry.topology_tolerance(e1)?,
                    geometry.topology_tolerance(e2)?,
                );
                if e1_box.are_disjoint(&e2_box, this_tol) {
                    jj += 1;
                    continue;
                }
                // We can also skip edges that share a vertex.
                let v11 = geometry.edge_vertex(e1)?;
                let v12 = geometry.edge_other_vertex(e1, v11)?;
                let v21 = geometry.edge_vertex(e2)?;
                if v21 == v11 || v21 == v12 {
                    jj += 1;
                    continue;
                }
                let v22 = geometry.edge_other_vertex(e2, v21)?;
                if v22 == v11 || v22 == v12 {
                    jj += 1;
                    continue;
                }

                let split = match test_and_split_two_edges(e1, e2, this_tol) {
                    Ok(s) => s,
                    Err(e) if e.is_unported() => return Err(e),
                    Err(e) => {
                        if ignore_problems {
                            jj += 1;
                            continue;
                        }
                        return Err(e);
                    }
                };
                let (dist_to_edge, min_unstitch, new_v, new_e1, new_e2) = split;
                if new_v.is_empty() {
                    jj += 1;
                    continue;
                }
                out.min_unstitched_edge_gap = out.min_unstitched_edge_gap.min(min_unstitch);
                out.max_edge_gap = out.max_edge_gap.max(dist_to_edge);
                out.new_verts.extend_from_slice(&new_v);
                out.new_edges1.extend_from_slice(&new_e1);
                out.new_edges2.extend_from_slice(&new_e2);
                if !out.split_edges1.contains(&e1) {
                    out.split_edges1.push(e1);
                }
                if !out.split_edges2.contains(&e2) {
                    out.split_edges2.push(e2);
                }
                edges1.push(Some(e1));
                edges1.extend(new_e1.iter().map(|&e| Some(e)));
                edges2.push(Some(e2));
                edges2.extend(new_e2.iter().map(|&e| Some(e)));
                if edges1.len() > MAX_EDGE_LIST || edges2.len() > MAX_EDGE_LIST {
                    return Err(undefined("edge list growth past the bound"));
                }
                num_edges1 = edges1.len();
                num_edges2 = edges2.len();
                // Update tolerances.
                update_object_tolerance(e1, dist_to_edge)?;
                update_object_tolerance(e2, dist_to_edge)?;
                for &v in &new_v {
                    update_object_tolerance(v, dist_to_edge)?;
                }
                for &e in new_e1.iter().chain(new_e2.iter()) {
                    update_object_tolerance(e, dist_to_edge)?;
                }
                jj += 1;
            }
            ii += 1;
        }
        Ok(out)
    }

    /// `EdgeFaceCoincidence`: embeds the edges joining two given vertices into the
    /// faces they are coincident with.
    ///
    /// # Errors
    ///
    /// Those of the topology reads; [`SolidError::Unported`] for the face domain,
    /// classification and imprint callees.
    pub fn edge_face_coincidence(
        &mut self,
        k: &Kernel,
        geometry: &impl IntersectorTopology<O>,
        given_verts: &[Option<O>],
        excluded_edges: &[O],
        ignore_problems: bool,
    ) -> Result<EdgeFaceOutputs<O>> {
        let mut out = EdgeFaceOutputs {
            imprinted_edges: Vec::new(),
            processed_faces: Vec::new(),
            new_faces: Vec::new(),
            split_faces: Vec::new(),
            max_edge_gap: 0.0,
            this_min_unstitched: BIG_DOUBLE,
        };
        for &v in given_verts {
            let Some(v) = v else {
                continue;
            };
            let edges_of_vtx = geometry.vertex_edges(v)?;
            for &e in &edges_of_vtx {
                if excluded_edges.contains(&e) {
                    continue;
                }
                let other_v = geometry.edge_other_vertex(e, Some(v))?;
                if !given_verts.contains(&other_v) {
                    continue;
                }
                let other_v = other_v.ok_or(undefined("NULL other vertex"))?;
                let faces1 = geometry.vertex_faces(v)?;
                let faces2 = geometry.vertex_faces(other_v)?;
                let common: Vec<O> = faces1
                    .iter()
                    .copied()
                    .filter(|f| faces2.contains(f))
                    .collect();
                let edge_faces = geometry.edge_faces(e)?;
                for &f in &common {
                    if edge_faces.contains(&f) {
                        continue;
                    }
                    let curve = geometry
                        .edge_curve(e)?
                        .ok_or(undefined("edge without curve"))?;
                    let surf = geometry
                        .face_surface_of(f)?
                        .ok_or(undefined("face without surface"))?;
                    let edge_dom = geometry.edge_interval(e)?;
                    let edge_zone_tol = geometry.edge_zone_tol_3d(k, e)?;
                    let face_dom = face_uv_domain_from_trim_curves(f)?;
                    // status ignored
                    let sol = ignored(surf.simple_coincidence_checker::<NoObject>(
                        k,
                        &face_dom,
                        curve,
                        &edge_dom,
                        self.this_approx_tol_3d,
                        false,
                        [None, None],
                    ))?
                    .flatten();
                    let Some(sol) = sol else {
                        continue;
                    };
                    if sol.start.solution_value > self.this_approx_tol_3d {
                        continue;
                    }
                    // Test two points of the curve against the face.
                    let mut pt_class = PointClassification::new(edge_zone_tol);
                    let crv_pt = curve.evaluate_point(edge_dom.evaluate(0.38))?;
                    ignore_status(face_point_3d_classify(
                        f,
                        crv_pt,
                        edge_zone_tol,
                        false,
                        &mut pt_class,
                    ))?;
                    if !on_face_class(pt_class.point_class()) {
                        continue;
                    }
                    if pt_class.deviation() > self.this_approx_tol_3d {
                        continue;
                    }
                    let crv_pt = curve.evaluate_point(edge_dom.evaluate(0.71))?;
                    let tol = self.this_approx_tol_3d;
                    ignore_status(face_point_3d_classify(f, crv_pt, tol, false, &mut pt_class))?;
                    if !on_face_class(pt_class.point_class()) {
                        continue;
                    }
                    if pt_class.deviation() > self.this_approx_tol_3d {
                        continue;
                    }
                    // The orientation and tolerance go into MakeEdgeInFace.
                    let _ = geometry.vertex_point(v)?;
                    let _ = geometry.topology_tolerance(e)?;
                    let (new_edge, _new_loop, new_face) = brep_make_edge_in_face(f, v, other_v)?;
                    if !out.processed_faces.contains(&f) {
                        out.processed_faces.push(f);
                    }
                    let dev = sol.start.solution_value;
                    if let Some(nf) = new_face {
                        out.new_faces.push(nf);
                        out.processed_faces.push(nf);
                        out.split_faces.push(f);
                        update_object_tolerance(nf, dev)?;
                    }
                    if let Some(ne) = new_edge {
                        update_object_tolerance(ne, dev)?;
                        update_object_tolerance(f, dev)?;
                        match brep_glue_edges_general(e, ne) {
                            Ok(()) => out.imprinted_edges.push(e),
                            Err(err) if err.is_unported() => return Err(err),
                            Err(_) => {
                                if !ignore_problems {
                                    return Err(SolidError::sm_err(OP));
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(out)
    }
}

/// A point class on any part of a face (face, edge or vertex).
fn on_face_class(t: PointClassType) -> bool {
    matches!(
        t,
        PointClassType::Face | PointClassType::Edge | PointClassType::Vertex
    )
}

/// `smos_3Max`.
fn three_max(a: f64, b: f64, c: f64) -> f64 {
    let ab = if a > b { a } else { b };
    if ab > c {
        ab
    } else {
        c
    }
}
