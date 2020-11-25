############## Gridap Setup ##################
model = GmshDiscreteModel("geometry.msh")

order = 1
diritags = ["BottomEdge","TopEdge", "BottomNodes", "TopNodes"]
neumanntags = ["LeftEdge","RightEdge","LeftNodes","RightNodes","LeftDEdge","RightDEdge","LeftDNodes","RightDNodes"]
sourcetags = ["SourceNodes","SourceEdge"]
designboundarytags = ["LeftDEdge","RightDEdge","TBDEdge","TBDNodes","LeftDNodes","RightDNodes"]
V = TestFESpace(
  reffe=:Lagrangian,
  order=1,
  valuetype=Float64,
  model=model,
  conformity=:H1,
  dirichlet_tags=diritags)

U = TrialFESpace(V,[0,0,0,0])

degree = 2
trian = Triangulation(model)
quad = CellQuadrature(trian,degree)

strian = BoundaryTriangulation(model,sourcetags)
squad = CellQuadrature(strian,degree)

btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,degree)
############### Get the design region and target region ################
labels = get_face_labeling(model)
dimension = num_cell_dims(model)
tags = get_face_tag(labels,dimension)
const design_tag = get_tag_from_name(labels,"Design")
const target1_tag = get_tag_from_name(labels,"Target1")
const target2_tag = get_tag_from_name(labels,"Target2")
cellmask_d = get_face_mask(labels,"Design",dimension)
cellmask_t1 = get_face_mask(labels,"Target1",dimension)
cellmask_t2 = get_face_mask(labels,"Target2",dimension)
# Subset of the mesh
trian_d = Triangulation(model,cellmask_d)
quad_d = CellQuadrature(trian_d,degree)
trian_t1 = Triangulation(model,cellmask_t1)
quad_t1 = CellQuadrature(trian_t1,degree)
trian_t2 = Triangulation(model,cellmask_t2)
quad_t2 = CellQuadrature(trian_t2,degree)
# Number of cells in design region (number of design parameters)
np = num_cells(trian_d)

# Boundary of the design region
dtrian = BoundaryTriangulation(model,designboundarytags)
dquad = CellQuadrature(dtrian,degree)
######################## Design Parameter FE Space #####################
# Piece-wise constant FE space for the design parameters p (before filter and threshold)
Q = TestFESpace(
  reffe=:Lagrangian,
  order=0,
  valuetype=Float64,
  conformity=:L2,
  model=model)

P = TrialFESpace(Q)

# FE function space for the filtered parameters pf

Qf = TestFESpace(
  model=model,
  reffe=:Lagrangian,
  order=1,
  conformity=:H1,
  valuetype=Float64)
Pf = TrialFESpace(Qf)


################### Assemble matrix and vector #####################
# Material distribution
# pt = Threshold(Filter(p))
@law ξ(pt,α) = ϵ1 + (ϵ2-ϵ1)*pt-α*1im*pt*(1-pt)
# Weak form of the Helmholtz equation : a(p,u,v)=Λ⋅∇v⋅Λ⋅∇u-k²μξ(p)v⋅u
a(u,v,pfh,α,β,η,k,σ,LH,d_pml,flag) = ((x->Λ(x,σ,k,LH,d_pml))⋅∇(v))⊙((x->Λ(x,σ,k,LH,d_pml))⋅∇(u)) - k^2*μ*ξ(Threshold(β,η,pfh,flag),α)*(v⊙u)
# Zero Neumann boundary condition for left and right sides
h(x) = 0
b_Γ(v) = v*h

# Source terms
fgauss(x,I,δ,pos) = I/(sqrt(2π)*δ)*exp(-((x[2]-pos)^2)/δ^2/2)
b_Ω(v) = v*(x->fgauss(x,I,δ,pos))
f(x) = I
b_Γs(v) = v*f



# Construct the finite element matrix and vector in Gridap
function GridapFEM(pf,α,β,η,k,σ,LH,dpml,trian,quad,strian,squad,btrian,bquad,linesource::Bool,flag_t::Bool,flag_f::Bool)
    if (flag_f)  
      pfh  = FEFunction(Pf,pf)
    else
      pfh = FEFunction(P,pf)
    end
    #t_Ω = AffineFETerm((u,v)->a(u,v,pfh,α,β,η,k,σ,LH,dpml),b_Ω,trian,quad)
    if (linesource)
        t_Ω = LinearFETerm((u,v)->a(u,v,pfh,α,β,η,k,σ,LH,dpml,flag_t),trian,quad)
        t_Γ = FESource(b_Γ,btrian,bquad)
        t_Γs = FESource(b_Γs,strian,squad)
        op = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},U,V,t_Ω,t_Γ,t_Γs)
    else
        t_Ω = AffineFETerm((u,v)->a(u,v,pfh,α,β,η,k,σ,LH,dpml,flag_t),b_Ω,trian,quad)
        t_Γ = FESource(b_Γ,btrian,bquad)
        op = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},U,V,t_Ω,t_Γ)
    end
    A = get_matrix(op)
    b = get_vector(op)
    A, b
end