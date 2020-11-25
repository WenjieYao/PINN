@law myabs2(x) = abs2(x)

#g=g_u(uvec)
function g_u(uvec)
    uh_temp = FEFunction(U,uvec)
    uh_t2 = restrict(uh_temp,trian_t2)
    uh_t1 = restrict(uh_temp,trian_t1)
    g2=sum(integrate(myabs2(myabs2(uh_t2)-It^2),trian_t2,quad_t2))
    g1=sum(integrate(myabs2(myabs2(uh_t1)),trian_t1,quad_t1))
    #g3=sum(integrate(myabs2(myabs2(uh_temp)-f2_target),trian,quad))
    g1+g2
end
#uvec = u_pf(pf)
function u_pf(pf)
    A,b = GridapFEM(pf,α,β,η,k,σ,LH,dpml,trian,quad,strian,squad,btrian,bquad,linesource,flag_t,flag_f)
    uvec = A\b
    uvec
end
#pf = pf_p(p)
function pf_p(p)
    pf = Filter(p,r,P,Pf,Qf,tags,design_tag,trian,quad,dtrian,dquad,flag_f)
    pf
end
# Chain Rule : dg/dp = dg/dg*dg/du*du/dpf*dpf/dp
# dg/du=dg/dg*dg/du
function rrule(::typeof(g_u),uvec)
  function g_pullback(dgdg)
    NO_FIELDS, dgdg*Dgdu(uvec)
  end
  g_u(uvec), g_pullback
end

function Dgdu(uvec)
  uh = FEFunction(U,uvec)
  uh_t1 = restrict(uh,trian_t1)
  uh_t2 = restrict(uh,trian_t2)
  t1 = FESource(du->(4*uh_t1*myabs2(uh_t1)*du),trian_t1,quad_t1)
  op1 = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},U,V,t1)
  dgdu1 = get_vector(op1)
  t2 = FESource(du->(4*uh_t2*(myabs2(uh_t2)-It^2)*du),trian_t2,quad_t2)
  op2 = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},U,V,t2)
  dgdu2 = get_vector(op2)
  return dgdu1+dgdu2
end

# dg/dpf=dg/du*du/dpf
function rrule(::typeof(u_pf),pf)
  uvec = u_pf(pf)
  function u_pullback(dgdu)
    NO_FIELDS, Dgdpf(dgdu,uvec,pf)
  end
  uvec, u_pullback
end

@law Dxidpf(pfh)=(ϵ2-ϵ1-α*1im*(1-2*Threshold(β,η,pfh,flag_t)))*(!flag_t+flag_t*β*(1.0-operate(tanh,β*(pfh-η))^2)/(tanh(β*η)+tanh(β*(1.0-η))))

function Dgdpf(dgdu,uvec,pf)

  dG(pfh,u,v,dp) = Dxidpf(pfh)*k^2*μ*(v⊙u)*dp

  A,b = GridapFEM(pf,α,β,η,k,σ,LH,dpml,trian,quad,strian,squad,btrian,bquad,linesource,flag_t,flag_f)
  λvec = A'\dgdu
  
  uh = FEFunction(U,uvec)
  λh = FEFunction(V,conj(λvec))

  
  if (flag_f)
    ph = FEFunction(Pf,pf)
    t = FESource((dp)->dG(ph,uh,λh,dp),trian,quad)
    op = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},Pf,Qf,t)
  else
    ph = FEFunction(P,pf)
    t = FESource((dp)->dG(ph,uh,λh,dp),trian,quad)
    op = AffineFEOperator(SparseMatrixCSC{ComplexF64,Int},Vector{ComplexF64},P,Q,t)
  end
  real(get_vector(op))
end

# dg/dp=dg/dpf*dpf/dp
function rrule(::typeof(pf_p),p)
  function pf_pullback(dgdpf)
    NO_FIELDS, Dgdp(dgdpf,flag_f)
  end
  Filter(p,r,P,Pf,Qf,tags,design_tag,trian,quad,dtrian,dquad,flag_f), pf_pullback
end

function Dgdp(dgdpf,flag_f::Bool)
  if (flag_f)
    t_Ω = LinearFETerm((u,v)->a_f(r,u,v),trian,quad)
    op = AffineFEOperator(Pf,Qf,t_Ω)
    A = get_matrix(op)
    λvec = A'\dgdpf

    λh = FEFunction(Pf,λvec)
    t = FESource(dp->(λh*dp),trian,quad)
    op2 = AffineFEOperator(P,Q,t)
    return extract_design(get_vector(op2),np,tags,design_tag)
  else
    return extract_design(dgdpf,np,tags,design_tag)
  end
end


# Final objective function
function g_p(p::Vector)
    #g(u_pt(Threshold(Filter(p))))
    pf = pf_p(p)
    uvec = u_pf(pf)
    g_u(uvec)
end

function g_p(p::Vector,grad::Vector)
    if length(grad) > 0
        dgdp, = Zygote.gradient(g_p,p)
        grad[:] = dgdp
    end
    @show g_value = g_p(p)
    return g_value
end