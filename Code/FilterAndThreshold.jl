################### Filter and Threshold #####################
# pf = Filter(p)
a_f(r,u,v) = r^2*(∇(v)⊙∇(u))+v⊙u
function Filter(p,r,P,Pf,Qf,tags,design_tag,trian,quad,dtrian,dquad,flag_f::Bool)
    pvec = p_vec(p,P,tags,design_tag)
    if (flag_f)
        ph = FEFunction(P,pvec)
        t_Ω = AffineFETerm((u,v)->a_f(r,u,v),v->(v*ph),trian,quad)
        t_Γ = FESource(v->(0.0*v),dtrian,dquad)
        op = AffineFEOperator(Pf,Qf,t_Ω,t_Γ)
        A = get_matrix(op)
        b = get_vector(op)
        pf = A\b
        return pf
    else
        return pvec
    end
end

# Threshold function
@law Threshold(β,η,pf,flag_t::Bool) = (!flag_t)*pf+flag_t*((tanh(β*η)+operate(tanh,β*(pf-η)))/(tanh(β*η)+tanh(β*(1.0-η))))