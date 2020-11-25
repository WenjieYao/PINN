# PML coordinate streching functions
function s_PML(x,σ,k,LH,d_pml)
    L = LH[1]
    h1 = LH[2]
    h2 = LH[3]
    ux = abs(x[1])-L/2.0
    uy = (x[2]>h1)*(x[2]-h1) + (x[2]<-h2)*(-x[2]-h2)
    u = [ux,uy]
    return @. ifelse(u > 0,  1-(1im*σ/k)*(u/d_pml)^2, $(1.0+0im))
end

function Λ(x,σ,k,LH,d_pml)  
    s_x,s_y = s_PML(x,σ,k,LH,d_pml)
    return TensorValue(1/s_x,0,0,1/s_y) 
end