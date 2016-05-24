weight.function <-          
function(
phi, 
grouped.fused,
vs,
spl, 
phis.v,
phis.gf,
phis.vs,
phis.sp,
...
)
{
phi.v  <- rep(1, length(phis.v))
phi.gf <- rep(1, length(phis.gf))
phi.vs <- rep(1, length(phis.vs))
phi.sp <- rep(1, length(phis.sp))

phi.v[which(phis.v==1)]  <- phi     # differences
phi.v[which(phis.v==-1)] <- 1-phi   # absolute values

phi.gf[which(phis.gf==1)]  <- grouped.fused     # differences
phi.gf[which(phis.gf==-1)] <- 1-grouped.fused   # grouped

phi.vs[which(phis.vs==1)]  <- vs    # splines
phi.vs[which(phis.vs==-1)] <- 1-vs  # grouped    

phi.sp[which(phis.sp==1)]  <- spl    # splines
phi.sp[which(phis.sp==-1)] <- 1-spl  # grouped

return(phi.v * phi.gf * phi.vs * phi.sp)
}

