TransformS2Fz <-
function(vs,zs,rp)
{
r=rp-zs

indies=which(r==0)
r[indies]=-1

zf=rp*log(rp/r+0i)
zf[indies]=abs(zf[indies])*Inf  
zf=Re(zf) 

vf=vs/r
vf=vf*rp
vf[indies]=vf[indies]*Inf 

return(list(vf=vf,zf=zf))

}

