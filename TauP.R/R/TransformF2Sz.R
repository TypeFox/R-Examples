TransformF2Sz <-
function(vf,zf,rp)
{
r=exp(-1*zf/rp) 


vs=r*vf 

zs=rp*(1-r)

return(list(vs=vs,zs=zs))

 
}

