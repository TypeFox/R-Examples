ptgpd_log <-
function(x,y,z,mar1=c(0,1,0.1),mar2=c(0,1,0.1),mar3=c(0,1,0.1),dep=1.5)
{
error=FALSE       
hxyz=NULL
                                  #frechet
                                  param   = as.numeric(c(mar1, mar2, mar3, dep))
                                  mux     = param[1]
                                  muy     = param[4]
                                  muz     = param[7]
                                  sigx    = param[2]
                                  sigy    = param[5]
                                  sigz    = param[8]
                                  gamx    = param[3]
                                  gamy    = param[6]
                                  gamz    = param[9]
                                  alpha   = param[10]

                                  if(gamx>0)
                                  {
                                  epx1=mux-sigx/gamx
                                  epx2=Inf
                                  } else {
                                  epx1=-Inf
                                  epx2=mux-sigx/gamx
                                  }

                                  if(gamy>0)
                                  {
                                  epy1=muy-sigy/gamy
                                  epy2=Inf
                                  } else {
                                  epy1=-Inf
                                  epy2=muy-sigy/gamy
                                  }

                                  if(gamz>0)
                                  {
                                  epz1=muz-sigz/gamz
                                  epz2=Inf
                                  } else {
                                  epz1=-Inf
                                  epz2=muz-sigz/gamz
                                  }

                                  #2Dbe is be kene tenni:
                                  if((min(x)<epx1)|(max(x)>epx2)) {error=T}#;print("Invalid parameters for x.")}
                                  if((min(y)<epy1)|(max(y)>epy2)) {error=T}#;print("Invalid parameters for y.")}
                                  if((min(z)<epz1)|(max(z)>epz2)) {error=T}#;print("Invalid parameters for z.")}

    #print(matrix(c(epx1,min(x),max(x),epx2,epy1,min(y),max(y),epy2,epz1,min(z),max(z),epz2),3,4,byrow=T))
    if (sigx < 0 | sigy < 0 | sigz < 0 | alpha < 1) {error=T}#;print("Invalid parameters.")}
    #itt teljesul a konvexitas.
    
    if(!error)
    {
    Hxyz    = NA
    tx      = tr (x, gamx, mux, sigx)
    ty      = tr (y, gamy, muy, sigy)
    tz      = tr (z, gamz, muz, sigz)
    tx0     = tr (0, gamx, mux, sigx)
    tz0     = tr (0, gamy, muy, sigy)
    ty0     = tr (0, gamz, muz, sigz)
    c0      = -mulog(tx0,ty0,tz0,alpha=alpha)
    #print(c0)
    #print(c(dep,alpha))
    ml1=mulog(tx, ty, tz, alpha=alpha)
    ml2=mulog(pmin(tx,rep(tx0,length(tx))), pmin(ty,rep(tx0,length(tx))),pmin(tz,rep(tz0,length(tz))),alpha=alpha)
    Hxyz    = 1/c0*(ml1-ml2)
    } else print("invalid parameter(s)")
Hxyz
}
