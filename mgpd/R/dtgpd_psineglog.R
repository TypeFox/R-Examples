dtgpd_psineglog <-
function(x, y, z, mar1=c(0,1,0.1), mar2=c(0,1,0.1), mar3=c(0,1,0.1), dep=1.5, A1=0, A2=0, B1=3, B2=3, checkconv=TRUE, ...)
{
error=FALSE
hxyz=NULL
                                  #frechet
                                  param   = as.numeric(c(mar1, mar2, mar3, dep, A1, A2, B1, B2))
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
                                  
if(checkconv)
          {
          conv          = function(x)
          {
          if(sum(x)<1)
                    {
                    conv.mtx      = matrix(NA,2,2)
                    conv.mtx[1,1] = d11psiAneglog(x[1],x[2], alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
                    conv.mtx[1,2] = d12psiAneglog(x[1],x[2], alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
                    conv.mtx[2,1] = d12psiAneglog(x[1],x[2], alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
                    conv.mtx[2,2] = d22psiAneglog(x[1],x[2], alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
                    #if(!is.nan(conv.mtx[1,1])) all( eigen(conv.mtx)$values >0 ) else NA
                    is.positive.definite(conv.mtx, method=c("chol"))
                    } else NA
          }
          
          conv.outer    = function(x,y) apply(cbind(x,y),1,conv)
          x1            = seq(0.01,1-0.01,0.01)
          y1            = seq(0.01,1-0.01,0.01)
          Conv          = outer(x1,y1,conv.outer)
          if (min(Conv,na.rm=T)==0) 
          {
          error=T
          image.plot(x1,y1,Conv,col=heat.colors(2))
          }
}


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

                                  if((min(x)<epx1)|(max(x)>epx2)) {error=T}#;print("Invalid parameters for x.")}
                                  if((min(y)<epy1)|(max(y)>epy2)) {error=T}#;print("Invalid parameters for y.")}
                                  if((min(z)<epz1)|(max(z)>epz2)) {error=T}#;print("Invalid parameters for z.")}

    #print(matrix(c(epx1,min(x),max(x),epx2,epy1,min(y),max(y),epy2,epz1,min(z),max(z),epz2),3,4,byrow=T))

    if(!error)
    {
    tx      = tr (x, gamx, mux, sigx)
    ty      = tr (y, gamy, muy, sigy)
    tz      = tr (z, gamz, muz, sigz)
    tx0     = tr (0, gamx, mux, sigx)
    tz0     = tr (0, gamy, muy, sigy)
    ty0     = tr (0, gamz, muz, sigz)
    dtx     = dtr(x, gamx, mux, sigx)
    dty     = dtr(y, gamy, muy, sigy)
    dtz     = dtr(z, gamz, muz, sigz)
    c0       = -1/mupsineglog (tx0,ty0,tz0, alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
    dddpsimu = d123mupsineglog(tx ,ty ,tz , alpha=alpha, A1=A1, A2=A2, B1=B1, B2=B2)
    Jc       = dtx*dty*dtz
    #c0
    #summary(Jc)
    #summary(dddpsimu)
    null     = (1 - ((tx < tx0) * (ty < ty0) * (tz < tz0)))
    hxyz     = -c0*dddpsimu*null*Jc
    } else print("invalid parameter(s)")
hxyz
}
