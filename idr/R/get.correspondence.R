get.correspondence <-
function(x1, x2, t, spline.df=NULL){
  psi.dpsi <- get.uri.2d(x1, x2, t, t, spline.df)

  psi <- list(t=psi.dpsi$tv[,1], value=psi.dpsi$uri, 
              smoothed.line=psi.dpsi$uri.spl, ntotal=psi.dpsi$ntotal, 
              jump.point=psi.dpsi$jump.left)

  dpsi <- list(t=psi.dpsi$t.binned, value=psi.dpsi$uri.slope,
               smoothed.line=psi.dpsi$uri.der, ntotal=psi.dpsi$ntotal, 
               jump.point=psi.dpsi$jump.left)

  psi.n <- list(t=psi$t*psi$ntotal, value=psi$value*psi$ntotal, 
                 smoothed.line=list(x=psi$smoothed.line$x*psi$ntotal, 
                                    y=psi$smoothed.line$y*psi$ntotal),
                 ntotal=psi$ntotal, jump.point=psi$jump.point) 

  dpsi.n <- list(t=dpsi$t*dpsi$ntotal, value=dpsi$value, 
                 smoothed.line=list(x=dpsi$smoothed.line$x*dpsi$ntotal, 
                                    y=dpsi$smoothed.line$y),
                 ntotal=dpsi$ntotal, jump.point=dpsi$jump.point) 

  return(list(psi=psi, dpsi=dpsi, psi.n=psi.n, dpsi.n=dpsi.n))
}

