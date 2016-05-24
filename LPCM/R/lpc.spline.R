

lpc.spline <- function(lpcobject, optimize=TRUE, compute.Rc=FALSE,  project=FALSE, ...){
  if (class(lpcobject)=="lpc.spline"){
     lpcobject <-lpcobject$lpcobject
  }   
  lpcsl    <- suppressWarnings(lpc.splinefun(lpcobject))
  spline   <- lpc.fit.spline(lpcsl,...)  
  if (project || compute.Rc){
       proj           <- lpc.project.spline(lpcsl, newdata=lpcobject$data, optimize=optimize,...)
       closest.coords <- proj$closest.coords
  } else {
       closest.coords<- "none"
     }   
  if (project){
       closest.pi     <- round(proj$closest.pi, digits=6)
       closest.dist   <- proj$closest.dist
       closest.branch <- proj$closest.branch
  } else {
      closest.pi <- closest.dist <-  closest.branch <-"none"
  }   
  if (compute.Rc){
       R <- Rc(lpcobject$data, closest.coords)
  } else {
       R<-"none"
  }

  fit <- list(knots.pi=spline$knots.pi, knots.coords= spline$knots.coords, closest.pi=closest.pi, closest.coords=closest.coords, closest.dist=closest.dist,closest.branch=closest.branch, Rc=R, project=project, lpcobject=lpcobject, splinefun=lpcsl)
  class(fit)<-"lpc.spline"
  return(fit)
}


  
  
 

