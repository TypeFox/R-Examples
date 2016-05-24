# Non-linear effect
# tp-spline : if knot < 0 , - (knot - x)+^d  (negative basis if x <= knot

nl <-function(x,
               Spline=c("b-spline","tp-spline", "tpi-spline"), 
               Knots=NULL, 
               Degree=3, 
               Intercept=FALSE,
               Boundary.knots = range(x),
               Keep.duplicates = TRUE,
               outer.ok=TRUE,
               ...){
  
  Spline <- match.arg(Spline)
  if (Spline=="b-spline") {
    xspline  <- MSplineBasis(knots=c(Boundary.knots[1], Knots, Boundary.knots[2]),
                             degree=Degree,
                             keep.duplicates=Keep.duplicates,
                             log=FALSE)

  }
  else if (Spline=="tp-spline") {
    xspline  <- TPSplineBasis(knots=Knots,
                              degree=Degree,  
                              min=Boundary.knots[1],
                              max=Boundary.knots[2],
                              log=FALSE,
                              type="standard")
  }
  else if (Spline=="tpi-spline") {
    xspline  <- TPSplineBasis(knots=Knots,
                              degree=Degree,  
                              min=Boundary.knots[1],
                              max=Boundary.knots[2],
                              log=FALSE,
                              type="increasing")
  }
evaluate(xspline, x, intercept=Intercept, outer.ok=outer.ok)
  
}


# Non-linear effect
# tp-spline : if knot < 0 , - (knot - x)+^d  (negative basis if x <= knot
# same as nl() but multiply by BETAt

nlbeta <-function(y, x,
                   Spline=c("b-spline","tp-spline", "tpi-spline"), 
                   Knots=NULL, 
                   Degree=3, 
                   Intercept=FALSE,
                   Boundary.knots = range(x),
                  Keep.duplicates = TRUE,
                   outer.ok=TRUE,
                   ...){
  
  if (length(x)!=length(y)) {
    stop("x and y have different length.")
  }
    
  nl(x,
      Spline=Spline,
      Knots=Knots,
      Degree=Degree,
      Intercept=Intercept,
      Boundary.knots=Boundary.knots,
      Keep.duplicates = Keep.duplicates,
      outer.ok=outer.ok, ...)*y
}


# Time dependent effect

td <-function(x,timevar,
               Spline=c("b-spline","tp-spline"), 
               Knots.t=NULL, 
               Degree.t=3, 
               Intercept.t=TRUE, 
               Boundary.knots.t = range(timevar), 
               Keep.duplicates.t = TRUE,
               outer.ok=TRUE,
               ...){
    
  Spline <- match.arg(Spline)
  if (Spline=="b-spline") {
    tspline  <- MSplineBasis(knots=c(Boundary.knots.t[1], Knots.t, Boundary.knots.t[2]),
                             degree=Degree.t,
                             keep.duplicates=Keep.duplicates.t,
                             log=FALSE)

  }
  else if (Spline=="tp-spline") {
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                              log=FALSE,
                              type="standard")
  }
  else if (Spline=="tpi-spline") {
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                            log=FALSE,
                              type="standard")
  }
  
  evaluate(tspline, timevar, intercept=Intercept.t, outer.ok=outer.ok)*x
  
}


# Time dependent effect
# same as td but do not multiply by x (return only t basis)

tdalpha <-function(x,timevar,
                   Spline=c("b-spline","tp-spline", "tpi-spline"), 
                   Knots.t=NULL, 
                   Degree.t=3, 
                   Intercept.t=TRUE, 
                   Boundary.knots.t = range(timevar),
                   Keep.duplicates.t = TRUE,
                   outer.ok=TRUE,
                   ...){
  
  Spline <- match.arg(Spline)

  if (Spline=="b-spline") {
    tspline  <- MSplineBasis(knots=c(Boundary.knots.t[1], Knots.t, Boundary.knots.t[2]),
                             degree=Degree.t,
                             keep.duplicates=Keep.duplicates.t,
                             log=FALSE)

  }
  else if (Spline=="tp-spline") {
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                              log=FALSE,
                              type="standard")
  }
 else if (Spline=="tpi-spline") {
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                              log=FALSE,
                              type="standard")
  }
  
  evaluate(tspline, timevar, intercept=Intercept.t, outer.ok=outer.ok)
}




# non-linear and time-dependent effect
# tp-spline : if knot < 0 , - (knot - x)+^d  (negative basis if x <= knot

nltd <- function(x,timevar,
                   model=c("additive","multiplicative"),
                   Spline=c("b-spline","tp-spline", "tpi-spline"), 
                   Knots=NULL, Degree=3,
                   Intercept=FALSE, 
                   Boundary.knots = range(x), 
                   Knots.t=NULL, Degree.t=3,
                   Intercept.t=(model=="multiplicative"), 
                   Boundary.knots.t = range(timevar),
                   outer.ok=TRUE,
                   Keep.duplicates = TRUE,
                   xdimnames=":XxXxXXxXxX ",
                   tdimnames=":TtTtTTtTtT ", ...) {
  
  
  Spline <- match.arg(Spline)
  model <- match.arg(model)
  
  if (Spline=="b-spline") {
    xspline  <- MSplineBasis(knots=c(Boundary.knots[1], Knots, Boundary.knots[2]),
                             degree=Degree,
                             keep.duplicates=Keep.duplicates,
                             log=FALSE)

    tspline  <- MSplineBasis(knots=c(Boundary.knots.t[1], Knots.t, Boundary.knots.t[2]),
                             degree=Degree.t,
                             keep.duplicates=Keep.duplicates,
                             log=FALSE)

  }
  else if (Spline=="tp-spline") {
    xspline  <- TPSplineBasis(knots=Knots,
                              degree=Degree,  
                              min=Boundary.knots[1],
                              max=Boundary.knots[2],
                              log=FALSE,
                              type="increasing")
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                              log=FALSE)
  }
  else if (Spline=="tpi-spline") {
    xspline  <- TPSplineBasis(knots=Knots,
                              degree=Degree,  
                              min=Boundary.knots[1],
                              max=Boundary.knots[2],
                              log=FALSE,
                              type="increasing")
    tspline  <- TPSplineBasis(knots=Knots.t,
                              degree=Degree.t,  
                              min=Boundary.knots.t[1],
                              max=Boundary.knots.t[2],
                              log=FALSE,
                              type="standard")
  }

  xx <- fevaluate(xspline, x, intercept=Intercept, outer.ok=outer.ok)
  tt <- fevaluate(tspline, timevar, intercept=Intercept.t, outer.ok=outer.ok)
  
  if(model == "additive"){
    zz <- cbind(xx, tt*x)
  }
  else {
    zz <- cbind(xx, tt)
  }

  dimnames(zz)[[2]]<-c(paste(xdimnames ,1:(Degree   + Intercept   + length(Knots)) , sep=""),
                      paste(tdimnames ,1:(Degree.t + Intercept.t + length(Knots.t)) , sep=""))

  zz
}






