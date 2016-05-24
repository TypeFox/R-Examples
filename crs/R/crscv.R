crscv <- function(K,
                  I,
                  basis,
                  basis.vec,
                  degree.max, 
                  segments.max, 
                  degree.min, 
                  segments.min, 
                  complexity,
                  knots,
                  degree,
                  segments,
                  restarts,
                  K.mat,
                  lambda,
                  lambda.mat,
                  cv.objc,
                  cv.objc.vec,
                  num.x,
                  cv.func,
                  tau) {

    tregcv = list(K=K,
                  I=I,
                  basis=basis,
                  basis.vec=basis.vec,    
                  degree.max=degree.max, 
                  segments.max=segments.max, 
                  degree.min=degree.min, 
                  segments.min=segments.min, 
                  complexity=complexity,
                  knots=knots,
                  degree=degree,
                  segments=segments,
                  restarts=restarts,
                  K.mat=K.mat,
                  lambda=lambda,
                  lambda.mat=lambda.mat,
                  cv.objc=cv.objc,
                  cv.objc.vec=cv.objc.vec,
                  num.x=num.x,
                  cv.func=cv.func,
                  tau=tau)

  class(tregcv) <- "crscv"

  tregcv
}

print.crscv <- function(x, ...){

  if(!is.null(x$lambda)&&is.null(x$I)) {
    cat("\nCategorical Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective function: ", format(x$cv.func), sep=""))        
    cat(paste("\nObjective function value: ",format(x$cv.objc),sep=""),sep="")

    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))    
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in 1:length(x$degree))
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")
    if(!is.null(x$I)) for(j in 1:length(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in 1:length(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(x$restarts>0) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")    
    cat("\n\n")
  } else if(!is.null(x$I)) {
    cat("\nFactor Regression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective function: ", format(x$cv.func), sep=""))        
    cat(paste("\nObjective function value: ",format(x$cv.objc),sep=""),sep="")
    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))    
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in 1:length(x$degree))
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")
    if(!is.null(x$I)) for(j in 1:length(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in 1:length(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(!is.null(x$restarts) && (x$restarts > 0)) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")    
    cat("\n\n")
  } else {
    cat("\nRegression Spline Cross-Validation",sep="")
    cat(paste("\n\nObjective Function Value: ",format(x$cv.objc),sep=""),sep="")

    cat(paste("\n\nKnot type: ", format(x$knots), sep=""))    
    cat(paste("\nModel complexity proxy: ", format(x$complexity), sep=""))

    for(j in 1:x$num.x)
      cat(paste("\nSpline degree/number of segments for x[", j, "]: ", format(x$degree[j]),"/",format(x$segments[j]),sep=""),sep="")

    if(!is.null(x$I)) for(j in 1:length(x$I))
      cat(paste("\nInclusion indicator for z[", j, "]: ",format(x$I[j]),sep=""),sep="")
    if(!is.null(x$lambda)) for(j in 1:length(x$lambda))
      cat(paste("\nBandwidth for  z[", j, "]: ",format(x$lambda[j]),sep=""),sep="")

    cat(paste("\n\nMaximum spline degree for search: ",format(x$degree.max),sep=""),sep="")
    cat(paste("\nBasis: ", x$basis,sep=""))
    if(!is.null(x$restarts) && (x$restarts > 0)) cat(paste("\nNumber of restarts = ", format(x$restarts),sep=""),sep="")    
    cat("\n\n")
  }
}

summary.crscv <- function(object, ...){
  print(object)
}
