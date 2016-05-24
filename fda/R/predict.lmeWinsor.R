predict.lmeWinsor <- function(object, newdata, level=Q, asList=FALSE,
      na.action=na.fail, ...){
##
## 0.  Set up
##
  if(any(grep('lme4', search())>0))
    stop("predict.lmeWinsor requires nlme, which can NOT be used",
         " with lme4 in the search path.") 
  library(nlme)
#  
  cl <- match.call()
# default 'level'  
  Q <- object$dims$Q
##
## 1.  Identify inputs and outputs
##
  mdly <- mdlx <- formula(object) 
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
#
  random <- eval(object$call$random)
  xFixed <- all.vars(mdlx)
#  xNames <- unique(c(all.vars(mdlx), all.vars(random)))
#  xNames <- unique(xFixed, all.vars(random)))
  yName <- all.vars(mdly)
  if(length(yName) != 1)
    stop("'fixed' must include a single column of 'data'")
  if(as.character(mdly[[2]]) != yName)
    stop("lmWinsor can not accept a fixed with a transformed 'y';",
         "  left hand side = ", mdly[[2]], ";  y = ", yName)
##
## 2.  Check limits 
##
  Lower <- object$lower
  got.yL <- (yName %in% names(Lower))
  if(!got.yL)
    stop("No lower limit for ", yName)
#  
  Upper <- object$upper
  got.yU <- (yName %in% names(Upper))
  if(!got.yU)
    stop("No upper limit for ", yName)
##
## 3.  Clip newdata[, xFixed] and predict 
##  
  {
    if(missing(newdata))
      pred <- fitted(object, level, asList)
    else{
      newDat <- newdata
#      
      xNum <- sapply(newDat[xFixed], is.numeric)
      xNms <- xFixed[xNum]
      got.xL <- (xNms %in% names(Lower))
      if(any(!got.xL))
        stop("No lower limit for ",
             paste(names(got.xL[!got.xL]), collapse=", "))
      got.xU <- (xNms %in% names(Upper))
      if(any(!got.xU))
        stop("No upper limit for ",
             paste(names(got.xU[!got.xU]), collapse=", "))
#          
      for(x. in xNms){
        xl <- pmax(newdata[[x.]], Lower[x.])
        newDat[[x.]] <- pmin(xl, Upper[x.])
      }
      cl0 <- as.list(cl)
#     remove reference to 'predict.lmeWinsor'       
      cl0[[1]] <- NULL
#     Prepare to call predict.lme  
      obj <- eval(cl0$object)
      class(obj) <- 'lme'
      cl0$object <- obj
#      
      cl0$newdata <- as.name('newDat')
      pred <- do.call('predict', cl0)      
    }
  }
##
## 4.  Clip pred
##
  yl <- pmax(pred, Lower[yName])
  y. <- pmin(yl, Upper[yName])
##
## 5.  done
##
  y.
}
