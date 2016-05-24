#
# transformations.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2012
# Contains: transform, donothing, msqrt,mlog, reciproce, mprobit, mlogit
#

donothing <- function(x, ...){ invisible(return(x)) }
msqrt     <- function(x, ...){ invisible(return(sqrt(x))) }
mlog      <- function(x, ...){ invisible(return(log(x, ...))) }
reciproce <- function(x, ...){ invisible(return(1/x)) }
mprobit   <- function(x, ...){ invisible(return(probit(x, ...))) }
mlogit    <- function(x, ...){ invisible(return(logit(x, ...))) }

transformation <- function(x, transformations=c("nothing","log","sqrt","reciprocal","probit","logit"), ... , verbose=TRUE){
  optionsAvailable <- c("nothing","log","sqrt","reciprocal","probit","logit")
  chosen  <- pmatch(transformations, optionsAvailable)
  methodsChosen <- c(donothing, mlog, msqrt, reciproce, mprobit, mlogit)

  res <- vector("list",length(methodsChosen))
  idx <- 1
  for(n in chosen){
    if(verbose)cat("Applying",n,"a",optionsAvailable[n],"transformation to the data\n")
    res[[idx]] <- apply(x,1,function(t){
      (methodsChosen[n][[1]])(t, ...)
    })
  idx <- idx+1
  }
  res
}

