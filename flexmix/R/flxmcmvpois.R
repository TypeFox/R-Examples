#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: flxmcmvpois.R 4880 2013-02-10 22:28:57Z gruen $
#

FLXMCmvpois <- function(formula=.~.)
{
  z <- new("FLXMC", weighted=TRUE, formula=formula,
                dist="mvpois", name="model-based Poisson clustering")
  
  z@preproc.y <- function(x){
    storage.mode(x) <- "integer"
    x
  }

  z@defineComponent <- expression({
    logLik <- function(x, y){
      colSums(dpois(t(y), lambda, log=TRUE))
    }

    predict <- function(x, ...){
      matrix(lambda, nrow = nrow(x), ncol=length(lambda),
             byrow=TRUE)
    }
    
    new("FLXcomponent", parameters=list(lambda=lambda), df=df,
        logLik=logLik, predict=predict)
  })
  z@fit <- function(x, y, w, ...){
    with(list(lambda = colSums(w*y)/sum(w), df = ncol(y)),
         eval(z@defineComponent))
  }
  z 
}

