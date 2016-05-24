setClass("FLXMCfactanal",
         contains = "FLXMC")

###**********************************************************

FLXMCfactanal <- function(formula=.~., factors = 1, ...)
{
    z <- new("FLXMCfactanal", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="mixtures of factor analyzers")
    
    z@fit <- function(x, y, w, ...){
      cov.weighted <- cov.wt(y, wt = w)[c("center","cov")]
      cov <- cov.weighted$cov; center <- cov.weighted$center
      fa <- factanal(covmat = cov, factors = factors, ...)
      Sigma <- fa$loadings %*% t(fa$loadings) + diag(fa$uniquenesses)
      df <- (factors + 2) * ncol(y)
        
      predict <- function(x)
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      
      logLik <- function(x, y){
        sds <- sqrt(diag(cov))
        mvtnorm::dmvnorm(y, mean = center, sigma = Sigma * (sds %o% sds), log = TRUE)
      }
      
      new("FLXcomponent", parameters=list(mu = center,
                            variance = diag(cov),
                            loadings = fa$loadings,
                            uniquenesses = fa$uniquenesses),
          df=df, logLik=logLik, predict=predict)
    }
    z
}


###**********************************************************

setMethod("rFLXM", signature(model = "FLXMCfactanal", components = "FLXcomponent"),
          function(model, components, class, ...) {
            FUN <- paste("r", model@dist, sep = "")
            Sigma <- components@parameters$loadings %*% t(components@parameters$loadings) + diag(components@parameters$uniquenesses)
            sds <- sqrt(components@parameters$variance)
            args <- list(n = nrow(model@x), mean = components@parameters$mu,
                      sigma =  Sigma * (sds %o% sds))
            return(do.call(FUN, args))
          })
