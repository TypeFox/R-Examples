
summary.relative.effect <- function(object,
                                    ...){

  results <- matrix(c(as.numeric(object$adj.treat.cov),
                      as.numeric(object$rel.eff.treat)),
                    
                    nrow=length(object$name.sel),
                    ncol=2,
                    dimnames=list(object$name.sel,
                      c("adjusted treatment effect",
                        "relative effect")))

  res <- list(resp   = object$name.resp,
              treat  = object$name.treat,
              sel    = object$name.sel,
              unadj  = object$unadj.treat,
              result = results)
  

  class(res) <- "summary.relative.effect"

  res

}



