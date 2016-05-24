
.onAttach <- function(library, pkg)
{
buildStartupMessage(pkg = "distrRmetrics", "",  library = library,
                    packageHelp = TRUE, 
                    VIGNETTE = gettext(
"Package \"distrDoc\" provides a vignette to this package as well as to several extension packages; try vignette(\"distr\")."
                                      )
                   )
  invisible()
} 


## Class: SkewNormParameter
setClass("SNormParameter",
          representation = representation(xi="numeric"),
          prototype = prototype(name =
                      gettext("Parameter of a Skewed Normal distribution"),
                       xi=1.5),
          contains = "UniNormParameter"
          )

## Class: STParameter
setClass("SSTdParameter",
          representation = representation(mean = "numeric", sd = "vector", nu="numeric", xi="numeric"),
          prototype = prototype(mean = 0, sd = 1, nu=5, xi=1.5, name =
                      gettext("Parameter of a Skewed-T-distribution")
                      ),
          contains = "UniNormParameter"
          )

## Class: Skewed normal distribution (snorm in fGarch)
setClass("SNorm",
          prototype = prototype(
                      r=function(n)rsnorm(n,mean=0,sd=1,xi=1.5),
                      d=function(x, log=FALSE){
                          d0 <- dsnorm(x,mean=0,sd=1,xi=1.5)
                          return(if(log) log(d0) else d0)
                          },
                      p=function(q, lower.tail=TRUE, log.p=FALSE){
                          p00 <- psnorm(q,mean=0,sd=1,xi=1.5)
                          p0  <- if(lower.tail) p00 else 1-p00
                          return(if(log.p) log(p0) else p0)
                          },
                      q=function(p, lower.tail=TRUE, log.p=FALSE){
                          p00 <- if(log.p) exp(p) else p
                          p0 <- if(lower.tail) p00 else 1-p00
                          return(qsnorm(p0,mean=0,sd=1,xi=1.5))
                          },
                      param = new("SNormParameter"),
                     .withArith = FALSE,
                     .withSim = FALSE,
                     .logExact = FALSE,
                     .lowerExact = FALSE),
          contains = "AbscontDistribution"
          )


## Class: Skewed T distribution (sstd in fGarch)
setClass("SSTd",
          prototype = prototype(
                      r=function(n)rsstd(n,mean=0,sd=1,nu=5,xi=1.5),
                      d=function(x, log=FALSE){
                          d0 <- dsstd(x,mean=0,sd=1,nu=5,xi=1.5)
                          return(if(log) log(d0) else d0)
                          },
                      p=function(q, lower.tail=TRUE, log.p=FALSE){
                          p00 <- psstd(q,mean=0,sd=1,nu=5,xi=1.5)
                          p0  <- if(lower.tail) p00 else 1-p00
                          return(if(log.p) log(p0) else p0)
                          },
                      q=function(p, lower.tail=TRUE, log.p=FALSE){
                          p00 <- if(log.p) exp(p) else p
                          p0 <- if(lower.tail) p00 else 1-p00
                          return(qsstd(p0,mean=0,sd=1,nu=5,xi=1.5))
                          },
                      param = new("SSTdParameter"),
                     .withArith = FALSE,
                     .withSim = FALSE,
                     .logExact = FALSE,
                     .lowerExact = FALSE
                      ),
          contains = "AbscontDistribution"
          )
