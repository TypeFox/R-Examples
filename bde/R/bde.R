#' Main user function to construct bounded densities
#'
#' This function allows simple access to all the bounded density estimators
#' implemented in the package
#'
#' @param dataPoints vector containing the data points from where the density will be estimated
#' @param dataPointsCache list of points where the density will be estimated. By default, 200 points are created uniformely in the [0,1] interval
#' @param estimator string indicating the estimator to be used. It can be one of the list "chen99" (see ...), "microbeta" (see ...), 
#' "macrobeta" (see ...), "vitale" (see ...), "bkernel" (see ...), "hirukawa" (see ...)
#' @param b bandwidth of the estimator. By default it is set at the length of the sample to -2/5
#' @param lower.limit lower bound of the interval where the function should be defined. By default it is 0, unles the smallest sample value
#' is negative; in that case the lower limit is set to 0.99 times the smallest value
#' @param upper.limit upper bound of the interval where the function should be defined. By default it is 1, unles the highest sample value
#' is greater than 1; in that case the upper limit is set to 1.01 times the highest value
#' @param options list of options for the estimator. The particular options will depend on the estimator used
#' @keywords density estimator
#' @export
#' @examples 

bde<-function (dataPoints,dataPointsCache=NULL,estimator,b=length(sample)^{-2/5},
               lower.limit=0, upper.limit=1,options=NULL){
  sample <- dataPoints
  x <- dataPointsCache
  if (min(sample)<min(0,lower.limit)) lower.limit<-0.99*min(sample)
  if (max(sample)>max(1,upper.limit)) upper.limit<-1.01*max(sample)
  if (is.null(x)) x<-seq(lower.limit,upper.limit,(upper.limit - lower.limit)/100)
  density<-switch(estimator,
                  "betakernel"= beta_aux(sample,x,b,options,lower.limit,upper.limit),
                  "vitale" = vitale_aux(sample,x,b,options,lower.limit,upper.limit),
                  "boundarykernel" = boundary_aux(sample,x,b,options,lower.limit,upper.limit),
                  "kakizawa" = kakizawa_aux(sample,x,b,options,lower.limit,upper.limit),
                  {
                    stop("The estimator is a mandatory parameter that has to be either 'betakernel', 'vitale', 'boundarykernel' or 'kakizawa'. Note that the parameter is case sensitive. For more information please type ?bde")
                  })
  density
}


# Default values ----------------------------------------------------------

default<-list()
default$beta<-list()
default$beta$mod=F
default$beta$normalization='none'
default$beta$mbc='none'
default$beta$c=0.5
default$vitale<-list()
default$vitale$biasreduction=F
default$vitale$M=1
default$boundary$mu<-1
default$boundary$kernel<-"muller94"
default$boundary$correct<-F


# Auxiliar functions ------------------------------------------------------

beta_aux<-function(sample,x,b,options,lower.limit,upper.limit){
  if (mode(sample)=="list") sample<-unlist(sample)
  sample<-as.numeric(sample)
  mod=default$beta$mod
  if (!is.null(options$modified)){mod=options$modified}
  if (!is.logical(mod)) stop("The 'modified' parameter used has to be a logical value. For more information please type ?bde")
  
  normalization=default$beta$normalization
  if (!is.null(options$normalization)) {normalization=options$normalization}
  
  mbc=default$beta$mbc
  if (!is.null(options$mbc)) mbc=options$mbc
  
  c=default$beta$c
  if (!is.null(options$c)) c=options$c
  if (c<0 | c>1) stop("The c parameter used in the TS multiplicative bias correction technique has to be a number between 0 and 1. For more information please type ?bde")
  
  switch(normalization,
         "none" = {
           switch(mbc,
                  "none"={
                    chen99Kernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                 b=b,modified=mod)
                  },
                  "jnl" = {
                    hirukawaJLNKernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                      b=b,modified=mod)
                  },
                  "ts" = {
                    hirukawaTSKernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                     b=b,modified=mod,c=c)
                  },
{
  stop("Unrecognized mbc option. Valid options are 'none', 'jnl' and 'ts'. Note that the parameter is case sensitive. For more information please type ?bde")
})
         },
"densitywise" = {
  switch(mbc,
         "none"={
           macroBetaChen99Kernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                 b=b,modified=mod)
         },
         "jnl" = {
           macroBetaHirukawaJLNKernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                      b=b,modified=mod)
         },
         "ts" = {
           macroBetaHirukawaTSKernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                     b=b,modified=mod,c=c)
         },
{
  stop("Unrecognized mbc option. Valid options are 'none', 'jnl' and 'ts'.Note that the parameter is case sensitive. For more information please type ?bde")
})
},
"kernelwise" = {
  microBetaChen99Kernel(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                        b=b,modified=mod)
},
{
  stop("Unrecognized normalization option. Valid options are 'none', 'densitywise' and 'kernelwise'. Note that the parameter is case sensitive. For more information please type ?bde")
})
}

vitale_aux<-function(sample,x,b,options,lower.limit,upper.limit){
  biasreduced=default$vitale$biasreduction
  if (!is.null(options$biasreduced)) biasreduced=options$biasreduced
  if (!is.logical(biasreduced)) stop("The 'biasreduction' parameter used has to be a logical value. For more information please type ?bde")
  if (!biasreduced){
    if (1/b<1) {
      warning(paste("Vitales m parameter has to be an integer value, but 1/",b," is not integer. Setting it at 1. For more information please type ?bde"))
      b=1
    }
    vitale(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,m=round(1/b))
  }else{
    M=1/(2*b)
    if (!is.null(options$M)) M=options$M
    if (M<1) {
      warning(paste("Vitale's M parameter has to be an integer value. Setting it at 1. For more information please type ?bde"))
      M=1
    }
    brVitale(dataPoints=sample,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,m=round(1/b),M=M)
  }
}

boundary_aux<-function(sample,x,b,options,lower.limit,upper.limit){
  mu=default$boundary$mu
  if (!is.null(options$mu)) mu=options$mu
  if (!is.numeric(mu) & !(mu %in% 0:3)){
    stop(paste("'mu' parameter has to be an integer between 0 and 3. For more information please type ?bde"))
  }
  kernel=default$boundary$kernel
  if (!is.null(options$kernel)) kernel=options$kernel
  correct=F
  if (!is.null(options$nonegative)) correct=options$nonegative
  if (!is.logical(correct)) stop(paste("'nonegative' has to be a logic value. For more information please type ?bde"))
  switch(kernel,
         "muller94" = {
           if (!correct){
             muller94BoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                    mu=mu)
           }else{
             jonesCorrectionMuller94BoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                                   mu=mu)
           }
         },
         "muller91" = {
           if (!correct){
             muller91BoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                    mu=mu)
           }else{
             jonesCorrectionMuller91BoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                                   mu=mu)
           }
         },
         "normalized" = {
           normalizedBoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                                    mu=mu)
         },
         "none" = {
           noBoundaryKernel(dataPoints=sample,b=b,dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                            mu=mu)
         },
         {
           stop(paste("'kernel' option unrecognized, valid options are 'muller94', 'muller91', 'normalized', 'none'. For more information please type ?bde"))
         })
}

kakizawa_aux<-function(sample,x,b,options,lower.limit,upper.limit){
  if (!is.null(options$estimator)) {
    estimator=options$estimator
  }else{
    estimator=muller94BoundaryKernel(dataPoints=sample,b=b,mu=1,lower.limit = lower.limit,upper.limit = upper.limit)
  }
  if (!inherits(estimator,"BoundedDensity")) stop ("Unrecognized 'estimator' parameter. It should be an object of class 'BoundedDensity'")
  gamma=0.5
  if (!is.null(options$gamma)) gamma=options$gamma
  if (gamma<=0 | gamma>1) stop("The 'gamma' parameter has to in the interval (0,1]")
  method="b3"
  if(!is.null(options$method)) method=options$method
  switch(method,
         "b1" = {
           kakizawaB1(dataPoints=sample,m=round(1/b),dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                      gamma=gamma,estimator=estimator)
         },
         "b2" = {
           kakizawaB2(dataPoints=sample,m=round(1/b),dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                      estimator=estimator)
         },
         "b3" = {
           kakizawaB3(dataPoints=sample,m=round(1/b),dataPointsCache=x,lower.limit = lower.limit,upper.limit = upper.limit,
                      estimator=estimator)
         },
         {
           stop("Unrecognized 'method' option. It should be either 'b1', 'b2' or 'b3'")  
         })
}

launchApp<-function (...){
  shiny::runApp(system.file("App",package="bde"),...)
}