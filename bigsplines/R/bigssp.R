bigssp <-
  function(formula,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,thetas=NULL,weights=NULL,
           random=NULL,remlalg=c("FS","NR","EM","none"),remliter=500,
           remltol=10^-4,remltau=NULL) {
    ###### Fits Smoothing Splines with Parametric effects
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: October 20, 2015
    
    if(class(data)[1]=="makessp"){
      sspfit <- sspwork(formula,data)
    } else{
      sspmk <- makessp(formula,data,type,nknots,rparm,
                       lambdas,skip.iter,se.fit,rseed,
                       gcvopts,knotcheck,thetas,weights,
                       random,remlalg,remliter,remltol,remltau)
      sspfit <- sspwork(formula,sspmk)
    }
    sspfit <- c(sspfit,list(call=formula))
    class(sspfit) <- "bigssp"
    return(sspfit)
    
  }