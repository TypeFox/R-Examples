bigssg <-
  function(formula,family,data=NULL,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.lp=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL,
           gcvtype=c("acv","gacv","gacv.old")) {
    ###### Fits Generalized Smoothing Spline ANOVA models
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: July 20, 2015
    
    if(class(data)[1]=="makessg"){
      ssgfit <- ssgwork(formula,data)
    } else{
      ssgmk <- makessg(formula,family,data,type,nknots,rparm,
                       lambdas,skip.iter,se.lp,rseed,gcvopts,
                       knotcheck,gammas,weights,gcvtype)
      ssgfit <- ssgwork(formula,ssgmk)
    }
    ssgfit <- c(ssgfit,list(call=formula))
    class(ssgfit) <- "bigssg"
    return(ssgfit)
    
  }