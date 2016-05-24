## Distributed Evaluation of Local Sensitivity Analysis
## These functions calculate first order parameter sensitivity using DELSA, (Rakovec et al., 2014, WRR): 
## Rakovec, O., M. C. Hill, M. P. Clark, A. H. Weerts, A. J. Teuling, R. Uijlenhoet (2014), Distributed 
##                         Evaluation of Local Sensitivity Analysis (DELSA), with application to hydrologic models, 
##                         Water Resour. Res., 50, 1-18, doi:10.1002/2013WR014063.
##
## Author: Oldrich Rakovec  
## email: oldrich.rakovec@ufz.de; oldrich.rakovec@gmail.com, 
## date: 12 May 2014
##
## Converted to functions compatible with sensitivity package
## Joseph Guillaume, josephguillaume@gmail.com, 23/5/2014
##
## Calling a model from R
##
## obj <- delsa(model=sobol.fun,par.ranges=replicate(8,c(0,1),simplify=FALSE),samples=100,method="sobol")
## print(obj)
## plot(obj)
##
## Calling the model separately
##
## obj <- delsa(model=NULL,par.ranges=replicate(8,c(0,1),simplify=FALSE),samples=100,method="sobol")
## ## Separately: Evaluate x$X and load the result in y
## ## Now calculate indices using those results:
## tell(obj,y)
## print(obj)
## plot(obj)

## suppress R CMD check NOTE 'no visible binding for global variable' for ggplot2
if (getRversion() >= "2.15.1") utils::globalVariables(c("value","variable","y"))

delsa <- function(
  model=NULL,
  perturb=1.01,
  #passed to parameterSets
  par.ranges,
  samples,
  method,
  ## Alternatively, directly define prior variance 
  ##  and locations to evaluate sensitivity
  X0,
  varprior,
  ## Additional arguments to model
  ...
){
  if(missing(X0)){
    X0=parameterSets(par.ranges,samples,method)
    ## Prior variance for compatibility with Sobol. Rakovec et al, 2014, WRR of Sect. 3.2 & Appendix B
    ## Only the ratio of variances matters
    varprior=sapply(par.ranges,diff)^2/12 
  } else if (!missing(X0) & missing(varprior)){
    stop("If X0 is specified, varprior must also be specified")
  }
  
  
  ## 1% perturbation of each par
  X=do.call(rbind,
                  lapply(1:ncol(X0),
                      function(i){
                        X2i=X0
                        X2i[,i]=X2i[,i]*perturb
                        X2i
                      }))
  ## Add the base parameters first, such that X= nrow(X0)*(1+ncol(X0)),ncol(X0)
  ## with base, and then perturbation on each parameter in turn
  X=rbind(X0,X)
    
  x <- list(model = model, X0 = X0, X=X,varprior=varprior,call = match.call())
  class(x)<-"delsa"
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
  return(x)
}

tell.delsa<-function(x,y=NULL,...){
  
  id <- deparse(substitute(x))
  if (!is.null(y)) {
    x$y <- y
  }
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  ## initialization
  Kpar=ncol(x$X0)
  Nsamp=nrow(x$X0)
  vartot = rep(0,Nsamp)
  delsafirst=deriv=varfir=matrix(NA, ncol=Kpar,nrow=Nsamp)
  
  out <- as.numeric(x$y)
  
  for (rsamp in 1:Nsamp){    ##looping over parameter sets 
    
    for (jpar in 1: Kpar){  ##looping over model parameters
      idx.pert=Nsamp*jpar+rsamp
      
      ## calculate derivative (see eq. 10)
      deriv[rsamp,jpar] = (out[idx.pert] -  out[rsamp]) / (x$X[idx.pert,jpar] - x$X[rsamp,jpar])
      
      ## calculate local first order variance (see nominator of eq. 13)
      varfir[rsamp,jpar] =  (deriv[rsamp,jpar]**2)*(x$varprior[jpar])
      
      ## calculate local equation total variance (see eq. 12)
      vartot[rsamp] = vartot[rsamp] + varfir[rsamp,jpar]
      
      if (jpar==Kpar){
        ## first order Local first order sensitivity index (see eq. 13 in Rakovec et al 2014)  
        for (jjpar in 1:Kpar)
          delsafirst[rsamp,jjpar] =  varfir[rsamp,jjpar] /  vartot[rsamp] 
      } ## end of if statements
      
    } ## End loop over model parameters
    
  }## End loop over parameter sets
  
  colnames(delsafirst)=colnames(x$X)
  
  x$delsafirst=delsafirst
  assign(id, x, parent.frame()) 
}

print.delsa<-function (x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nLocations calculated:", nrow(x$X0), "\n")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nSummary of first order indices across parameter space:\n")
    print(summary(x$delsafirst))
  }
}

plot.delsa <- function(x,which=1:3,
                       ask = dev.interactive(),
                       ...){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)){ 
    stop('The package ggplot2 is missing, but is required for the plots')
  }
  if (!requireNamespace("reshape2", quietly = TRUE)){ 
    stop('The package reshape2 is missing, but is required for the plots')
  }
  if (requireNamespace("ggplot2", quietly = TRUE)){ 
    if (requireNamespace("reshape2", quietly = TRUE)){ 
      
  obj=x ##must be called x in the arguments for compatibility
  
  show <- rep(FALSE, 3)
  show[which] <- TRUE
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  if(show[1L]){
    # Figure 11 - cdfs
    temp=as.data.frame(obj$delsafirst)
    temp$id<-1:nrow(temp)
    temp=reshape2::melt(temp,id.var="id")
    print(
      ggplot2::ggplot(data=temp,ggplot2::aes(x=value,colour=variable))+
        ggplot2::stat_ecdf()+
        ggplot2::scale_x_continuous("DELSA results for first order sensitivity")+
        ggplot2::scale_y_continuous("Cum. frequency")+
        ggplot2::labs(title="CDF of first order sensitivity across parameter space")
    )
  }
  
  if(any(show[2:3L])){
    temp=as.data.frame(obj$delsafirst)
    temp$id<-1:nrow(temp)
    temp=reshape2::melt(temp,id.var="id")
    stopifnot(identical(obj$X0,obj$X[1:nrow(obj$X0),]))
    temp$y<-obj$y[temp$id] #first set of evaluations is at X0 position
    ## Add parameter value
    temp2=as.data.frame(obj$X0)
    temp2$id<-1:nrow(temp2)
    temp2=reshape2::melt(temp2,id.var="id")
    temp2$x <- temp2$value
    temp2$value<-NULL
    temp=merge(temp,temp2)
    #columns: id,variable,value,y,x
  }
  if(show[2L]){
    # Similar to Figure 12 - related to model output
    print(
      ggplot2::ggplot(data=temp)+
        ggplot2::geom_point(ggplot2::aes(x=value,y=y))+
        ggplot2::scale_x_continuous(name="DELSA first order sensitivity")+
        ggplot2::scale_y_continuous(name="Model output")+
        ggplot2::facet_wrap(~variable,scales="free")+
        ggplot2::labs(title="First order sensitivity as related to model response")
    )
  }
  if(show[3L]){
    # Relating sensitivity to parameter value
    print(
      ggplot2::ggplot(data=temp)+
        ggplot2::geom_point(ggplot2::aes(y=value,x=x,colour=y))+
        ggplot2::scale_y_continuous(name="DELSA first order sensitivity")+
        ggplot2::scale_x_continuous(name="Parameter value")+
        ggplot2::scale_color_continuous(name="Model response")+
        ggplot2::facet_wrap(~variable,scales="free")+
        ggplot2::labs(title="First order sensitivity as as related to parameter value")
    )
  }
  invisible()
  }
  }
  
}
