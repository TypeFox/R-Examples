MImix<-function(results, ...) UseMethod("MImix")

### Function for obtaining quantiles of mixture summary distribution
### given a set of complete data parameter estimates and associated
### variances for a number of imputed datasets.  Syntax is similar to
### the default MIcombine function in the mitools package.
###
### Arguments:
###
### results :  Parameter estimates from M completed data sets
###
### variances: Squared complete data standard errors for parameter estimates
### in results
###
### weights:  Optional weights for PMDA-2 estimates. Default setting is
### equal weights (1/M) for each imputed data set.
###
### percentiles:  Desired percentile values.  Default is for 2.5%ile, 50%ile,
### and 97.5%ile.

MImix.one<-function(results,variances,weights=1/length(results),
 		     percentiles=c(0.025,0.5,0.975),...){

  temp<-rep(0,length(percentiles))
  means<-unlist(results)
  sds<-sqrt(unlist(variances))
  for(i in 1:length(percentiles)){
    limits<-c(min(qnorm(percentiles[i],means,sds)),
    	       max(qnorm(percentiles[i],means,sds)))
     temp[i]<-uniroot(mixcdf.target,limits,means=means,sds=sds,weights=weights,
		     percentile=percentiles[i])$root
  }
  rval<-c(unlist(temp))
  names(rval)<-paste(as.character(percentiles*100),"%ile",sep="")
  rval

}


mixcdf.target<-function(x,means,sds,weights=1/length(means),percentile){
 sum(weights*pnorm(x,means,sds)) - percentile
}




MImix.default<-function(results, variances,weights=1/length(results),
 		     percentiles=c(0.025,0.5,0.975),...){

  m<-length(results)

  if (missing(variances)){
      thetas<-lapply(results, coef)
      coefnames<-names(thetas[[1]])
      numparams<-length(thetas[[1]])
      thetas<-matrix(unlist(thetas),ncol=numparams,byrow=T)
      vars<-suppressWarnings(lapply(results, vcov))
      vars<-matrix(unlist(lapply(vars,diag)),ncol=numparams,byrow=T)
    
      if(numparams==1){
        rval<-MImix.one(thetas, vars,...)
      }
      else{
       rval<-list()
       
       for(i in 1:numparams){
         rval[[i]]<-MImix.one(thetas[,i], vars[,i],...)
       }       	
       names(rval)<-coefnames
       rval
       }
  }
  else{
     rval<-MImix.one(results, variances,...)
     rval
  }
}

print.MImixresult<-function(x,...){
  cat("Multiple imputation results:\n")
  lapply(x$call, function(a) {cat("      ");print(a)})
  out<-data.frame(results=coef(x), se=sqrt(diag(vcov(x))))
  print(out)
}

summary.MImixresult<-function(object,...,alpha=0.05, logeffect=FALSE){
  cat("Multiple imputation results:\n")
  lapply(object$call, function(a) {cat("      ");print(a)})
  out<-data.frame(results=coef(object), se=sqrt(diag(vcov(object))))
  crit<-qt(alpha/2,object$df, lower=FALSE)
  out$"(lower"<-out$results-crit*out$se
  out$"upper)"<-out$results+crit*out$se
  if (logeffect){
    out$results<-exp(out$results)
    out$se<-out$se*out$results
    out$"(lower"<-exp(out$"(lower")
    out$"upper)"<-exp(out$"upper)")
  }
  out$"missInfo" <-paste(round(100*object$missinfo), "%")
  print(out,...)

}

