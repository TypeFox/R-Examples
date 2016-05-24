
"ahazisis"<-function(surv,X,weights,standardize=TRUE,nsis=floor(nobs/1.5/log(nobs)),
                   do.isis=TRUE, maxloop=5, penalty=sscad.control(),
                   tune=cv.control(), rank=c("FAST","coef","z","crit"))
{
  this.call<-match.call()
  nobs<-nrow(X)

  detail.pickind<-detail.ISISind<-detail.ISIScoef<-list()

  if(nsis<0)
    stop("invald 'nsis'")
  if(nsis>ncol(X))
    nsis<-ncol(X)
  if(maxloop<1)
    stop("invalid 'maxloop'")
  maxloop<-as.integer(maxloop)
  
  rank<-match.arg(rank)
  # INITIAL RECRUITMENT
  m<-ahaz(surv=surv,X=X,weights=weights,univariate=TRUE);
  initRANKorder<-switch(rank,
                        FAST=order(abs(m$d),decreasing=TRUE),
                        coef=order(abs(m$d/m$D),decreasing=TRUE),
                        z=order(m$d^2/m$B,decreasing=TRUE),
                        crit=order(m$d^2/m$D,decreasing=TRUE))

                                        # VARIABLE SELECTION STEP
  pick.ind<-sort(initRANKorder[1:floor(2*nsis/3)])
  detail.pickind[[1]]<-pick.ind
  SISind<-initRANKorder[1:nsis]
  
  if(!do.isis){
    out<-list(call=this.call,initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
              detail.ISISind = detail.ISISind, SISind = SISind, ISISind = NULL, 
              nsis = nsis,do.isis=do.isis,maxloop=maxloop)
    class(out)<-"ahazisis"
    return(out)
  }

  # Switch to 't'-ranking if initial
  if(rank=="FAST")
    rank<-"z"

  pentune<-tune.ahazpen(surv=surv,X=X[,pick.ind],weights=weights,standardize=standardize,penalty=penalty,tune=tune)
  coef<-as.numeric(coef(pentune))

  ISISind<-pick.ind[coef!=0]
  ISIScoef<-coef[coef!=0]
  
  detail.ISISind[[1]]<-ISISind
  detail.ISIScoef[[1]]<-ISIScoef

    
  # Exit gracefully if no variables selected
  if(length(ISISind)==0){
            warnings("no covariates selected by ISIS.")
                  out<-list(call=this.call,initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
              detail.ISISind = detail.ISISind, detail.ISIScoef = detail.ISIScoef, SISind = SISind, ISISind = ISISind, 
              ISIScoef = ISIScoef, nsis = nsis,do.isis=do.isis,maxloop=maxloop)
                class(out)<-"ahazisis"
                return(out)
  }

  # REPEAT VARIABLE SELECTION/RE-RECRUITMENT
  if(maxloop>1)
    {
      oldISISind<-NULL
      for(i in 2:maxloop)
        {        
          if(length(ISISind)!=nsis && !setequal(oldISISind, ISISind))
            {
              oldISISind<-ISISind
        
        
                                        # CALCULATE 'ADJUSTED' RELEVANCE
              s<-ahaz.adjust(surv=surv,X=X,weights=weights,idx=oldISISind,method=rank)
              s$adj[oldISISind]<-0
         
                                        # CONSTRUCT NEW SET OF POTENTIALLY RELEVANT FEATURES
              vv<-setdiff(1:ncol(X),ISISind)
              new.pickind<-sort(vv[order(abs(s$adj[vv]),decreasing=TRUE)[1:(nsis-length(ISISind))]])
              pick.ind = sort(c(ISISind, new.pickind))

                                        # SELECT FEATURES
              pentune<-tune.ahazpen(surv=surv,X=X[,pick.ind],weights=weights,standardize=standardize,penalty=penalty,tune=tune)
         coef<-as.numeric(coef(pentune))
             
              
              ISISind<-pick.ind[coef!=0]
              ISIScoef<-coef[coef!=0]
            
              detail.pickind[[i]]<-pick.ind
              detail.ISISind[[i]] <- ISISind
              detail.ISIScoef[[i]]<-ISIScoef
              
              if(all(coef==0)){
                warnings("No covariates selected by ISIS.")
                  out<-list(call=this.call,initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
                            detail.ISISind = detail.ISISind, detail.ISIScoef=detail.ISIScoef, SISind = SISind, ISISind = ISISind, 
                            ISIScoef = ISIScoef, nsis = nsis,do.isis=do.isis,maxloop=maxloop)
                class(out)<-"ahazisis"
                return(out)
              }
            }
     }
    }
  out<-list(call=this.call,initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
            detail.ISISind = detail.ISISind,detail.ISIScoef=detail.ISIScoef, SISind = SISind, ISISind = ISISind, 
            ISIScoef = ISIScoef, nsis = nsis,do.isis=do.isis,maxloop=maxloop)
  class(out)<-"ahazisis"
  return(out)
}

"print.ahazisis"<-function(x,  digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCall: ", deparse(x$call), "\n\n")
    if(!x$do.isis){
      cat("Number of covariates selected by SIS:    ",x$nsis,"\n")
    } else {
          cat("Number of covariates selected by SIS:   ",x$nsis,"\n\n")
          cat("ISIS iterations actually used:          ", length(x$detail.ISISind),"\n")
          cat("Number of covariates selected by ISIS:  ",length(x$ISISind),"\n\n")
        }
    invisible(x)
  }
