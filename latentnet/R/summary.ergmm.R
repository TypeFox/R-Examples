summary.ergmm <- function (object, point.est=c(
                                     if(!is.null(object[["mle"]])) "mle",
                                     if(!is.null(object[["sample"]])) c("pmean","mkl")),
                           quantiles=c(.025,.975),se="mle"%in%point.est,
                           bic.eff.obs=c("ties", "dyads", "actors"), ...)
{
  extraneous.argcheck(...)
  ## Just for convenience.
  sample<-object[["sample"]]
  model<-object[["model"]]
  control<-object[["control"]]
  d<-model[["d"]]
  p<-model[["p"]]
  G<-model[["G"]]
  n<-network.size(model[["Yg"]])
  sample.size<-control[["sample.size"]]
  

  summ<-list(ergmm=object,model=model)

  ## Compute the two-stage MLE point estimates.
  if("mle" %in% point.est){
    if(is.null(object[["mle"]])){
      stop("MLE was not computed for this fit.")
    }
    if(is.null(object[["mle"]][["cov"]])){
      if(model[["sender"]] || model[["receiver"]] || model[["sociality"]])
        warning("Fitting random effects as fixed effects.")
      if(se){
        ## Refit the MLE (mostly for the Hessian)
        mle<-find.mle(model,object[["mle"]],control=control,hessian=TRUE)
      }
      else mle<-object[["mle"]]

      if(d>0) mle[["Z"]]<-scale(mle[["Z"]],scale=FALSE)
      if(p>0){
        if(se){
          beta.hess<-mle[["hessian"]][1:p,1:p]
          
          beta.cov <- try(robust.inverse(-beta.hess), silent=TRUE)
          if(inherits(beta.cov,"try-error")){
            warning("Coefficient Hessian appears to be singular. Using a less accurate estimate.")
            beta.cov <- diag(1/diag(-beta.hess))
          }
        
        colnames(beta.cov) <- rownames(beta.cov) <- model[["coef.names"]]
        
        mle[["cov"]]<-beta.cov
        mle[["cor"]]<-cov2cor(beta.cov)
        
        z.values<-mle[["beta"]]/sqrt(diag(mle[["cov"]]))
        coef.table<-data.frame(mle[["beta"]],sqrt(diag(mle[["cov"]])),
                               z.values,2*pnorm(abs(z.values),0,1,lower.tail=FALSE),
                             row.names=model[["coef.names"]])
        }
        else coef.table<-data.frame(mle[["beta"]],
                                    row.names=model[["coef.names"]])
        colnames(coef.table)<-c("Estimate",if(se) "Std. Error",if(se) "z value",if(se) "Pr(>|z|)")
        mle[["coef.table"]]<-coef.table
      }

      if(model[["sender"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["sender"]])
          mle[["sender"]]<-mle[["sender"]]-mean(mle[["sender"]])
        }
        mle[["sender.var"]]<-mean(mle[["sender"]]^2)
      }
        
      if(model[["receiver"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["receiver"]])
          mle[["receiver"]]<-mle[["receiver"]]-mean(mle[["receiver"]])
        }
        mle[["receiver.var"]]<-mean(mle[["receiver"]]^2)
      }

      if(model[["sociality"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["sociality"]])
          mle[["sociality"]]<-mle[["sociality"]]-mean(mle[["sociality"]])
        }
        mle[["sociality.var"]]<-mean(mle[["sociality"]]^2)
      }

      if(model[["dispersion"]]){
        mle[["dispersion"]]<-mle[["dispersion"]]
      }
      
      object[["mle"]]<-mle
    }
    summ[["mle"]]<-object[["mle"]]
  }

  ## Compute the posterior mean point estimates.
  if("pmean" %in% point.est){
    if(is.null(object[["sample"]])){
      stop("MCMC was not was not run for this fit.")
    }
    if(is.null(object[["pmean"]])){
      pmean<-list()
      for(name in names(sample)){
        if(is.null(sample[[name]])) next
        name.dim<-length(dim(sample[[name]]))
        pmean[[name]]<-{
          if(name.dim<2) mean(sample[[name]])
          else if(name.dim==2) apply(sample[[name]],2,mean)
          else if(name.dim==3) apply(sample[[name]],2:3,mean)
        }
      }
      if(G>0){
        pmean[["Z.pZK"]]<-t(apply(sample[["Z.K"]],2,tabulate,G))/sample.size
        pmean[["Z.K"]]<-apply(pmean[["Z.pZK"]],1,which.max)
      }
    
      beta.cov<-cov(sample[["beta"]])
      colnames(beta.cov) <- rownames(beta.cov) <- model[["coef.names"]]
      
      pmean[["cov"]]<-beta.cov
      pmean[["cor"]]<-cov2cor(beta.cov)
      
      beta.q0<-apply(sample[["beta"]],2,function(x) min(mean(x<=0),mean(x>=0))*2)
      
      coef.table<-data.frame(pmean[["beta"]],
                           t(apply(sample[["beta"]],2,function(x)quantile(x,quantiles))),
                           beta.q0,
                           row.names=model[["coef.names"]])
      colnames(coef.table)<-c("Estimate",paste(quantiles*100,"%",sep=""),"2*min(Pr(>0),Pr(<0))")
      pmean[["coef.table"]]<-coef.table
      object[["pmean"]]<-pmean
    }
    summ[["pmean"]]<-object[["pmean"]]
  }
  ## Compute the MKL point estimates.
  if("mkl" %in% point.est){
    if(is.null(object[["mkl"]])){
      stop("MKL was not produced for this fit.")
    }
    mkl<-summ[["mkl"]]<-object[["mkl"]]
    coef.table<-data.frame(mkl[["beta"]],
                           row.names=model[["coef.names"]])
    colnames(coef.table)<-c("Estimate")
    summ[["mkl"]][["coef.table"]]<-coef.table
  }

  if("pmode" %in% point.est){
    if(is.null(object[["pmode"]])){
      stop("Conditional posterior mode was not computed for this fit.")
    }
    summ[["pmode"]]<-object[["pmode"]]
  }

  if(!is.null(object[["mkl"]]) && !is.null(bic.eff.obs)){
    summ[["bic"]]<-bic.ergmm(object, eff.obs = bic.eff.obs, ...)
  }

  class(summ)<-'summary.ergmm'
  summ
}

print.summary.ergmm<-function(x,...){
  ## For convenience
  model<-x[["model"]]
  control<-x[["ergmm"]][["control"]]
  
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")

  cat("Formula:   ")
  print(model[["formula"]])
  cat("Attribute: ")
  if(is.null(model[["response"]])) cat("edges") else cat(model[["response"]])
  cat("\n")
  cat("Model:    ",model[["family"]],"\n")
  
  digits = max(3, getOption("digits") - 3)
  
  if(!is.null(x[["pmean"]])) cat ("MCMC sample of size ", control[["sample.size"]], ", draws are ",
       control[["interval"]]," iterations apart, after burnin of ",control[["burnin"]], " iterations.\n",sep="")
       
  if(!is.null(x[["pmean"]])){
    cat("Covariate coefficients posterior means:\n")
    printCoefmat(as.matrix(x[["pmean"]][["coef.table"]]),P.values=TRUE,has.Pvalue=TRUE)
    cat("\n")
    if(!is.null(x[["pmean"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["pmean"]][["dispersion"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["sender.var"]]))
      cat("Sender effect variance: ",x[["pmean"]][["sender.var"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["receiver.var"]]))
      cat("Receiver effect variance: ",x[["pmean"]][["receiver.var"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["sociality.var"]]))
      cat("Sociality effect variance: ",x[["pmean"]][["sociality.var"]],".\n", sep="")
  }

  if(!is.null(x[["mle"]])){
    cat("Covariate coefficients MLE:\n")
    printCoefmat(as.matrix(x[["mle"]][["coef.table"]]),P.values=length(names(x[["mle"]][["coef.table"]]))>1)
    cat("\n")
    if(!is.null(x[["mle"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["mle"]][["dispersion"]],".\n", sep="")
  }
  if(!is.null(x[["bic"]])){
    cat("Overall BIC:       ", x[["bic"]][["overall"]],"\n")
    cat("Likelihood BIC:    ", x[["bic"]][["Y"]],"\n")
    if(model[["d"]]>0){
      cat("Latent space/clustering BIC:    ", x[["bic"]][["Z"]],"\n")
    }
    if(model[["sender"]]){
      cat("Sender effect BIC:    ", x[["bic"]][["sender"]],"\n")
    }
    if(model[["receiver"]]){
      cat("Receiver effect BIC:    ", x[["bic"]][["receiver"]],"\n")
    }
    if(model[["sociality"]]){
      cat("Sociality effect BIC:    ", x[["bic"]][["sociality"]],"\n")
    }
    cat("\n")
  }

  if(!is.null(x[["mkl"]])){
    cat("Covariate coefficients MKL:\n")
    print(x[["mkl"]][["coef.table"]])
    cat("\n")
    if(!is.null(x[["mkl"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["mkl"]][["dispersion"]]," (probably invalid).\n", sep="")
    cat("\n")
  }
  if(!is.null(x[["pmode"]])){
    cat("Covariate coefficients posterior mode:\n")
    print(x[["pmode"]][["coef.table"]])
    cat("\n")
    if(!is.null(x[["pmode"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["pmode"]][["dispersion"]],".\n", sep="")
    cat("\n")
  }
}

bic.ergmm<-function(object, eff.obs=c("ties", "dyads", "actors"), ...){
  extraneous.argcheck(...)
  if(is.null(object[["mkl"]])){
    stop("MKL estimates were not computed for this fit.")
  }

  model<-object[["model"]]

  n<-network.size(model[["Yg"]])
  
  if(is.character(eff.obs)){
    eff.obs <- switch(match.arg(eff.obs),
                      ties = {if(!all(match(model[["Ym"]], c(0,1)) > 0, na.rm=TRUE)) warning('Number of "ties" in a valued network may not be well-defined.'); network.edgecount(model[["Yg"]])},
                      dyads = network.dyadcount(model[["Yg"]]),
                      actors = n)
  }
  
  
  condZRE<-with(object,find.mle(model,mkl,given=list(Z=mkl[["Z"]],sender=mkl[["sender"]],receiver=mkl[["receiver"]],sociality=mkl[["sociality"]]),control=object[["control"]]))

  bic<-with(model,list(Y = -2*condZRE[["lpY"]] + (p)*log(eff.obs),
                              Z =
                              if(d>0){
                                if(G>0){
                                  mbc.llk<-Inf
                                  Gsub<--1
                                  while(!is.finite(mbc.llk)){
                                    Gsub<-Gsub+1
                                    mbc.llk<-mbc.VII.EM(G-Gsub,object[["mkl"]][["Z"]])[["llk"]]
                                  }
                                  if(Gsub) warning(paste("Bad clustering: treating",Gsub,"clusters as empty."))
                                  -2*mbc.llk+((G-Gsub)-1 + d*(G-Gsub) + (G-Gsub))*log(n)
                                } else {
                                  -2*sum(dnorm(object[["mkl"]][["Z"]],0,sqrt(mean(condZRE[["Z"]]^2)*d),log=TRUE))+1*log(n*d)
                                }
                              } else 0,
                              sender=if(sender) -2*sum(dnorm(condZRE[["sender"]],0,sqrt(mean(condZRE[["sender"]]^2)),log=TRUE))+1*log(n) else 0,
                              receiver=if(receiver) -2*sum(dnorm(condZRE[["receiver"]],0,sqrt(mean(condZRE[["receiver"]]^2)),log=TRUE))+1*log(n) else 0,
                              sociality=if(sociality) -2*sum(dnorm(condZRE[["sociality"]],0,sqrt(mean(condZRE[["sociality"]]^2)),log=TRUE))+1*log(n) else 0
                              )
            )

  if(!.latentnetEnv$BIC.warned){
    message("NOTE: It is not certain whether it is appropriate to use latentnet's BIC to select latent space dimension, whether or not to include actor-specific random effects, and to compare clustered models with the unclustered model.")
    .latentnetEnv$BIC.warned <- TRUE
  }
  
  bic[["overall"]]<-sum(unlist(bic))
  
  bic
}
