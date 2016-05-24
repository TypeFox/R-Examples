extractAIC.logistf<-function(fit, scale, k=2, ...){
  dev<--2*diff(fit$loglik)
  AIC<-dev+k*fit$df
  edf<-fit$df
  return(c(edf,AIC))
}

nobs.logistf<-function(object, ...){
  return(object$n)
}

drop1.logistf<-function(object, scope, test="PLR", ...){
   variables<-attr(terms(object$formula),"term.labels")
   nvar<-length(variables)
   mat<-matrix(0,nvar,3)
   for(i in 1:nvar){
      newform<-as.formula(paste("~",variables[i]))
      res<-anova(object, formula=newform)
      mat[i,1]<-res$chisq
      mat[i,2]<-res$df
      mat[i,3]<-res$pval
    }
    rownames(mat)<-variables
    colnames(mat)<-c("ChiSq","df","P-value")
    return(mat)
}



terms.logistf<-function(x, ...){
   object<-x
   options(warn=-1)
   fakeglm<-try(glm(formula=object$formula,data=object$data, family=binomial,maxit=1))
   options(warn=0)
   return(terms(fakeglm))
}

add1.logistf<-function(object, scope, test="PLR", ...){
   if(missing(scope)) scope<-colnames(object$data)
   if(is.numeric(scope)) scope<-colnames(object$data)[scope]
   whereisresponse<-(match(colnames(model.frame(object$formula,data=object$data))[1],scope))
   if(!is.na(whereisresponse)) scope<-scope[-whereisresponse]
   scope<-scope[is.na(match(scope, attr(terms(object),"term.labels")))]
   variables<-scope
   
   nvar<-length(variables)
   mat<-matrix(0,nvar,3)
   for(i in 1:nvar){
      newform<-as.formula(paste("~.+",variables[i]))
      res<-anova(object, update(object,formula=newform))
      mat[i,1]<-res$chisq
      mat[i,2]<-res$df
      mat[i,3]<-res$pval
    }
    rownames(mat)<-variables
    colnames(mat)<-c("ChiSq","df","P-value")
    return(mat)
}




forward<-function(object, scope, steps=1000, slentry=0.05, trace=TRUE, printwork=FALSE, pl=TRUE, ...){
  istep<-0
  working<-object
  if(missing(scope)) scope<-colnames(object$data)
  if(is.numeric(scope)) scope<-colnames(object$data)[scope]
  scope<-scope[-(match(colnames(model.frame(object$formula,data=object$data))[1],scope))]
  if(trace){
    cat("Step ", istep, ": starting model\n")
    if(printwork) {
        print(working)
        cat("\n\n")
        }
    }
  if(missing(scope)) stop("Please provide scope (vector of variable names).\n")
  inscope<-scope
  while(istep<steps & length(inscope)>=1){
    istep<-istep+1
    mat<-add1(working, scope=inscope)
    #if(printwork)   print(mat)
    if(all(mat[,3]>slentry)) break
    index<-(1:nrow(mat))[mat[,3]==min(mat[,3])]
    if(length(index)>1) index<-index[mat[index,1]==max(mat[index,1])]
    addvar<-rownames(mat)[index]
    newform=as.formula(paste("~.+",addvar))
    working<-update(working, formula=newform, pl=FALSE)
    newindex<-is.na(match(inscope, attr(terms(object),"term.labels")))
    if(all(is.na(newindex)==TRUE)) break
    inscope<-inscope[-index]
    if(trace){
      cat("Step ", istep, ": added ", addvar, " (P=", mat[addvar,3],")\n")
      if(printwork) {
        print(working)
        cat("\n\n")
        }
      }
    }
   if(pl) working<-update(working, pl=TRUE)
   if(trace) cat("\n")
   return(working)
   }





backward<-function(object, scope, steps=1000, slstay=0.05, trace=TRUE, printwork=FALSE, ...){
  istep<-0
  working<-object
  if(trace){
    cat("Step ", istep, ": starting model\n")
    if(printwork) {
        print(working)
        cat("\n\n")
        }
    }
  if(missing(scope)) scope<-attr(terms(working),"term.labels")
#    if(missing(scope)) scope<-names(coef(working))
  
  while(istep<steps & working$df>=1){
    istep<-istep+1
    mat<-drop1(working)
    if(all(mat[,3]<slstay)) break
    inscope<-match(scope,rownames(mat))
    inscope<-inscope[!is.na(inscope)]
    removal<-rownames(mat)[mat[,3]==max(mat[inscope,3])]
    newform=as.formula(paste("~.-",removal))
    if(working$df==1 | working$df==mat[mat[,3]==max(mat[,3]),2])  working<-update(working, formula=newform, pl=FALSE)
    else working<-update(working, formula=newform)
    if(trace){
      cat("Step ", istep, ": removed ", removal, " (P=", max(mat[,3]),")\n")
      if(printwork) {
        print(working)
        cat("\n\n")
        }
      }
    }
   if(trace) cat("\n")
   return(working)
   }
  
  
  
  
