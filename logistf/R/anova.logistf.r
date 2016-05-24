anova.logistf<-function(object,  fit2, formula, method="nested", ...){
 # methods: "PLR": take difference in PLR, "nested": proper method for nested models
 # needed in logistf class: $firth, $data
 fit1<-object
 if(missing(formula)){
   if(fit1$df<fit2$df) {
    ff0<-fit1
    fit1<-fit2
    fit2<-ff0
    }
   }
 
 if(method=="PLR"){
   if(fit1$df==fit2$df) stop("Models not comparable (equal df).\n")
   df<-abs(fit1$df-fit2$df)
   PLR1<-2*diff(fit1$loglik)
   PLR2<-2*diff(fit2$loglik)
   chisq<-PLR1-PLR2
   if (chisq<0) chisq<-0
   pval<-1-pchisq(chisq,df)
   model2<-as.character(fit2$formula)
 }
 if(method=="nested"){
    f1<-fit1$formula
    a<-attr(terms(f1),"term.labels")
    if(missing(formula)){  
      f2<-fit2$formula
      b<-attr(terms(f2),"term.labels")
      # find out about which model is nested in the other
      upper<-f1
      lower<-f2
      if(!any(is.na(match(a,b)))) {
        lower<-f1
        upper<-f2
        ab<-a
        a<-b
        b<-ab
        } else if(any(is.na(match(b,a)))) stop("Models are not nested. Try method=PLR on non-nested models.\n")
      aug<-a[is.na(match(a,b))]
      f3<-paste("~",aug[1])
      if(length(aug)>1) for(j in 2:length(aug)) f3<-paste(f3, aug[j], sep="+")
      f3<-paste(f3, "-1")
      f3<-as.formula(f3)
    } else f3<-as.formula(paste(paste(as.character(formula), collapse=""),"-1",collapse=""))
      
    test<-logistftest(object=fit1, test=f3, firth=fit1$firth, weights=fit1$weights)
    chisq<-2*diff(test$loglik)
    PLR1<-2*diff(fit1$loglik)
    PLR2<-PLR1-chisq
    df<-test$df
    pval<-test$prob
    model2<-as.character(f3)
 }
 res<-list(chisq=chisq, df=df, pval=pval, call=match.call(), method=method, model1=as.character(fit1$formula), model2=model2, PLR1=PLR1, PLR2=PLR2)
 attr(res,"class")<-"anova.logistf"
 return(res)
}

  print.anova.logistf<-function(x,...){
    cat("Comparison of logistf models:\n")
    obj<-x
    pdat<-data.frame(Formula=c(paste(as.character(obj$model1)[c(2,1,3)],collapse=" "),
      if(length(as.character(obj$model2))==3) paste(as.character(obj$model2)[c(2,1,3)],collapse=" ")
      else paste("- [",paste(as.character(obj$model2)[c(2)],collapse=" "),"]")
      )
    
    
    , ChiSquared=c(obj$PLR1,obj$PLR2))
    pdatrows<-capture.output(print(pdat))
    for(i in 1:3) cat(pdatrows[i],"\n")
    cat("\nMethod: ", obj$method, "\n")
    cat("Chi-Squared: ", obj$chisq, "  df=",obj$df,"  P=", obj$pval,"\n")
  }

