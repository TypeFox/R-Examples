ma.aps.reg <-
function(dv,iv,min=1,max,mad=FALSE,aic=FALSE,bic=FALSE,model.sig=TRUE,coeff.sig=TRUE,coeff.vif=TRUE,coeff.cor=FALSE){
  
  # Initial checks
  if(ncol(as.matrix(dv))!=1){stop('The number of dependent variable is not equal to one.')}
  if(nrow(as.matrix(dv))!=nrow(as.matrix(iv))){stop('The number of independent variables and dependent variable is not the same.')}
  if(min<1){stop('The minimum number of independent variables must be greater than one.')}
  dv<-as.matrix(dv)
  
  # Load libraries
  # library(combinat);library(car)
  
  # Parameters
  fit<-NULL # Nullifying variable to avoid non-binding issue
  r<-sum(mad,aic,bic)
  
  # Run
  for(k in min:max){
    c<-combn(ncol(as.matrix(iv)),k);t<-dim(c)[2];index<-t(c)
    results.k<-matrix(NA,t,r+2);rownames(results.k)<-c(1:t)  
    colnames(results.k)<-c("Model","Adj.R2",if(mad)"MAD",if(aic)"AIC",if(bic)"BIC")
    for(i in 1:t){
      if(i==1){j<-1}
      model<-paste0("fit<-lm(dv[,1] ~ iv[,",paste(index[i,],collapse="] + iv[,"),"])")
      model.check<-paste0("alias(",model,")$Complete")
      if(is.null(eval(parse(text=model.check)))){
        eval(parse(text=model))
        
        # is the model significant?
        ifelse(model.sig,model.sig.t<-anova(fit)$Pr[1]<0.05,model.sig.t<-TRUE)
        # are coefficients significant?
        ifelse(coeff.sig,coeff.sig.t<-sum((summary(fit)$coefficients[,4]<0.05)==FALSE)==0,coeff.sig.t<-TRUE) 
        # is multicollinearity allowable?
        if(coeff.vif && ncol(as.matrix(iv))>1){coeff.vif.t<-ifelse(k>1,sum(vif(fit)>10)==0,TRUE)}else{coeff.vif.t<-TRUE}
        # are coefficidents consistent with bivariate cor?
        ifelse(coeff.cor,coeff.cor.t<-sum(cor(iv[,index[i,]],dv)*fit$coefficients[-1]<0)==0,coeff.cor.t<-TRUE)
        
        if(model.sig.t && coeff.sig.t && coeff.cor.t && coeff.vif.t){
          rownames(results.k)[j]<-c(paste0("Model-",k,"-",i,"/",t))  
          if(!is.null(colnames(dv)) && !is.null(colnames(iv))){results.k[j,1]<-paste(colnames(dv),"~",paste0(colnames(iv)[index[i,]],collapse=" + "))}
          else{results.k[j,1]<-substr(model,9,nchar(model)-1)}
          results.k[j,2]<-round(summary(fit)$adj.r.squared,4)
          if(r>0){
            for(p in 1:r){
              if(colnames(results.k)[2+p]=="MAD"){results.k[j,2+p]<-round(mean(abs(fit$fitted.values-dv[,1])),4)}
              if(colnames(results.k)[2+p]=="AIC"){results.k[j,2+p]<-round(AIC(fit),4)}
              if(colnames(results.k)[2+p]=="BIC"){results.k[j,2+p]<-round(BIC(fit),4)}
            }  
          }
          j<-j+1
        }
      }
      if(i==t){results.k<-results.k[1:(j-1),,drop=FALSE]}
    }
    if(k==min){results<-results.k}else if(!is.na(results.k[1,1])){results<-rbind(results,results.k)}
    if(k==max){results<-noquote(results)}
  }
  if(is.na(results[1,1])){stop('There is no model meeting given criteria.')}
  return(results)
}
