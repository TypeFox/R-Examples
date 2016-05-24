##########################################################################
#library(glmnet)
#library(psych)
#library(gplots)
###########################################################################
SubLasso<-function(X,y,subset,nfold){
  xnames=row.names(X)
    if (!missing(subset)) {
        jd = match(subset, xnames, 0)
        if (!all(jd > 0)) 
            stop("Some subsetd variables out of range")
    } else {jd = as.integer(0)}

    if (missing(nfold)) { nfold=5}
  x=t(X)
  n=nrow(x)
  p=ncol(x)
  fact=rep(1,p)
  fact[jd]=0
  x=as.matrix(x)

   cv.fit <- cv.glmnet(x, y, family = "binomial",nfolds=length(y))
   lambda<- cv.fit$lambda.min
   fit <- glmnet(x, y, family = "binomial",penalty.factor =fact)
   predrawy=predict(fit,x,s=lambda,type="response")
   Coef <- coef(fit, s = lambda)[-1]
   Active.Index <- which(Coef!= 0)
   Active.Coef <- Coef[Active.Index]   
   newdata=cbind(y,x[,Active.Index])

   num_tp=num_fp=num_fn=num_tn=NULL
   sn=sp=acc=avc=mcc=NULL
   set.seed(1010)
   id <- sample(1:nfold,nrow(newdata),replace=TRUE)
 for (i in 1:nfold){
   y_train=newdata[id!=i,1]
   y_test=newdata[id==i,1]
   x_train=as.matrix(newdata[id!=i,-1])
   x_test=as.matrix(newdata[id==i,-1])
   cv.fit.train<- cv.glmnet(x_train, y_train, family = "binomial")
   fit.train<- glmnet(x_train, y_train, family = "binomial")
   pfit= predict(fit.train,x_test,s = cv.fit.train$lambda.min,type="response")
   predy<-ifelse(pfit>0.5,1,0) 
   num_tp[i]=sum(y_test==1 & predy==1)
   num_fn[i]=sum(y_test==1 & predy==0)
   num_fp[i]=sum(y_test==0 & predy==1)
   num_tn[i]=sum(y_test==0 & predy==0)
   }
  sn=sum(num_tp)/sum(y==1)
  sp=sum(num_tn)/sum(y==0)
  acc=sum(num_tp+num_tn)/n
  mcc=(sum(num_tp)*sum(num_tn)-sum(num_fp)* sum(num_fn))/
      (sqrt((sum(num_tp)+sum(num_fp))*(sum(num_tp)+sum(num_fn))*(sum(num_tn)+sum(num_fp))*(sum(num_tn)+sum(num_fn))))
  valid=cbind(sn,sp,acc,mcc)

  color.map <- function(mol.biol)
    { if (mol.biol==0) "#FF0000"
      else {if (mol.biol==1) "#000FFF" 
          else "#F000FF"
           }
     }
  patientcolors <- unlist(lapply(y, color.map))
  dev.new()
  heatmap.2(t(x[,Active.Index]),col=redgreen(75), scale="row", ColSideColors=patientcolors,
           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7)

  description=describeBy(x[,Active.Index],y)
  correlation=cor(x[,Active.Index])
  dev.new()
  boxplot(predrawy~y, col="gold", xlab="Group", ylab="Linear predictor") 

 return(list(valid=valid,selname=xnames[Active.Index],w=Active.Coef,fit=fit,lambda=lambda,description=description,correlation=correlation)) 
}


