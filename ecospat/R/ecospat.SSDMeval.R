ecospat.SSDMeval <- function (eval,pred,proba,ntir)
{
     if(proba==F){
          ntir<-1
          pred2<-pred
     }

     Rdev<-over<-under<-succ<-sens<-spec<-kappa<-tss<-sim<-jac<-
          matrix(nrow=nrow(eval),ncol=ntir,dimnames=list(rownames(eval),c(1:ntir)))

     for(i in 1:ntir){
          if(proba==T){
               pred2<-matrix(nrow=nrow(pred),ncol=ncol(pred))
               for(k in 1:nrow(pred)){
                    pred2[k,]<-rbinom(n=ncol(pred),size=1,prob=as.numeric(pred[k,]))
               }
          }
          errors<-2*pred2+eval
          # if 0 true negative
          # if 1 false negative
          # if 2 false positive
          # if 3 true positive
          for (j in 1:nrow(errors)){
               a<-length(which(errors[j,]==3))
               b<-length(which(errors[j,]==2))
               c<-length(which(errors[j,]==1))
               d<-length(which(errors[j,]==0))
               n<-ncol(errors)

               Rdev[j,i]<-sum(eval[j,])-sum(pred2[j,])
               over[j,i]<- b/(b+d)
               under[j,i]<- c/(a+c)
               succ[j,i]<- (a+d)/n
               sens[j,i]<-  a/(a+b)
               spec[j,i]<-  d/(c+d)
               kappa[j,i]<-  (((a+d)/n)-(((a+c)*(a+b)+(b+d)*(d+c))/(n^2)))/(1-(((a+c)*(a+b)+(b+d)*(d+c))/(n^2)))    ## Kappa
               tss[j,i]<-  sens[j,i]+spec[j,i]-1
               sim[j,i]<-(2*a)/(2*a+b+c)
               jac[j, i]<- a/(a+b+c)
          }
     cat("trial",i,"on",ntir,"\n")
     }
     res<-list(Rdev,over,under,succ,sens,spec,kappa,tss,sim,jac)
     names(res)<-c("deviation.rich.pred",
                   "overprediction",
                   "underprediction",
                   "prediction.success",
                   "sensitivity",
                   "specificity",
                   "kappa",
                   "TSS",
                   "similarity",
                   "Jaccard")
     return(res)
}