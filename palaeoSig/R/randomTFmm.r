randomTFmm <-
function(fos,modelList,col, condition, ord=rda,...){#reconstruct random variables
  partial<-!missing(condition)
  
  if(!partial){
    PC<-ord(fos)
  }else{
    form1<-formula(paste("fos~1+Condition(",paste(names(condition), collapse="+"),")"))
    PC<-ord(form1, data=condition)
  }
  MAX<- PC$CA$eig[1]/PC$tot.chi
  
  nenv<-attr(modelList,"Nenv")
  obs<-lapply(modelList[1:nenv],function(mod){
    Pred<-predict(mod,fos)
    if(is.list(Pred))Pred<-Pred$fit[,col]
    if(!partial){
      RDA<-ord(fos~Pred)
    }else{
      form<-formula(paste("fos~Pred+Condition(",paste(names(condition), collapse="+"),")"))
      RDA<-ord(form, data=condition)
    }
    EX<- RDA$CCA$tot.chi/RDA$tot.chi
    EIG1<-RDA$CA$eig[1]/RDA$tot.chi
    list(EX=EX, pred=Pred, EIG1=EIG1, clas=class(Mod), mod=Pred)
  })

  sim.ex<-sapply(modelList[-(1:nenv)],function(m){
       p<-predict(m,fos)
       if(is.list(p))p<-p$fit[,col]
       if(!partial){
          r<-ord(fos~p)
       }else{
          form<-formula(paste("fos~p+Condition(",paste(names(condition), collapse="+"),")"))
          r<-ord(form, data=condition)
       }
       r$CCA$tot.chi/r$tot.chi
  })
  preds<-lapply(obs,function(x)x$pred)
  EXs<-sapply(obs,function(x)x$EX)
  eig1<-sapply(obs,function(x)x$EIG1)
  sig<-sapply(EXs,function(E)mean(E<=c(E,sim.ex)))
  res<-list(PCA=PC,preds=preds, MAX=MAX, EX=EXs, eig1=eig1,sim.ex=sim.ex, sig=sig)
  class(res)<-"palaeoSig"
  return(res)
}

