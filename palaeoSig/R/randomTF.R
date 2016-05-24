randomTF <-
function(spp,env,fos,n=99,fun,col, condition, autosim, ord=rda,...){#reconstruct random variables
  if(!is.list(env))env<-list(env=env)
  rownames(spp)<-1:nrow(spp)
  if(identical(names(formals(MAT)),names(formals(fun)))){
      mod1<-predict(MAT(spp,env[[1]],...),fos)
      analogues<-unique(as.vector(as.numeric(mod1$match.name)))
      spp<-spp[analogues,]
      env<-env[analogues,,drop=FALSE]
      rownames(spp)<-1:nrow(spp)
  }
  partial<-!missing(condition)
  
  if(!partial){
    PC<-ord(fos)
  }else{
    form1<-formula(paste("fos~1+Condition(",paste(names(condition), collapse="+"),")"))
    PC<-ord(form1, data=condition)
  }
  MAX<- PC$CA$eig[1]/PC$tot.chi
  
  obs<-lapply(env,function(ev){
    Mod<-fun(spp, ev,...)
    Pred<-predict(Mod,fos)
    if(is.list(Pred))p<-Pred$fit[,col]
    else p<-Pred
    if(!partial){
      RDA<-ord(fos~p)
    }else{
      form<-formula(paste("fos~p+Condition(",paste(names(condition), collapse="+"),")"))
      RDA<-ord(form, data=condition)
    }
    EX<- RDA$CCA$tot.chi/RDA$tot.chi
    EIG1<-RDA$CA$eig[1]/RDA$tot.chi
    list(EX=EX, pred=p, EIG1=EIG1, clas=class(Mod), mod=Pred)
  })
  if(missing(autosim))rnd<-matrix(runif(nrow(spp)*n),ncol=n)
  else rnd<-autosim
  if(obs[[1]]$clas=="MAT"){
         p<-t(apply( apply(obs[[1]]$mod$m,2,as.numeric),1,function(n)colMeans(rnd[n,])))
         sim.ex<-apply(p,2,function(pp){
            if(!partial){
                r<-ord(fos~pp)
            }else{
                form<-formula(paste("fos~pp+Condition(",paste(names(condition), collapse="+"),")"))
                r<-ord(form, data=condition)
            }
            r$CCA$tot.chi/r$tot.chi                      
         })
  }
  else{
    sim.ex<-apply(rnd,2,function(sim){
         m<-fun(spp, sim,...)
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
  }
  preds<-lapply(obs,function(x)x$pred)
  EXs<-sapply(obs,function(x)x$EX)
  eig1<-sapply(obs,function(x)x$EIG1)
  sig<-sapply(EXs,function(E)mean(E<=c(E,sim.ex)))
  res<-list(PCA=PC,preds=preds, MAX=MAX, EX=EXs, eig1=eig1,sim.ex=sim.ex, sig=sig)
  class(res)<-"palaeoSig"
  return(res)
}

