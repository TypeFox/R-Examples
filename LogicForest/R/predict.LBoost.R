predict.LBoost <-
function(object, newdata,...)
{ 
 if (class(object)!="LBoost")
   stop("object not of class LBoost")
 if (missing(newdata)) {newdata<-object$data}
 mods<-object$CVmod
 alphas<-object$alphas
 max.vote<-sum(unlist(alphas))
 nCV<-length(alphas)#number of CV runs
 ntrees<-nCV*length(mods[[1]]$AllFits)
 n<-nrow(newdata)
 pred.mat<-matrix(0, nrow=n, ncol=ntrees)
 for (i in 1:nCV)
   {
   fits<-mods[[i]]$AllFits
   alpha<-alphas[[i]]
   for (j in 1:length(fits))
     {
     prd<-predict(fits[[j]], newbin=newdata)
     pred<-ifelse(prd==0, -1, prd)
     pred.mat[,(j+(i-1)*length(fits))]<-alpha[j]*pred
     }
   }
 prd.sign<-rowSums(pred.mat)
 prdct<-ifelse(prd.sign<0, 0, 1)
 wt.vote<-prd.sign/max.vote
 ans<-list(prediction=prdct, weighted.prop=wt.vote)
 class(ans)<-"LBoost.prediction"
 ans
}
