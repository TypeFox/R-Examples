PlusMinus.PIimp <-
function(resp, Xs, fullPIdata, mod, wts, ntrees)
{
 All.PInms<-colnames(fullPIdata)
 n.obs<-nrow(fullPIdata) 
 nCV<-length(mod)
 test.ids<-CV.ids(n=n.obs, kfold=nCV)
 loc.ids<-CV.ids(n=ntrees, kfold=nCV)
 APIchange<-matrix(0, nrow=nCV, ncol=length(All.PInms))
 SPIchange<-matrix(0, nrow=nCV, ncol=length(All.PInms))
 colnames(APIchange)<-All.PInms
 colnames(SPIchange)<-All.PInms
 for(i in 1:nCV)
   {
   ids<-test.ids[[i]]
   chg.ids<-loc.ids[[i]]
   TR.resp<-resp[[i]]
   TR.Xs<-Xs[[i]]
   fits<-mod[[i]]$AllFits
   alphas<-wts[[i]]
   Orig.pred<-ada.pred(fits=fits, alphas=alphas, Xs=TR.Xs)
   Orig.miss<-sum(abs(Orig.pred-TR.resp))
   tPIdata<-fullPIdata[ids,]
   c.PIdat<-mod[[i]]$PI.TSdata
   for(j in 1:length(All.PInms))
     {
     PI<-All.PInms[j]
     Anewpred.mat<-matrix(0, nrow=nrow(TR.Xs), ncol=length(c.PIdat))
     Snewpred.mat<-matrix(0, nrow=nrow(TR.Xs), ncol=length(c.PIdat))
     for(k in 1:length(c.PIdat))
       {
       tree.PIdat<-c.PIdat[[k]]
       if(is.vector(tree.PIdat)){Otree.pred<-alphas[k]*ifelse(tree.PIdat==1, 1, -1)}
       if(is.matrix(tree.PIdat)){Otree.pred<-alphas[k]*ifelse(rowSums(tree.PIdat)>0, 1, -1)}
       tree.PInms<-colnames(tree.PIdat)
       if(PI%in%tree.PInms==TRUE) 
         {
         Anewpred.mat[,k]<-Otree.pred
         loc<-which(tree.PInms%in%PI)
         new.treePIdat<-as.matrix(tree.PIdat[,-loc])
         if(dim(new.treePIdat)[2]==0) {Ntree.pred<-alphas[k]*rep(-1, length(TR.resp))}
         if(dim(new.treePIdat)[2]>0) {Ntree.pred<-alphas[k]*ifelse(rowSums(new.treePIdat)>0, 1, -1)}
         Snewpred.mat[,k]<-Ntree.pred
         }
       if(PI%in%tree.PInms==FALSE)
         {
         Snewpred.mat[,k]<-Otree.pred
         new.treePIdat<-cbind(tree.PIdat, tPIdata[,j])
         Ntree.pred<-alphas[k]*ifelse(rowSums(new.treePIdat)>0, 1, -1)
         Anewpred.mat[,k]<-Ntree.pred
         }
       }
       Apred<-ifelse(rowSums(Anewpred.mat)>0, 1, 0)  
       Spred<-ifelse(rowSums(Snewpred.mat)>0, 1, 0)
       Amiss<-sum(abs(Apred-TR.resp))
       Smiss<-sum(abs(Spred-TR.resp))
       Adiff<-Orig.miss-Amiss
       Sdiff<-Smiss-Orig.miss
       APIchange[i,j]<-Adiff
       SPIchange[i,j]<-Sdiff
     }
   }
 APIimp<-colSums(APIchange)#/nrow(APIchange)
 SPIimp<-colSums(SPIchange)#/nrow(SPIchange)
 Imp<-rowSums(cbind(APIimp, SPIimp))
 ans<-list(APIimp=APIimp, APIchange=APIchange, SPIimp=SPIimp, SPIchange=SPIchange, Imp=Imp)
}
