ada.pred <-
function(fits, alphas, Xs)
#fits is list of boosted LR trees
#alphas is vector of alphas for each tree
{
 pred.mat<-matrix(0, nrow=nrow(Xs), ncol=length(fits))
 for (i in 1:length(fits))
   {
   prd<-predict(fits[[i]], newbin=Xs)
   pred<-ifelse(prd==0, -1, prd)
   pred.mat[,i]<-alphas[i]*pred
   }
 prd.sign<-rowSums(pred.mat)
 prdct<-ifelse(prd.sign<0, 0, 1)
 prdct
}
