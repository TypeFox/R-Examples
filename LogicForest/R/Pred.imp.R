Pred.imp <-
function(CVmod, alphas, pred, nperm, CVmiss)
{
 M.CVPpred<-matrix(0, nrow=nperm, ncol=pred)
 for (m in 1:nperm)
   {
   CVPpred<-matrix(0, nrow=length(CVmod), ncol=pred)
   for (i in 1:length(CVmod))
     {
     testdata<-CVmod[[i]]$testdata #matrix of test data for this CV run
     truth<-testdata[,pred+1]
     fits<-CVmod[[i]]$AllFits  #vector of LR trees for this CV run
     c.alphas<-alphas[[i]]  #Vector of weights for fits in this CV run
     perm.ids<-sample(1:nrow(testdata), nrow(testdata), replace=F)
     miss.mat<-matrix(0, nrow=nrow(testdata), ncol=pred)  
     for (j in 1:pred)
       {
       p.testdata<-testdata
       p.testdata[,j]<-testdata[perm.ids,j]
       vote.mat<-matrix(0, nrow=nrow(p.testdata), ncol=length(fits)) 
       for (k in 1:length(fits))
         {
         fit<-fits[[k]]
         prd<-predict.logreg(fit, newbin=as.matrix(p.testdata[,1:pred]))
         pm.prd<-(prd+(prd-1))
         vote.mat[,k]<-c.alphas[k]*pm.prd
         }
       vote.sum<-rowSums(vote.mat)
       vote<-ifelse(vote.sum>0, 1, 0)
       miss.mat[,j]<-abs(vote-truth)
       }
     CVPpred[i,]<-colMeans(miss.mat) #proportion missed for each permuted predictor
     }
   M.CVPpred[m,]<-colMeans(CVPpred)
   }
 Pred.imp<-colMeans(M.CVPpred)-CVmiss
 names(Pred.imp)<-colnames(CVmod[[1]]$testdata[,1:pred])
 Pred.imp
}
