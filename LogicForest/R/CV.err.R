CV.err <-
function(CVmod, alphas, pred)
{
  vote<-c()
  truth<-c()
  CVpred<-vector("list", length(CVmod))
  for (i in 1:length(CVmod))
   {
   testdata<-CVmod[[i]]$testdata #matrix of test data for this CV run
   truth<-append(truth, testdata[,pred+1])
   fits<-CVmod[[i]]$AllFits  #vector of LR trees for this CV run
   c.alphas<-alphas[[i]]  #Vector of weights for fits in this CV run
   vote.mat<-matrix(0, nrow=nrow(testdata), ncol=length(fits))  #matrix to collect "votes" for this CV run
   for (j in 1:length(fits))
    {
    fit<-fits[[j]]
    prd<-predict.logreg(fit, newbin=as.matrix(testdata[,1:pred]))
    pm.prd<-(prd+(prd-1))
    vote.mat[,j]<-c.alphas[j]*pm.prd
    }
   vote.sum<-rowSums(vote.mat)
   vote<-append(vote, ifelse(vote.sum>0, 1, 0))
   CVpred[[i]]<-ifelse(vote.sum>0, 1, 0)
   }
 CVmiss<-sum(abs(vote-truth))/length(truth)
 ans<-list(CVmiss=CVmiss, CVpred=CVpred)
 ans
}
