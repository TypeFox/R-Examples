LBoost.PIs <-
function(ids, resp, Xs, currentfit, pred)
{
 n<-nrow(Xs)
 n2<-length(ids)
 data<-cbind(Xs, resp)
 if (sum(currentfit$model$trees[[1]]$trees[,3])>0)
   {
   PIinfo<-prime.imp(tree=currentfit$model$trees[[1]], data=data, Xs=c(1:pred))
   PIs<-pimp.mat(pimps.out=PIinfo, testdata=data)
   PImat<-PIs$pimp.datamat
   PITSdata<-as.matrix(PImat[ids,])
   colnames(PITSdata)<-colnames(PImat)
   pred.ids<-PIinfo$vec.pimpvars
   }
 if (sum(currentfit$model$trees[[1]]$trees[,3])==0)
   {
   PImat<-as.matrix(rep(1, n))
   PITSdata<-as.matrix(rep(1, n2))
   colnames(PImat)<-colnames(PITSmat)<-"Null"
   pred.ids<-"Null"
   }
 ans<-list(Alldata=PImat, PITSdata=PITSdata, pred.ids=pred.ids) 
}
