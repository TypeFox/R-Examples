
aml.pred.outside<-function(marker, response, kin, which.pred, numkeep, selectvar){
## in this function, adaptive mixed lasso with a subset of data (training set of lines), the result is then used
## to calculate genetic values (prediction) for the remaining lines (validation set of lines)
## marker is the genetic effect matrix for all lines (both for training and validation)
## y is the response (phenotype) of all lines 
## kin is the relationship matrix for all lines, again both for training and validation sets.
## which.pred is a vector to designate which lines should be used for validation,
##      it corresponding to the row numbers in the genetic effect matrix.
## The lines not designated by which.pred is used as the training set
## nkeep(for numkeep) and sv (for selectvar) are passed to amltest()
## The value is a list of two items
##     pred.vl is the vector of predicted values (genetic values) for lines in the validation sample
##     y.vl is the actual observed phenotype values in the validation sample.


 try<- response[- which.pred]
 trmarker<- marker[- which.pred,]
 trmark<- sweep(trmarker,2, apply(trmarker,2,mean),"-")
 trkin<- kin[- which.pred, -which.pred]
 vly<- response[which.pred]
 vlmarker<- marker[which.pred,]
 vlmarker<- sweep(vlmarker,2, apply(vlmarker,2,mean),"-")
 fit<- amltest(try, trmarker, trkin,numkeep=numkeep, selectvar=selectvar )
 selvlmarker<-vlmarker[,as.vector(fit$estimate[,1])]
 seltrmarker<-trmarker[,as.vector(fit$estimate[,1])]
 fixpredvl<- mean(try)+as.matrix(selvlmarker)%*% as.vector(fit$estimate[,2])
 fixpredtr<- mean(try)+as.matrix(seltrmarker)%*% as.vector(fit$estimate[,2])
 vars<-fit$vars
 sigma<- vars[1]*kin+ vars[1]*vars[2]*diag(rep(1,length(response)))
 sigma11<- sigma[-which.pred, -which.pred]
 sigma21<- sigma[which.pred, -which.pred]

 predvl <- fixpredvl+sigma21 %*% solve(sigma11) %*% (try-fixpredtr)

 res<-list(pred.vl=predvl, response.vl=vly )
 return(res)

}


