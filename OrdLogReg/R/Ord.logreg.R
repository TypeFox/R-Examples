Ord.logreg <-
function(resp, Xs, nleaf, use.cv=TRUE, kfold=5, anneal.params) 
{
 if(missing(anneal.params)) {anneal.params<-logreg.anneal.control(start=2, end=-1, iter=50000)}
 if(use.cv==FALSE) {CV<-"Cross-validation not used"}
 else {CV<-paste(kfold, "-fold cross-validation", sep="")}
 maxcat<-max(resp)#defining highest category
 nms<-paste("Category", c(1:maxcat), sep=" ")
 model<-vector("list", maxcat)#matrix to collect dummy response variables
 mod.dat<-vector("list", maxcat)#list of each dataset used to construct model
 names(model)<-names(mod.dat)<-nms
 mod.preds<-c()#vector of predictor ids
 pos<-c()#vector of indicators for preds if compliment or not
 leaves<-c()
 n<-nrow(data)
 m<-0
 Ys<-vector("list", maxcat)
 for (i in maxcat:1)#building model from highest category down
   {
   mod.dat[[i]]<-Xs
   Y<-ifelse(resp>=i, 1, 0)#constructing dummy response, includes all cats >= current i
   Ys[[i]]<-Y
   if (use.cv==TRUE) {
     if(missing(nleaf)) {nleaf<-c(1,8)}
     sz.fit<-logreg(resp = Y, bin = Xs, type = 1, select = 3, ntrees = 1, nleaves = nleaf,  
       kfold=5, anneal.control = anneal.params)
     sz.cvresp<-sz.fit$cvscores[c(5,10,15,20,25,30,35,40),8]
     nleaves<-order(sz.cvresp)[1]
     leaves<-append(leaves, nleaves)
     fit <- logreg(resp = Y, bin = Xs, type = 1, select = 1, ntrees = 1, nleaves = nleaves,  
         anneal.control = anneal.params)
     Ps<-fit$model$tree[[1]]$trees[,3]#identifying predictors used in model
     p.loc<-which(Ps!=0)
     x.loc<-Ps[c(p.loc)]
     if (length(p.loc)<=1) {mod.ps<-colnames(Xs)[x.loc]}
     else {mod.ps<-colnames(Xs[,c(x.loc)])}
     mod.preds<-append(mod.preds, mod.ps)
     pos<-append(pos, fit$model$tree[[1]]$trees[c(p.loc),4])
     mod.pred<-predict.logreg(fit)
     correct.ids<-which(mod.pred==1 & Y==1)#identifying which observations predicted correctly
     if (length(correct.ids)==0) {Xs<-Xs; resp<-resp}
     if (length(correct.ids)>0) {Xs<-Xs[-c(correct.ids),]; resp<-resp[-c(correct.ids)]}
     model[[i]]<-fit
     }
   else {
     if(missing(nleaf) & i==maxcat) {nleaf<-rep(8, maxcat); leaves<-append(leaves, nleaf)}
     if(length(nleaf)==maxcat & sum(nleaf)!=24 & i==maxcat) {leaves<-append(leaves, nleaf)}
     if(length(nleaf)==1 & i==maxcat) {nleaf<-rep(nleaf, maxcat); leaves<-append(leaves, nleaf)}
     if(length(nleaf)!=maxcat & i==maxcat) {stop("The number of leaves specified does not equal the number of trees to be constructed")}
     fit <- logreg(resp = Y, bin = Xs, type = 1, select = 1, ntrees = 1, nleaves = nleaf[i], anneal.control = anneal.params)
     Ps<-fit$model$tree[[1]]$trees[,3]#identifying predictors used in model
     p.loc<-which(Ps!=0)
     x.loc<-Ps[c(p.loc)]
     if (length(p.loc)<=1) {mod.ps<-colnames(Xs)[x.loc]}
     else {mod.ps<-colnames(Xs[,c(x.loc)])}
     mod.preds<-append(mod.preds, mod.ps)
     pos<-append(pos, fit$model$tree[[1]]$trees[c(p.loc),4])
     mod.pred<-predict.logreg(fit)
     correct.ids<-which(mod.pred==1 & Y==1)#identifying which observations predicted correctly
     if (length(correct.ids)==0) {Xs<-Xs; resp<-resp}
     if (length(correct.ids)>0) {Xs<-Xs[-c(correct.ids),]; resp<-resp[-c(correct.ids)]}
     model[[i]]<-fit
     }
   }
 ans<-list(mod.dat=mod.dat, model=model, Ys=Ys, mod.preds=mod.preds, pos=pos, leaves=leaves, CV=CV)
 class(ans)<-"Ord.logreg"
 return(ans)
}
