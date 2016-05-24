LBoost <-
function(resp, Xs, anneal.params, nBS=100, kfold=5, nperm=1, PI.imp=NULL, pred.imp=FALSE)
{
  pred<-ncol(Xs)
  if (missing(anneal.params)) {anneal.params<-logreg.anneal.control(start=2, end=-1, iter=50000)}
  nBS2<-nBS/kfold
  if (is.null(colnames(Xs))) {x.nms<-paste("X", 1:pred, sep="")}
  else {x.nms<-colnames(Xs)}
  CVdata<-CV.data(resp=resp, Xs=Xs, kfold=kfold) 
  fullPIdata<-matrix()
  allPInms<-c()
  predmat<-matrix(0, nrow=nBS, ncol=pred+1)
  colnames(predmat)<-c(x.nms, "Null")
  test.resp<-vector("list", kfold)
  test.Xs<-vector("list", kfold)
  CVmod<-vector("list", kfold)
  cvids<-CV.ids(n=length(resp), kfold=kfold)
  wt.mat<-vector("list", kfold)
  alphas<-vector("list", kfold)
  for (k in 1:kfold)
    {
    TRdata<-CVdata[[k]]$traindata
    TRXs<-TRdata[,1:pred]
    TRresp<-TRdata[,(pred+1)]
    TSdata<-CVdata[[k]]$testdata
    test.Xs[[k]]<-TSdata[,1:pred]
    test.resp[[k]]<-TSdata[,(pred+1)]
    n<-nrow(TRdata)    
    iprob<-rep(1/n, n)
    prob.list<-matrix(0, nrow=nBS2, ncol=n)
    prob.list[1,]<-iprob
    alpha<-c()
    fitlist<-vector("list", nBS2)
    PI.TSdata<-vector("list", nBS2)
    PIdata<-matrix()
    for(b in 1:nBS2)
      {
      prob<-prob.list[b,]
      nleaves<-sample(2:8, 1, replace=FALSE) 
      FinalX<-TRXs[,1:pred]
      colnames(FinalX)<-x.nms[1:pred]
      fit <- logreg(resp = TRresp, bin = FinalX, wgt=prob, type=1, select=1,
             ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
      fitlist[[b]]<-fit
      wt.est<-ada.weights(fit=fit, cwts=prob, Xs=FinalX, resp=TRresp)
      if (b<nBS2) {prob.list[b+1,]<-wt.est$nwts}
      alpha<-append(alpha, wt.est$alpha)
      test.ids<-CV.ids(n=nrow(Xs), kfold=kfold)
      PIdata.info<-LBoost.PIs(ids=test.ids[[k]], resp=resp, Xs=Xs, currentfit=fit, pred=pred)
      pred.ids<-PIdata.info$pred.ids
      if (pred.ids[1]=="Null") {predmat[((k-1)*nBS2 +b), pred+1]<-1}
      else {predmat[((k-1)*nBS2 +b), pred.ids]<-1}
      PI.TSdata[[b]]<-PIdata.info$PITSdata#data matrix describing PIs from tree b for test data for this CV run!
      if(b==1){PIdata<-PIdata.info$Alldata}
      if(b>1){PInms<-colnames(PIdata)
         new.PInms<-colnames(PIdata.info$Alldata)
         all.nms<-unique(c(PInms, new.PInms))
         c.ids<-which(new.PInms%in%PInms)
         if(length(c.ids)==0) {PIdata<-cbind(PIdata, PIdata.info$Alldata)}
         if(length(c.ids)>0) {PIdata<-cbind(PIdata, PIdata.info$Alldata[,-c.ids])}
         colnames(PIdata)<-all.nms
         }
      }
    if(k==1){fullPIdata<-PIdata
             allPInms<-append(allPInms, colnames(PIdata))}
    if(k>1){allPInms<-append(allPInms, colnames(PIdata))
        c.PInms<-colnames(fullPIdata)
        n.PInms<-colnames(PIdata)
        c.ids<-which(n.PInms%in%c.PInms)
        if(length(c.ids)==0) {fullPIdata<-cbind(fullPIdata, PIdata)}
        if(length(c.ids)>0) {fullPIdata<-cbind(fullPIdata, PIdata[,-c.ids])}
        }
    mod<-list(AllFits=fitlist, weights=prob.list, nBS=nBS, traindata=TRdata, testdata=TSdata, 
              PI.TSdata=PI.TSdata)
    alphas[[k]]<-alpha
    wt.mat[[k]]<-prob.list
    CVmod[[k]]<-mod
    }
  pred.freq<-colSums(predmat)
  PI.freq<-table(allPInms)
  CVmiss<-CV.err(CVmod=CVmod, alphas=alphas, pred=pred)
  if(pred.imp==TRUE) {pred.import<-Pred.imp(CVmod=CVmod, alphas=alphas, pred=pred, nperm=nperm, CVmiss=CVmiss$CVmiss)}
  else {pred.import<-"Not measured"}
  if(is.null(PI.imp)) {primeimp.import1<-primeimp.import2<-"Not measured"}
  if(PI.imp=="AddRemove") {primeimp.import1<-PlusMinus.PIimp(resp=test.resp, Xs=test.Xs, fullPIdata=fullPIdata, mod=CVmod, wts=alphas, ntrees=nBS)$Imp
                 primeimp.import2<-"Not measured"}
  if(PI.imp=="Permutation") {primeimp.import2<-Perm.PIimp(resp=test.resp, fullPIdata=fullPIdata, mod=CVmod, wts=alphas, ntrees=nBS, nperm=nperm)
                 primeimp.import1<-"Not measured"}
  if(PI.imp=="Both") {primeimp.import1<-PlusMinus.PIimp(resp=test.resp, Xs=test.Xs, fullPIdata=fullPIdata, mod=CVmod, wts=alphas, ntrees=nBS)$Imp
                    primeimp.import2<-Perm.PIimp(resp=test.resp, fullPIdata=fullPIdata, mod=CVmod, wts=alphas, ntrees=nBS, nperm=nperm)}
  ans<-list(CVmod=CVmod, CVmisclass=CVmiss, AddRemove.PIimport=primeimp.import1, Perm.PIimport=primeimp.import2, 
            Pred.import=pred.import, Pred.freq=pred.freq, PI.frequency=PI.freq, wt.mat=wt.mat, alphas=alphas, data=Xs, PIimp=PI.imp, PredImp=pred.imp)
  class(ans)<-"LBoost"
  ans
}
