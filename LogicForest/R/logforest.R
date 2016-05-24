logforest <-
function(resp, Xs, nBSXVars, anneal.params, nBS=100, h=0.5, norm=TRUE, numout=5)
{
  pred<-ncol(Xs)
  if (missing(anneal.params)) {anneal.params<-logreg.anneal.control(start=2, end=-1, iter=50000)}
  if (missing(nBSXVars)) {nBSXVars<-pred}
  n<-nrow(Xs)
  if (is.null(colnames(Xs))) {x.nms<-paste("X", 1:pred, sep="")}
  else {x.nms<-colnames(Xs)}
  fitlist<-vector("list", nBS)
  IBdata<-vector("list", nBS)
  OOBdata<-vector("list", nBS)
  OOBpred<-matrix(nrow=n, ncol=nBS)
  single.predimport<-vector("list",nBS)
  vimp.import<-vector("list", nBS)
  treepreds.list<-vector("list", nBS)
  for(b in 1:nBS)
     { 
     nleaves<-sample(2:8, 1, replace=FALSE) 
     BSindices<-sample(1:n, n, replace=TRUE) 
     OOBindices<-(1:n)[!is.element(1:n, BSindices)]
     BS<-Xs[BSindices, ]
     OOB<-Xs[OOBindices, ]
     BSY<-resp[BSindices] 
     OOBY<-resp[OOBindices]
     XVarIndices<-sort(sample(1:pred, nBSXVars, replace=FALSE)) 
     rsBSX<-BS[ ,XVarIndices]
     rsOOBX<-OOB[,XVarIndices]
     FinalBS<-cbind(rsBSX, BSY) 
     FinalOOB<-cbind(rsOOBX, OOBY)
     colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
     fit <- logreg(resp = BSY, bin = FinalBS[,1:nBSXVars], 
            type = 1, select = 1, ntrees = 1, nleaves = nleaves,   anneal.control = anneal.params)
     OOBpred[as.vector(OOBindices),b]<-predict.logreg(fit, newbin=as.matrix(FinalOOB[,1:nBSXVars]))
     if (sum(fit$model$tree[[1]]$trees[,3])!=0) {
        pred.import<-pimp.import(fit=fit, data=FinalBS, testdata=FinalOOB, BSpred=length(XVarIndices), 
        pred=pred, Xs=XVarIndices) 
        vimp.import[[b]]<-pred.import$pimp.vimp
        single.predimport[[b]]<-pred.import$single.vimp
        treepreds.list[[b]]<-pred.import$Xids
        }
     else {single.predimport[[b]]= 0 ;vimp.import[[b]]= 0; treepreds.list[[b]]=0}
     fitlist[[b]]<-fit
     IBdata[[b]]<-BSindices
     OOBdata[[b]]<-OOBindices
     }
  OOB.pred<-matrix(0, nrow=n, ncol=2)
  for (i in 1:n)
     {
     pred.ids<-which(OOBpred[i,]==1|OOBpred[i,]==0)
     pred.vec<-OOBpred[i,c(pred.ids)]
     OOBprop<-sum(pred.vec)/length(pred.vec)
     OOBi.pred<-ifelse(OOBprop>h, 1, 0)
     OOB.pred[i,]<-c(OOBi.pred, OOBprop)
     }
  colnames(OOB.pred)<-c("predicted_class","proportion")
  OOBmisclass<-sum(abs(OOB.pred[,1]-resp))/n
  pred.importmat<-matrix(0, nrow=nBS, ncol=pred)
  colnames(pred.importmat)<-x.nms
  for (i in 1:nBS)
    {
    pred.ind<-treepreds.list[[i]]
    m<-length(pred.ind)
    for (j in 1:m)
      {
      col<-pred.ind[j]
      pred.importmat[i,col]<-single.predimport[[i]][j]
      }
    }
  pred.imp<-colSums(pred.importmat)
  names(pred.imp)<-x.nms
  freq.table<-table(names(unlist(single.predimport)))
  all.pimps<-unique(names(unlist(vimp.import)))
  npimps<-length(unique(names(unlist(vimp.import))))
  pimp.freqtable<-table(names(unlist(vimp.import)))
  pimptable<-matrix(0, nrow=nBS, ncol=npimps)
  colnames(pimptable)<-all.pimps
  for (i in 1:nBS)
    {
    npimps.tree<-length(vimp.import[[i]])
    col.ids<-which(colnames(pimptable)%in%names(vimp.import[[i]]))
    n.ids<-length(col.ids)
    for (j in 1:n.ids)
      {
      pimptable[i,col.ids[j]]<- vimp.import[[i]][j]
      }
    }
  pimpsum<-colSums(pimptable)
  t5PIs<-names(sort(pimpsum, decreasing=TRUE)[1:5]) 
  ans<-list(AllFits=fitlist, Top5.PI=t5PIs, Predictor.importance=pred.imp, Predictor.frequency=freq.table,
            PI.frequency=pimp.freqtable, PI.importance=pimpsum,  ModelPI.import=vimp.import,
            OOBmisclass=OOBmisclass, OOBprediction=OOB.pred, IBdata=IBdata, OOBdata=OOBdata, norm=norm, 
            numout=numout, predictors=pred, Xs=Xs)
 class(ans)<-"logforest"
 ans
}
