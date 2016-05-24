`multi.mat` <-
function (training.spp, envs, core.spp,noanalogues=10, method="sq-chord", run="both"){
  ests<-function(d.mat, nRow, nSamp){
    lapply(noanalogues,function(ana){
      res<-sapply(1:nRow,function(s){
        analogues<-(1:nSamp)[order(d.mat[,s], decreasing=F)][1:ana]
        nDWmean<-colMeans(envs[analogues,,drop=F])
        c(nDWmean=nDWmean)
      })
      t(res)
    })
  }
  if(missing(core.spp)) run="jack"
  nSamp<-nrow(training.spp)
  if(run=="both"||run=="core")    spp<-rbind(training.spp, core.spp)
  else spp<-training.spp
  dist.mat<-make.dist(spp,method=method)
  diag(dist.mat)<-Inf
  if(!run=="core"){
    jack.dist.mat<-dist.mat[1:nSamp,1:nSamp]
    jack<-ests(d.mat=jack.dist.mat, nRow=nrow(training.spp), nSamp=nSamp)
  }
  if(!run=="jack"){
    core.dist.mat<-dist.mat[1:nSamp,-(1:nSamp)]
    core<-ests(d.mat=core.dist.mat, nRow=nrow(core.spp), nSamp=nSamp)
  }
  if(run=="both")list(jack=jack,core=core)
  else if(run=="jack") jack
  else if(run=="core") core
}

