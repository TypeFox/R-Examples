#' Create the design matrix, variance-covariance matrix, the variance of each
#' pairwise comparison and the efficicency of each pairwise comparison for a
#' cross-over design
#' 
#' Function to read in a cross-over design and create the design matrix X, the
#' variance-covariance matrix of the parameter estimates (X'X)^-1, the variance
#' of each pairwise comparison and the efficicency of each pairwise comparison.
#' 
#' This is done for a model with fixed subject effects, period effects and
#' treatments and for the above model with first-order carry-over effects
#' added.
#' 
#' See the vignette of this package for further details.
#' 
#' @param design Cross-over design.
#' @return list(xmat.no.subjects=xmat.no.subjects,
#' var.trt.pair=var.trt.pair,eff.trt.pair=eff.trt.pair,av.eff.trt.pair=av.eff.trt.pair,
#' var.trt.pair.adj=var.trt.pair.adj,eff.trt.pair.adj=eff.trt.pair.adj,av.eff.trt.pair.adj=av.eff.trt.pair.adj)
#' @author Byron Jones, 18/11/2011
#' @references Jones, B., & Kenward, M. G. (2003). Design and analysis of
#' cross-over trials (Vol. 98). Chapman & Hall.
#' @keywords misc
#' @examples
#' 
#' data(fletcher)
#' design.efficiency(fletcher1)
#' 
#' @export design.efficiency
design.efficiency.old <- function(design) {
  if (class(design)=="CrossoverSearchResult") design <- design@design
  if (class(design)=="CrossoverDesign") {
    #model <- design@model
    design <- design@design
  }
  #model <- getModelNr(model)
  design <- t(design)
  
  ## input
  ## nseq   is the number of sequences
  ## nper   is the number of periods
  ## ntrt   is the number of treatments
  ## nrep   is the nseq x 1 vector containing the number of subjects allocated to each sequence
  ## design is an nseq x nper array containing the nseq treatment sequences
  
  ## create the factors for groups, subjects, periods, treatments and carry-over effects
  
  nper <- dim(design)[2]
  nseq <- dim(design)[1]
  ntrt <- length(levels(as.factor(design)))	
  n<-rep(1,nseq)
  
  ## group factor
  group<-rep(1:nseq,nper*n)
  ##subject factor
  subject<-rep(1:sum(n), rep(nper,sum(n)))
  ##period factor
  per<-rep(1:nper,  sum(n))
  
  
  trt<-NULL
  for(i in 1:nseq){
    trt<-c(trt,rep(design[i,],n[i]))
  }
  ## treatment design matrix
  trt.mat<-matrix(trt,sum(n),nper,byrow=T)
  
  ## carry-over design matrix
  car.mat<-matrix(0,sum(n),nper,byrow=T)
  car.mat[,1]<-1
  car.mat[,2:nper]<-trt.mat[,1:(nper-1)]
  
  car<-NULL
  for(i in 1:sum(n)){
    car<-c(car,car.mat[i,])
  }
  
  
  ## create factors for design matrix
  nsubj<-sum(n)
  ndata<-nsubj*nper
  
  mean01<-matrix(1,ndata,1)
  subj01<-matrix(0,ndata,nsubj)
  per01<-matrix(0,ndata,nper)
  trt01<-matrix(0,ndata,ntrt)
  car01<-matrix(0,ndata,ntrt)
  ## subjects
  for(i in 1:nsubj){
    subj01[,i][subject==i]<-1
  }
  ##periods
  for(i in 1:nper){
    per01[,i][per==i]<-1
  }
  ## treatments
  for(i in 1:ntrt){
    trt01[,i][trt==i]<-1
  }
  ##carry-over
  car01[,1:ntrt][per==1]<-0
  for(j in 2:nper){
    car01[,1:ntrt][per==j]<-trt01[,1:ntrt][per==(j-1)]
  }
  
  ## X matrix without subjects - returned for use in simulating data - not used further here
  xmat.no.subjects<-matrix(cbind(mean01,per01,trt01,car01),ndata,(1+nper+2*ntrt),byrow=F)
  ## X matrix
  xmat<-matrix(cbind(mean01,subj01,per01,trt01,car01),ndata,(1+nsubj+nper+2*ntrt),byrow=F)
  ## X'X matrix
  xtx<-t(xmat)%*%xmat
  ## inverse of X'X matrix
  xtx.inv<-ginv(xtx)
  
  ## extract part of inverse for carry-over effects 
  n.cols.all<-1+nsubj+nper+2*ntrt
  n.cols.car.upp<-n.cols.all
  n.cols.car.low<-n.cols.car.upp-(ntrt-1)
  n.cols.trt.upp<-n.cols.car.low-1
  n.cols.trt.low<-n.cols.trt.upp-(ntrt-1)
  ## extract part of inverse for treatments adjusted for carry-over effects 
  xtx.inv.car.adj<-xtx.inv[n.cols.car.low:n.cols.car.upp,n.cols.car.low:n.cols.car.upp]
  xtx.inv.trt.adj<-xtx.inv[n.cols.trt.low:n.cols.trt.upp,n.cols.trt.low:n.cols.trt.upp]
  
  ## repeat the above for the model without carry-over effects
  ## X matrix
  xmat<-matrix(cbind(mean01,subj01,per01,trt01),ndata,(1+nsubj+nper+ntrt),byrow=F)
  ## X'X matrix
  xtx<-t(xmat)%*%xmat
  ## inverse of X'X matrix
  xtx.inv<-ginv(xtx)
  ## extract part of inverse for carry-over effects 
  n.cols.all<-1+nsubj+nper+ntrt
  n.cols.trt.upp<-n.cols.all
  n.cols.trt.low<-n.cols.trt.upp-(ntrt-1)
  ## extract part of inverse for treatments unadjusted for carry-over effects 
  xtx.inv.trt<-xtx.inv[n.cols.trt.low:n.cols.trt.upp,n.cols.trt.low:n.cols.trt.upp]
  
  ## variances of pairwise comparisons
  ## adjusted for carry-over effects
  var.trt.pair.adj<-matrix(0,ntrt,ntrt)
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      if (i!=j) {
        var.trt.pair.adj[i,j]<-xtx.inv.trt.adj[i,i]+xtx.inv.trt.adj[j,j] - 2*xtx.inv.trt.adj[i,j]
      }}}
  ## variances of pairwise comparisons
  ## unadjusted for carry-over effects
  var.trt.pair<-matrix(0,ntrt,ntrt)
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      if (i!=j) {
        var.trt.pair[i,j]<-xtx.inv.trt[i,i]+xtx.inv.trt[j,j] - 2*xtx.inv.trt[i,j]
      }}}
  
  ## efficiency calculations
  ## replication of each treatment in the design
  rep<-matrix(0,ntrt)
  for(i in 1:ntrt){
    rep[i]<-sum(trt.mat==i)
  }
  ## variances of pairwise comparisons in ideal design
  var.trt.pair.ideal<-matrix(0,ntrt,ntrt)
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      if (i!=j) {
        var.trt.pair.ideal[i,j]<-1/rep[i] + 1/rep[j]
      }}}
  ## matrix of efficiencies for pairwise treatment comparisons unadjusted for carry-over effects
  eff.trt.pair<-matrix(0,ntrt,ntrt)
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      if (i!=j) {
        eff.trt.pair[i,j]<-var.trt.pair.ideal[i,j]/var.trt.pair[i,j]
      }}}
  ## matrix of efficiencies for pairwise treatment comparisons adjusted for carry-over effects
  eff.trt.pair.adj<-matrix(0,ntrt,ntrt)
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      if (i!=j) {
        eff.trt.pair.adj[i,j]<-var.trt.pair.ideal[i,j]/var.trt.pair.adj[i,j]
      }}}
  ## average efficiency of pairwise treatment comparisons unadjusted for carry-over effects
  av.eff.trt.pair<-sum(eff.trt.pair)/(ntrt*(ntrt-1))
  ## average efficiency of pairwise treatment comparisons adjusted for carry-over effects
  av.eff.trt.pair.adj<-sum(eff.trt.pair.adj)/(ntrt*(ntrt-1))
  
  return(list(xmat.no.subjects=xmat.no.subjects,
              var.trt.pair=var.trt.pair,eff.trt.pair=eff.trt.pair,av.eff.trt.pair=av.eff.trt.pair,
              var.trt.pair.adj=var.trt.pair.adj,eff.trt.pair.adj=eff.trt.pair.adj,av.eff.trt.pair.adj=av.eff.trt.pair.adj))
  
}

test.eff <- function() {
  f <- stop
  # Test old versus new function:
  path <- system.file("data", package="Crossover")
  for (file in dir(path=path)) {     
    if (file %in% c("clatworthy1.rda", "clatworthyC.rda", "pbib2combine.rda", "exampleSearchResults2t.rda")) next
    designs <- load(paste(path, file, sep="/"))
    for (designS in designs) {
      design <- get(designS)
      r1 <- design.efficiency(design)
      r2 <- design.efficiency.old(design)
      if(!isTRUE(all.equal(r1$var.trt.pair.adj, r2$var.trt.pair.adj))) {
        f(paste("Unequal variances for",designS," (",max(abs(r1$var.trt.pair.adj - r2$var.trt.pair.adj)),")!\n"))
      }
      if(!isTRUE(all.equal(r1$eff.trt.pair.adj, r2$eff.trt.pair.adj))) {
        f(paste("Unequal efficiencies for",designS," (",max(abs(r1$eff.trt.pair.adj - r2$eff.trt.pair.adj))," - ",getCounts(design),")!\n"))
      }
      if(!isTRUE(all.equal(r1$xmat, r2$xmat.no.subjects, check.names=FALSE))) {
        #f(paste("Unequal design matrices for",designS," (",max(abs(r1$eff.trt.pair.adj - r2$eff.trt.pair.adj))," - ",getCounts(design),")!\n"))
      }
    }
  }
}
