##tests with continuous, binary, and survival response, with and
##without parallelization.

library(pensim)
data(beer.exprs)
data(beer.survival)

gene.quant <- apply(beer.exprs,1,quantile,probs=0.75)
dat.filt <- beer.exprs[gene.quant>log2(100),]
gene.iqr <- apply(dat.filt,1,IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr>0.5,])
dat.filt <- t(dat.filt)
dat.filt <- dat.filt[,1:50]

library(survival)
surv.obj <- Surv(beer.survival$os,beer.survival$status)
cont.obj <- beer.survival$os
bin.obj <- factor(beer.survival$status)

##
##
## tests without parallelization
##
##


##--------------------------------
## survival
##--------------------------------
set.seed(1)
opt.nested.crossval(outerfold=5,nprocessors=1,  #opt.nested.crossval arguments
             optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
             setpen="L1",nsim=1,                         #opt1D arguments
             response=surv.obj,
             penalized=dat.filt,
             fold=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


##--------------------------------
## continuous
##--------------------------------
set.seed(1)
opt.nested.crossval(outerfold=5,nprocessors=1,  #opt.nested.crossval arguments
             optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
             setpen="L1",nsim=1,                         #opt1D arguments
             response=cont.obj,
             penalized=dat.filt,
             fold=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


##--------------------------------
## binary
##--------------------------------
set.seed(1)
opt.nested.crossval(outerfold=5,nprocessors=1,  #opt.nested.crossval arguments
             optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
             setpen="L1",nsim=1,                         #opt1D arguments
             response=bin.obj,
             penalized=dat.filt,
             fold=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


##
##
## tests with parallelization
##
##

if(require(parallel)){

  ##--------------------------------
  ## survival
  ##--------------------------------
  set.seed(1)
  opt.nested.crossval(outerfold=5,nprocessors=2,  #opt.nested.crossval arguments
                      optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
                      setpen="L1",nsim=2,                         #opt1D arguments
                      response=surv.obj,
                      penalized=dat.filt,
                      fold=5,
                      positive=FALSE,
                      standardize=TRUE,
                      trace=FALSE)


  ##--------------------------------
  ## continuous
  ##--------------------------------
  set.seed(1)
  opt.nested.crossval(outerfold=5,nprocessors=2,  #opt.nested.crossval arguments
                      optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
                      setpen="L1",nsim=2,                         #opt1D arguments
                      response=cont.obj,
                      penalized=dat.filt,
                      fold=5,
                      positive=FALSE,
                      standardize=TRUE,
                      trace=FALSE)


  ##--------------------------------
  ## binary
  ##--------------------------------
  set.seed(1)
  opt.nested.crossval(outerfold=5,nprocessors=2,  #opt.nested.crossval arguments
                      optFUN="opt1D",scaling=FALSE,               #opt.splitval arguments
                      setpen="L1",nsim=2,                         #opt1D arguments
                      response=bin.obj,
                      penalized=dat.filt,
                      fold=5,
                      positive=FALSE,
                      standardize=TRUE,
                      trace=FALSE)
}
