###########################################################################
## This 'tangle' R script was created from an RSP document.
## RSP source document: 'babel.tex.rsp'
## Metadata 'title': 'The babel vignette'
## Metadata 'keywords': 'example, R, package, vignette'
## Metadata 'author': 'Adam B. Olshen, Richard A. Olshen, and Barry S. Taylor'
## Metadata 'engine': 'R.rsp::rsp'
###########################################################################

R.utils::use("R.utils, R.devices")
options(width=100)
getOption("devEval/args/path", "figures/")
withCapture({
library(babel)
data(ribo.prof)
})
withCapture({
test.rna <- ribo.prof$test.rna
print(test.rna[1:5,])
})
withCapture({
test.rp <- ribo.prof$test.rp
print(test.rp[1:5,])
})
toPDF("scatter", tags="1A", {
   plot(test.rna[,1]+1,test.rp[,1]+1,xlab="mRNA counts",ylab="RPF counts",pch=16,log="xy",xlim=c(1,10000),ylim=c(1,10000),font.lab=2)
   })
withCapture({
test.group <- c("A","B","A","B")
})
withCapture({
options(mc.cores=1)
})
withCapture({
set.seed(12345)
})
withCapture({
test.babel <- babel(test.rna,test.rp,group=test.group,nreps=1e+05,min.rna=10)
})
withCapture({
within.babel <- test.babel$within
print(within.babel[[1]][1:5,])
})
toPDF("scatter2", tags="1A", {
   which.025 <- which(within.babel[[1]]$"P-value (one-sided)"<0.025)
   which.975 <- which(within.babel[[1]]$"P-value (one-sided)">0.975)	
   plot(test.rna[,1]+1,test.rp[,1]+1,xlab="mRNA counts",ylab="RPF counts",pch=16,log="xy",xlim=c(1,10000),ylim=c(1,10000),font.lab=2)
   points(test.rna[which.025,1]+1,test.rp[which.025,1]+1,pch=16,col=2)
   points(test.rna[which.975,1]+1,test.rp[which.975,1]+1,pch=16,col=3)
   })
withCapture({
combined.babel <- test.babel$combined
print(combined.babel[[1]][1:5,])
})
toPDF("scatter3", tags="1A", {
   which.025 <- which(combined.babel[[1]]$"P-value"<0.025)
   plot(test.rna[,1]+1,test.rp[,1]+1,xlab="mRNA counts",ylab="RPF counts",pch=16,log="xy",xlim=c(1,10000),ylim=c(1,10000),font.lab=2)
   points(test.rna[which.025,1]+1,test.rp[which.025,1]+1,pch=16,col=2)
   })
toPDF("scatter4", tags="1A", {
   which.025 <- which(combined.babel[[1]]$"P-value"<0.025)
   plot(test.rna[,3]+1,test.rp[,3]+1,xlab="mRNA counts",ylab="RPF counts",pch=16,log="xy",xlim=c(1,10000),ylim=c(1,10000),font.lab=2)
   points(test.rna[which.025,3]+1,test.rp[which.025,3]+1,pch=16,col=2)
   })
withCapture({
between.babel <- test.babel$between
print(between.babel[[1]][1:5,])
})
toLatex(sessionInfo())
