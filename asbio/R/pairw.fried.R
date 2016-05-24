pairw.fried<-function(y, x, blocks, nblocks, conf = 0.95)
{block.ranks<-matrix(ncol=nlevels(x),nrow=nblocks)
for(i in 1:nblocks){
block.ranks[i,]<-rank(y[blocks==i])
}	
mean.ranks<-apply(block.ranks,2,mean)

R1<-tapply(y,x,length)
r<-length(R1)
fitted<- tapply(y, x, median)
dif.mat<-outer(mean.ranks,mean.ranks,"-")
diffs<-dif.mat[upper.tri(dif.mat)]

SE.diff<-sqrt((r*(r+1))/(6*nblocks))
B<-qnorm(1-((1-conf)/(2*(r^2-r)/2)))
p.val<-2*pnorm(abs(diffs)/SE.diff,lower.tail=FALSE)
p.adj<-ifelse(p.val*((r^2-r)/2)>=1,1,round(p.val*((r^2-r)/2),6))
hwidths<-B*SE.diff
val<-round(cbind(diffs,diffs-hwidths,diffs+hwidths),5)
Decision<-ifelse((val[,2]>0&val[,3]>0)|val[,2]<0&val[,3]<0,"Reject H0","FTR H0")
mr<-as.matrix(mean.ranks)
row.names(mr)<-levels(x)
val<-as.data.frame(cbind(val,Decision,p.adj))
lvl<-outer(levels(x),levels(x),function(x1,x2){paste(paste("Avg.rank",x1,sep=""),paste("Avg.rank",x2,sep=""),sep="-")})
dimnames(val)<-list(lvl[upper.tri(lvl)],
c("Diff","Lower","Upper","Decision","Adj. P-value"))
head<-paste(paste(as.character(conf*100),"%",sep=""),c("confidence intervals for Friedman's comparisons"))
###
res <- list()
res$head <- head
res$conf <- conf
comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
res$comp <- comp[upper.tri(comp)]
res$summary <- val
res$band <- cbind(diffs-hwidths, diffs+hwidths)
res$fitted <- fitted
res$mean.rank.in.blocks<-mr
res$method <- "mBonferroni"
res$x <- x
res$y <- y
class(res)<-"pairw"
res
}
