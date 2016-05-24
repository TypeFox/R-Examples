
 

spautocor <- function(gen.m,eucl.m, shuffle=FALSE, bins = 10)
{

gd <- gen.m
ed <- eucl.m


if (shuffle==TRUE)
{

gdd <- as.dist(gd)
 gdsample <- sample(1:length(gdd), length(gdd))
gd[lower.tri(gd)] <- gdd[gdsample]
gd[upper.tri(gd)] <- gdd[gdsample]
diag(gd) <- 0
}




 cdmat <- function(gd)
 {
 dimen <- nrow(gd)
 sgd <- sum(gd, na.rm=TRUE)
 cscd <- matrix(colSums(gd, na.rm=TRUE), dimen, dimen) 
 rscd <- matrix(rowSums(gd, na.rm=TRUE), dimen, dimen, byrow=TRUE) 
 cd <- 0.5*( -gd + 1/dimen*( cscd + rscd) - 1/dimen^2*(sgd ))
 cd 
 }
 cd <- cdmat(gd)
#remove upper triangel to speed things up....
ed[upper.tri(ed)] <-NA
diag(ed) <- NA
r<- NA
distance <- NA
N<- NA

steps <- signif(diff(range(ed, na.rm=TRUE))/bins,4)
for (d in 1:bins )
{
index <- which(ed<=d*steps & ed >(d-1)*steps, arr.ind=TRUE)
cx <- sum(cd[index])
cxii<-sum(diag(cd)[index[,1]])
cxjj<-sum(diag(cd)[index[,2]])
r[d] <-  2 * cx /(cxii+cxjj)

distance[d] <- steps*d
N[d] <- length(index)
}
if (shuffle==FALSE) res <- data.frame(bin = distance,  N=N, r =r)
else res <- data.frame(r=r)

res

}





#b<- redpossums[1:100]
#b@other$xy <- b@other$latlong
#

#xy <- read.csv("D:\\Bernd\\Projects\\aprasia\\apfinal\\apxy.csv")
#
#aprasia@other$xy <- xy
#
#gen.m<-as.matrix(gd_smouse(cats, verbose=FALSE))
#eucl.m <- as.matrix(dist(cats@other$xy))
#reps=1000
#bins=10
#
#splist<- spautocor(gen.m, eucl.m, bins=20)
#
#
#system.time(
#bssplist <- replicate(reps, spautocor(gen.m, eucl.m,shuffle=TRUE, bins=bins))
#)
#
#bs <-matrix(unlist(bssplist), nrow=reps, ncol=bins, byrow=TRUE)
#
#bs.l <- apply(bs,2, quantile, probs=0.025, na.rm=TRUE)
#bs.u <- apply(bs,2, quantile, probs=0.975, na.rm=TRUE)
#
#
#
#matplot(cbind(splist$r,bs.u, bs.l), type="l", lty=c(1,2,2), lwd=c(2,1,1), ylab="Spatial autocorrelation r", axes=FALSE, col=c(1,3,3), xlab="distance")
#axis(2)
#axis(1, at=1:nrow(splist), labels=signif(splist$bin,3))
#axis(1, at=1:nrow(splist), labels=splist$N, line=1, tick=FALSE)
#box()
#mtext("N=",1,line=2, at=0)
#mtext("Bins",1,line=1, at=0)
#
