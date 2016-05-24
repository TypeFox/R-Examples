"HWE2sided.table" <-
function(maf,n,ylim=c(0,1),xlim=NULL){

## n is no. genotypes (not alleles)

na <- maf*n*2
nb <- 2*n-na

if(na%%2==0) nba <- seq(0,na,2)
else nba <- seq(1,na,2)

nbb <- (nb - nba)/2
naa <- (na - nba)/2

xtable <- matrix(rep(0,4*length(nba)),ncol=4)

for(i in 1:length(nba)){
 geno <- c(nba[i],nbb[i],naa[i])
 dum <- HWE2sided(geno)
 xtable[i,1] <- nba[i]
 #xtable[i,2] <- nbb[i]
 #xtable[i,3] <- naa[i]
 xtable[i,2] <- dum$pval.inbreed
 xtable[i,3] <- dum$pval.H
 xtable[i,4] <- dum$pval.cond
}

xtable <- data.frame(xtable)
names(xtable) <- c("no. het","pval.inbred","pval.Hald","pval.cond")
#print(round(xtable,d=4))

#plot(nba,xtable[,3],type="b",pch=19,ylim=ylim)
#lines(nba,xtable[,4],lty=2,type="b",pch=4,lwd=2)

if(is.null(xlim)){
	plot(nba,xtable[,3],type="p",pch="o",ylim=ylim,ylab="p-		value",xlab="no. heterozygotes")
	points(nba,xtable[,4],lty=2,pch="x")
}
else{
	plot(nba,xtable[,3],type="p",pch="o",ylim=ylim,xlim=xlim,
	ylab="p-value",xlab="no. heterozygotes")
	points(nba,xtable[,4],lty=2,pch="x")
}
title(paste("maf = ",maf,",  n = ",n,sep=""))

##return(xtable)
print(round(xtable,d=4))

}

