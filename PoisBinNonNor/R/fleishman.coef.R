fleishman.coef <-
function(n.C, skewness.vec=NULL, kurtosis.vec=NULL){

fleishman.poly<- function(dd,gamma1,gamma2) {
 r <- c(NA,length(dd))
 r[1] <- dd[1]+dd[3]
 r[2] <- (dd[2]^2)+(6*dd[2]*dd[4])+(2*(dd[3]^2))+(15*(dd[4]^2))-1
 r[3] <- 2*dd[3]*((dd[2]^2)+(24*dd[2]*dd[4])+(105*(dd[4]^2))+2)-gamma1
 r[4] <- 24*((dd[2]*dd[4])+(dd[3]^2)*(1+(dd[2]^2)+(28*dd[2]*dd[4]))+(dd[4]^2)*(12+(48*dd[2]*dd[4])+(141*(dd[3]^2))+(225*(dd[4]^2))))-gamma2
 r
}#myfunc

scm=rbind(0,1,skewness.vec, kurtosis.vec)

coef.mat=matrix(0,4,n.C)
for ( i in 1:n.C){
p0 <- matrix(rnorm(25*4), 25, 4) 
log <- capture.output({
ans <- multiStart(par=p0, fn=fleishman.poly, gamma1=scm[3,i], gamma2=scm[4,i], control=list(trace=FALSE), quiet=FALSE) 
})
if(any(ans$conv==TRUE)==FALSE) stop(cat("The algorithm did not converge for continuous variable",i,"!","\n"))
soln1 <- ans$par[which(ans$converged==TRUE),]
soln <- round(ans$par[which(ans$converged==TRUE),],5)
pats <- do.call(paste, c(as.data.frame(soln), sep = "\r"))
pats <- factor(pats, levels = unique(pats))
amat=cbind(unique(soln), Freq = as.vector(table(pats))) 
index=which(amat[,"Freq"]==max(amat[,"Freq"]))
if(length(index)==1) {coef.mat[,i]<-round(soln1[which(apply(round(soln1,5), 1, function(x) all(x==amat[index,1:4])) )[1],],7)} else
{coef.mat[,i]<-round(soln1[which(apply(round(soln1,5), 1, function(x) all(x==amat[sample(index,1),1:4])) )[1],],7)
}
}#for
rownames(coef.mat)=c("a","b","c","d")
colnames(coef.mat)=paste("N",seq(1:n.C), sep="")
return(coef.mat)
}
