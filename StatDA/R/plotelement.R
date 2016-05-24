"plotelement" <-
function(da.object){

# Element plot for Discriminant Analysis

# number of groups
ngroups=nrow(da.object$means)
# number of variables
p=ncol(da.object$means)
# variable names
varnam=dimnames(da.object$means)[[2]]

par(mfrow=c(ngroups,1),mar=c(2,4,4,2))

lim=max(1,0.2+max(abs(da.object$means)))

for (i in 1:ngroups){
  plot(1:p,da.object$means[i,],ylim=c(-lim,lim),type="n",xaxt="n",xlab="",
       ylab="Group centre",cex.lab=1.2)
  title(paste("Group",i),cex=1.5)
  text(1:p,da.object$means[i,],labels=varnam,col=1,cex=1.2)
  abline(h=0)
  abline(h=c(-0.5,0.5),lty=2)
} 

}





