matl<-function(x){
#
# take data in list mode and store it in a matrix
#
J=length(x)
nval=NA
for(j in 1:J)nval[j]=length(x[[j]])
temp<-matrix(NA,ncol=J,nrow=max(nval))
for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
temp
}
