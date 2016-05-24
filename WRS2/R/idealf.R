idealf<-function(x,na.rm=FALSE){
#
# Compute the ideal fourths for data in x
#
if(na.rm)x<-x[!is.na(x)]
j<-floor(length(x)/4 + 5/12)
y<-sort(x)
g<-(length(x)/4)-j+(5/12)
ql<-(1-g)*y[j]+g*y[j+1]
k<-length(x)-j+1
qu<-(1-g)*y[k]+g*y[k-1]
list(ql=ql,qu=qu)
}
