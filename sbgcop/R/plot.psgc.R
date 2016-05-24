"plot.psgc" <-
function(x,...){ 

Y<-x$Y.pmean
p<-dim(Y)[2]
vnames<-colnames(Y)
nc<-1
if(!is.null(vnames)) {nc<-max( nchar(vnames)) }


par(mar=c(.6,2,.5,1),mgp=c(.6,.5,0),oma=c(.75*nc,2.85,1,1),xpd=FALSE,cex=.8)
par(mfcol=c(p,3))
###
for(j in 1:p){
y<-Y[,j]
if(length(unique(y))<=10 ){
barplot(table(y)/sum(!is.na(y)),xlab="",ylab=vnames[j],xaxt="n",yaxt="n")
                           }
if(length(unique(y))>10 ){
hist(y,prob=TRUE,main="",xlab="",ylab=vnames[j],,xaxt="n",yaxt="n",col="gray")
                          }
              }
###
plotci.sA(x$C.psamp,ylabs=rep("",p),mgp=c(.6,.5,0))
plotci.sA(sR.sC(x$C.psamp),ylabs=rep("",p),mgp=c(.6,.5,0))
###
                         }

