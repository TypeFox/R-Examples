superpc.predictionplot <- function (train.obj, data, data.test,  threshold, n.components=3,  n.class=2, shrinkage=NULL, call.win.metafile=FALSE)
{
  



BIG=1000
  
this.call=match.call()
  
if(n.components>3){
     stop("n.components cannot be bigger than 3")
}

if(is.null(shrinkage)){ fit.cts=superpc.predict(train.obj, data, data.test, threshold=threshold, n.components=n.components, prediction.type="continuous")


  pred.cts=superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred[,1:n.components], print=FALSE)$results 


pred.cts.1df=superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred.1df, print=FALSE)$results

                               }
else  {

     fit.cts=superpc.predict.red(train.obj, data, data.test, threshold=threshold, n.components=n.components,  shrinkages=shrinkages)


pred.cts=superpc.fit.to.outcome(train.obj, data.test,fit.cts$v.test[,,1:n.components], print=FALSE)$results


pred.cts.1df=superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.test.1df[1,], print=FALSE)$results

   }



if(call.win.metafile){win.metafile()}


if(train.obj$type=="survival"){

if(n.components==1){layout(matrix(c(1,2),1,2, byrow = TRUE),widths=c(.80,.20),heights=c(1,1))}

if(n.components==2){layout(matrix(c(1,4,2,5,3,6),3,2, byrow = TRUE),widths=c(.80,.20),heights=c(.34,.33,.33))}

if(n.components==3){layout(matrix(c(1,5,2,6,3,7,4,8),4,2, byrow = TRUE),widths=c(.80,.20),heights=rep(.25,4))}


if(is.null(shrinkage)){
  fit.groups<- superpc.predict(train.obj, data, data.test, threshold=threshold, n.components=n.components, prediction.type="discrete", n.class=n.class)


pred.groups=superpc.fit.to.outcome(train.obj,  data.test, fit.groups$v.pred, print=FALSE)$results


pred.groups.1df=superpc.fit.to.outcome(train.obj,  data.test, fit.groups$v.pred.1df, print=FALSE)$results


}
 else{
   
 fit.groups<- superpc.predict.red(train.obj, data, data.test, threshold=threshold, n.components=n.components, shrinkages=shrinkages, prediction.type="discrete", n.class=n.class)


pred.groups=superpc.fit.to.outcome(train.obj,  data.test,fit.groups$v.test[,,1], print=FALSE)$results

pred.groups.1df=superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.test.1df[1,], print=FALSE)$results


}

#plot survival curves

lastc=2+(n.class-1)
par(mar=c(5,4,2,1))
xmax=max(data.test$y)*1.5

for(i in 1:n.components){
if(is.null(shrinkage)){

plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred[,i]), col=2:lastc, xlab="time", ylab="Prob survival", xlim=c(0,xmax))
}
else{
 plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.test[1,,i]), col=2:lastc, xlab="time", ylab="Prob survival", xlim=c(0,xmax))
}


if(i==1 & n.class==2){legend(1.05*max(data.test$y), .8,lty=c(1,1),col=2:lastc,c("low score","high score"), cex=.8)}

if(i==1 & n.class==3){legend(1.05*max(data.test$y), .8, lty=c(1,1,1),col=2:lastc,c("low score","med score", "high score"), cex=.8)}


title(main=paste("Component",as.character(i),sep=" "))
}

# if  number of  components >1, plot combined predictor curves

if(n.components>1){
  if(is.null(shrinkage)){
plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred.1df), col=2:lastc, xlab="time", ylab="Prob survival" , xlim=c(0,xmax))
}
else{
 plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.test.1df[1,]), col=2:lastc, xlab="time", ylab="Prob survival" , xlim=c(0,xmax))
}


  title(main=" Combined predictor")
}

# output results on right ride of plot

res=NULL
rownames=NULL

for( ii in 1:n.components){
if(is.null(shrinkage)){
junk=superpc.fit.to.outcome(train.obj,  data.test, fit.groups$v.pred[,ii], print=FALSE)$results
}
else{
 junk=superpc.fit.to.outcome(train.obj,  data.test,fit.groups$v.test[,,ii], print=FALSE)$results
}

lr=round(2*(junk$loglik[2]-junk$loglik[1]),4)
likrat=c(lr, n.class-1,round(1-pchisq(lr,df=n.class-1),5))
res=c(res,likrat)
rownames=c(rownames,"LR stat", "df", "pvalue")
}





lr=round(2*(pred.groups.1df$loglik[2]-pred.groups.1df$loglik[1]),4)
likrat=c(lr, n.class-1,round(1-pchisq(lr,df=n.class-1),5))
res=c(res,likrat)
rownames=c(rownames,"LR stat", "df", "pvalue")



par(mar=c(1,1,1,0))
par(cex=.8)
leftcol= 0
midcol=.5
nrows=length(res)
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")
xc=c(midcol)
if(n.components==1){ yinc=.03}
if(n.components==2){ yinc=.07}
if(n.components==3){ yinc=.11}
yc=1-(1:nrows)*yinc

# write out first block of results
for(j in 1:3){
 text(leftcol,yc[j],rownames[j], pos=4, col=4)
}
#text(midcol,1, "Linear", pos=4, col=4)
for(j in 1:3){
  text(xc[1],yc[j],labels=as.character(res[j]),  pos=4)
}

if(n.components>1){
# write out 2nd block of results
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")  
 for(j in 4:6){
 text(leftcol,yc[j-3],rownames[j], pos=4, col=4)
}
for(j in 4:6){
  jj=j-3
  text(xc[1],yc[jj],labels=as.character(res[j]),  pos=4)
}}

 if(n.components>2){
# write out 3rd block of results
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")  
 for(j in 7:9){
 text(leftcol,yc[j-6],rownames[j], pos=4, col=4)
}
for(j in 7:9){
  jj=j-6
  text(xc[1],yc[jj],labels=as.character(res[j]),  pos=4)
}}

if(n.components>1){
  # write out combined 1df predictor results
#  for(i in 1:(n.components-1)){
#  plot(0,0,type="n",axes=F,xlab="",ylab="")
#}
 plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="") 
  
  nrows=3
  yc=1-(1:nrows)*yinc
  for(j in 1:nrows){
jj=length(res)-3+j
    text(leftcol,yc[j],rownames[jj], pos=4, col=4)
  }
  for(j in 1:nrows){
jj=length(res)-3+j
  text(xc[1],yc[j],labels=as.character(res[jj]),  pos=4)
}
}
}

if(train.obj$type=="regression"){
layout(matrix(c(1,2),1,2, byrow = TRUE),widths=c(.6,.40),heights=c(1,1))

par(mar=c(5,4,2,1))

plot(data.test$y, fit.cts$v.pred.1df, xlab="outcome",ylab="predicted outcome")


res=round(t(summary(pred.cts)$coef),3)

rownames=c("coef", "se","T stat","pvalue")

par(mar=c(5,1,2,1))

par(cex=.7)
leftcol= 0
midcol=.21
midcol2= .42
rightcol=.63
rightcol2=.84
nrows=nrow(res)
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")
xc=c(midcol,midcol2,rightcol, rightcol2)
yinc=.05
yc=1-(1:nrows)*yinc

for(j in 1:nrows){
 text(leftcol,yc[j],rownames[j], pos=4, col=4)
}
   
text(midcol,1, "Intcpt", pos=4, col=4)

if(n.components==1){
text(midcol2,1, "Comp", pos=4, col=4)
}
if(n.components==2){
text(midcol2,1, "Comp1", pos=4, col=4)
text(rightcol,1, "Comp2", pos=4, col=4)
}
if(n.components==3){
text(midcol2,1, "Comp1", pos=4, col=4)
text(rightcol,1, "Comp2", pos=4, col=4)
text(rightcol2,1, "Comp3", pos=4, col=4)
}

for(j in 1:nrows){
for(i in 1:ncol(res)){
  text(xc[i],yc[j],labels=as.character(res[j,i]),  pos=4)
}}
ylast=yc[nrows]

 junk=summary(pred.cts.1df)$fstat

fstat=junk[1]
pvalue=1-pf(fstat,junk[2],junk[3])
 

text(leftcol,ylast-2*yinc, "F-stat",pos=4,col=4)
text(midcol2,ylast-2*yinc, round(fstat,5),pos=4)
text(leftcol,ylast-3*yinc, "pvalue",pos=4, col=4)
text(midcol2,ylast-3*yinc, round(pvalue,5),pos=4)
par(cex=1)
}


if(call.win.metafile){dev.off()}


}


