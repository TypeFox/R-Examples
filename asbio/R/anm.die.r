
anm.die<-function(reps=300,interval=0.1,show.die=TRUE,p=c(1/6,1/6,1/6,1/6,1/6,1/6),cl=TRUE){
old.par <- par(no.readonly = TRUE)
if(sum(p)!=1)stop("sum of die probabilities must = 1")
p<-p*100
die<-sample(c(rep(1,p[1]),rep(2,p[2]),rep(3,p[3]),rep(4,p[4]),rep(5,p[5]),rep(6,p[6])),reps,replace=TRUE)
p1<-seq(1:reps);p2<-seq(1:reps);p3<-seq(1:reps);p4<-seq(1:reps);p5<-seq(1:reps);p6<-seq(1:reps)
for(i in 1:reps){
p1[i]<-length(die[0:i][die[0:i]==1])/i;p2[i]<-length(die[0:i][die[0:i]==2])/i
p3[i]<-length(die[0:i][die[0:i]==3])/i;p4[i]<-length(die[0:i][die[0:i]==4])/i
p5[i]<-length(die[0:i][die[0:i]==5])/i;p6[i]<-length(die[0:i][die[0:i]==6])/i}

for(i in 1:reps){
dev.hold()
if(show.die==TRUE){
layout(matrix(c(rep(1,6),0,2,0), 3, 3, byrow = TRUE)) 
par(mar=c(5.5,4.5,2,2.1))
exp = 1.3
}

else exp <- 1.2

plot(seq(1,reps),seq(0,reps-1)/(reps-1),type="n",xlab="Trials",ylab="Cumulative proportion", cex.lab = exp + .1, cex.axis = exp + .1)
points(seq(1:i),p1[1:i],lty=ifelse(cl==TRUE,1, 1),col=ifelse(cl==TRUE,1, 1),type="l", lwd=1.3)
points(seq(1:i),p2[1:i],lty=ifelse(cl==TRUE,1, 2),col=ifelse(cl==TRUE,2, gray(.2)),type="l",lwd=exp)
points(seq(1:i),p3[1:i],lty=ifelse(cl==TRUE,1, 3),col=ifelse(cl==TRUE,3, gray(.4)),type="l",lwd=exp)
points(seq(1:i),p4[1:i],lty=ifelse(cl==TRUE,1, 4),col=ifelse(cl==TRUE,4, gray(.6)),type="l",lwd=exp)
points(seq(1:i),p5[1:i],lty=ifelse(cl==TRUE,1, 1),col=ifelse(cl==TRUE,5, gray(.8)),type="l",lwd=exp)
points(seq(1:i),p6[1:i],lty=ifelse(cl==TRUE,1, 2),col=ifelse(cl==TRUE,4, gray(.6)),type="l",lwd=exp)

if(cl == TRUE){ 
lcol <- 1:6
llty <- rep(1, 6)
}
else{
llty <- c(1,2,3,4,1,2)
lcol <- gray(c(0,.2,.4,.6,.8,.6))
}

legend("topright",col=lcol, lty=llty, lwd=exp, legend=c("1","2","3","4","5","6"), title="Die result", bg="white", cex = exp-.1)
if(show.die==TRUE){
par(mar=c(0.2,0,0,0))
if(die[i]==1)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(0.5,0.5,pch=19,cex=8)}
if(die[i]==2)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(c(0.1,0.9),c(0.1,0.9),pch=19,cex=8)}
if(die[i]==3)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(c(0.1,0.5,0.9),c(0.1,0.5,0.9),pch=19,cex=8)}
if(die[i]==4)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(c(0.1,0.1,0.9,0.9),c(0.1,0.9,0.1,0.9),pch=19,cex=8)}
if(die[i]==5)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(c(0.1,0.1,0.5,0.9,0.9),c(0.1,0.9,0.5,0.1,0.9),pch=19,cex=8)}
if(die[i]==6)
{plot(seq(0,1),seq(0,1),xlab="",ylab="",xaxt="n",yaxt="n",type="n");points(c(0.1,0.5,0.9,0.1,0.5,0.9),c(0.1,0.1,0.1,0.9,0.9,0.9),pch=19,cex=8)}}
dev.flush()
Sys.sleep(interval)
}
on.exit(par(old.par))
invisible()	
}