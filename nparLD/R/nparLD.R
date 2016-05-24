nparLD<-function(formula,data=NULL,subject,
description=TRUE,time1.order=NULL,time2.order=NULL,group1.order=NULL, group2.order=NULL,
plot.CI=FALSE,alpha=0.05,show.covariance=FALSE,order.warning=TRUE){

#-----------------Store the Input Argument------------------------------------#
input.list<-list(formula=formula,data=data,subject=subject,description=description,
time1.order=time1.order,time2.order=time2.order,group1.order=group1.order,
group2.order=group2.order,plot.CI=plot.CI,alpha=alpha,show.covariance=show.covariance,
order.warning=order.warning)

#-----------------Transform the subject parameter into a suitable form--------#
if(is.character(subject))
{
 subject<-data[,subject]
}

if(is.numeric(subject))
{
 if(length(subject)==1)
 {
  subject<-data[,subject]
 }
}

#-----------------Create a Model Framae-----------------------------------------#
dat <- model.frame(formula, data)
dat2 <- data.frame(dat,subject=subject)
dat2<-dat2[order(dat2$subject),]
nc <- ncol(dat2)
nadat2 <- names(dat2)

#---------------------LD-F1-Design---------------------------------------------#

if (nc==3){
a1<-ld.f1(dat2[,1], dat2[,2], subject=dat2[,3], time.name=nadat2[2],
description=description, time.order=time1.order,plot.RTE=FALSE,
show.covariance=show.covariance,order.warning=order.warning)
a2<-ld.ci(dat2[,1], dat2[,2], dat2[,3], alpha=alpha, time.name=nadat2[2], group.name="Group",
description=FALSE, time.order=time1.order, group.order=NULL,order.warning=FALSE)
a1$Conf.Int<-a2
a1$input<-input.list

if (plot.CI==TRUE){
t1<-factor(dat2[,2])
nt1<-nlevels(t1)
uu<-"-"
Lower<-a1$Conf.Int$Lower
Upper<-a1$Conf.Int$Upper
RTE<-a1$Conf.Int$RTE
plot(rep(1,2),c(Lower[1],Upper[1]),pch=uu,cex=3,lwd=2,type="o",xlim=c(1,nt1),
ylim=c(0,1),xaxt="n",ylab="",xlab=nadat2[2],cex.lab=1.5)
for (ii in 2:nt1){
points(rep(ii,2),c(Lower[ii],Upper[ii]),pch=uu,cex=3,lwd=2,type="o")}
points(c(1:nt1),RTE,pch=15,type="o",lwd=2)
axis(1,at=1:nt1, labels=a1$Conf.Int$Time)
title(main=paste("Relative Effects"))}
class(a1)<-"nparLD"
return(a1)}
#--------------------LD-F2-Design----------------------------------------------#
if (nc==4) {
t11 <- factor(dat2[,2])
t22 <- factor(dat2[,3])
subf <- factor(dat2[,4])
nt11<-nlevels(t11)
nt22<-nlevels(t22)
hldf2 <-length(dat[,1])/(nt11*nt22)
if (hldf2 == nlevels(subf)) {
a1<-ld.f2(dat2[,1], dat2[,2], dat2[,3], subject=dat2[,4], description=description,time1.name=nadat2[2],time2.name=nadat2[3], time1.order=time1.order,
time2.order=time2.order,plot.RTE=FALSE,show.covariance=show.covariance,order.warning=order.warning)
ht2<-c(paste(dat2[,2],dat2[,3]))
dat3<-data.frame(dat2,TG=ht2)
a2<-ld.ci(dat3[,1],dat3[,5],dat2[,4],alpha=alpha,description=FALSE,order.warning=FALSE)
a2<-a2[order(a2$Time),]
a1$Conf.Int<-data.frame(a2,Time1=sort(rep(levels(t11),nt22)),Time2=rep(levels(t22),nt11))
a1$input<-input.list
if (plot.CI==TRUE){
uu<-"-"
samples<-split(a1$Conf.Int,a1$Conf.Int$Time2)
plot(1:nt11,samples[[1]]$RTE,pch=15,type="o",ylim=c(0,1),xlim=c(0,nt11+1),xaxt="n",ylab="",xlab=nadat2[2],cex.lab=1.5,lwd=2)
for (hh in 0:(nt22-1)){
for (ss in 1:nt11){
points(rep(ss+hh/10,2), c(samples[[hh+1]]$Lower[ss],samples[[hh+1]]$Upper[ss]),pch=uu,cex=3,lwd=2,type="o",col=hh+1)
points((1:nt11)+hh/10,samples[[hh+1]]$RTE,pch=15,lwd=2,type="o",col=hh+1)}}
axis(1,at=1:nt11, labels=levels(t11))
title(main=paste("Relative Effects"))
legend("top", col=c(1:nt22),paste(nadat2[3],levels(t22)),pch=c(rep(10,nt22)),lwd=c(rep(2,nt22)) ) }
class(a1)<-"nparLD"
return(a1)}
#------------------------F1-LD-F1-Design---------------------------------------#
if(hldf2 != nlevels(subf)){
h1 <- which(dat2[,4]==dat2[,4][1])
h2 <- which(dat2[,2]==dat2[,2][1])
h3 <- which(dat2[,3]==dat2[,3][1])
sh1h2 <- sum(h1==h2[h1])
sh1h3 <- sum(h1==h3[h1])
pgf1=ptf1=0
if (sh1h3==length(h1)){
pgf1<-3
ptf1<-2}
if (sh1h2==length(h1)){
pgf1<-2
ptf1<-3}
#------------------------------------------------------------------------------#
a1<-f1.ld.f1(dat2[,1], dat2[,ptf1], dat2[,pgf1], dat2[,4], time.name=nadat2[ptf1],
group.name=nadat2[pgf1], description=description,
time.order=time1.order, group.order=group1.order,plot.RTE=FALSE,show.covariance=show.covariance,order.warning=order.warning)
a2<-ld.ci(dat2[,1], dat2[,ptf1], dat2[,4], group=dat2[,pgf1], alpha=0.05,
time.name=nadat2[ptf1],
group.name=nadat2[pgf1], description=FALSE, time.order=NULL,order.warning=FALSE)
a2<-a2[order(a2$Group),]
a1$Conf.Int<-a2
a1$input<-input.list
#------------------Implementiere die Grafik--------------------#
if(plot.CI==TRUE){
uu<-"-"
samples<-split(a1$Conf.Int,a1$Conf.Int$Group)
t11<-factor(dat2[,ptf1])
t22<-factor(dat2[,pgf1])
nt11<-nlevels(t11)
nt22<-nlevels(t22)
plot(1:nt11,samples[[1]]$RTE,pch=15,type="o",ylim=c(0,1),xlim=c(0,nt11+1),xaxt="n",ylab="",xlab=nadat2[ptf1],cex.lab=1.5,lwd=2)
for (hh in 0:(nt22-1)){
for (ss in 1:nt11){
points(rep(ss+hh/10,2), c(samples[[hh+1]]$Lower[ss],samples[[hh+1]]$Upper[ss]),pch=uu,cex=3,lwd=2,type="o",col=hh+1)
points((1:nt11)+hh/10,samples[[hh+1]]$RTE,pch=15,lwd=2,type="o",col=hh+1)}}
axis(1,at=1:nt11, labels=a1$Conf.Int$Time[1:nt11])
title(main=paste("Relative Effects"))
legend("top", col=c(1:nt22),paste(nadat2[pgf1],levels(t22)),pch=c(rep(10,nt22)),lwd=c(rep(2,nt22)))}
#---------------------------End of the Graphic Code----------------------------#
class(a1)<-"nparLD"
return(a1)}}
#--------------------------------F1-LD-F2--------------------------------------#
if (nc==5){
h11 <- which(dat2[,5]==dat2[,5][1])
lh11 <- length(h11)
h22 <- c(which(dat2[,2] == dat2[,2][1]))[h11]
h33 <- c(which(dat2[,3] == dat2[,3][1]))[h11]
h44 <- c(which(dat2[,4] == dat2[,4][1]))[h11]
sh11h22 <- sum(h11==h22)
sh11h33 <- sum(h11==h33)
sh11h44 <- sum(h11==h44)
#------------------------Group is the second column----------------------------#
ht2=pg=pt1=pt2=pg1=pg2=0
if (sh11h22==lh11 && sh11h33 != lh11 && sh11h44!=lh11 ){
pg<-2
pt1<-3
pt2<-4
ht2<-c(paste(dat2[,3],dat2[,4]))
model<-"f1ldf2"}
#-----------------------Group is the third column------------------------------#
if (sh11h33==lh11 && sh11h22 != lh11 && sh11h44!=lh11 ){
pg<-3
pt1<-2
pt2<-4
ht2<-c(paste(dat2[,2],dat2[,4]))
model<-"f1ldf2"}
#----------------------Group is the fourth column------------------------------#
if (sh11h44==lh11 && sh11h33 != lh11 && sh11h22!=lh11){
pg<-4
pt1<-2
pt2<-3
ht2<-c(paste(dat2[,2],dat2[,3]))
model<-"f1ldf2"}
#--------------------------F2-LD-F1--------------------------------------------#
#----------------------Time is the second column-------------------------------#
if (sh11h22 != lh11 && sh11h33 == lh11 && sh11h44==lh11 ){
pt1<-2
pg1<-3
pg2<-4
ht2<-c(paste(dat2[,3],dat2[,4]))
model<-"f2ldf1"}
#---------------------Time is the third column---------------------------------#
if (sh11h33 != lh11 && sh11h22 == lh11 && sh11h44==lh11 ){
pt1<-3
pg1<-2
pg2<-4
ht2<-c(paste(dat2[,2],dat2[,4]))
model<-"f2ldf1"}
#-------------------Time is the fourth column----------------------------------#
if (sh11h44 != lh11 && sh11h22 == lh11 && sh11h33==lh11 ){
pt1<-4
pg1<-2
pg2<-3
ht2<-c(paste(dat2[,2],dat2[,3]))
model<-"f2ldf1"}
if (model=="f1ldf2"){
dat3<-data.frame(dat2,TG=ht2)
a1<-f1.ld.f2(dat2[,1], dat2[,pt1], dat2[,pt2], dat2[,pg], subject=dat2[,5],
group.name=nadat2[pg],
time1.name=nadat2[pt1],time2.name=nadat2[pt2],description=description,
time1.order=time1.order,
time2.order=time2.order, group.order=group1.order,plot.RTE=FALSE,show.covariance=show.covariance,
order.warning=order.warning)
a2<-ld.ci(dat3[,1],dat3[,6],dat3[,5],group=dat3[,pg],alpha=0.05,
time.name=nadat2[3],
group.name=nadat2[pg], description=FALSE, time.order=levels(factor(dat3$TG)),
group.order=levels(factor(dat3[,pg])),order.warning=FALSE)
a1$Conf.Int<-a2
a1$input<-input.list
#----------------------Include the Graphic to the Output-----------------------#
if (plot.CI==TRUE){
uu<-"-"
g1<-factor(dat3[,pg])
t11<-factor(dat2[,pt1])
t22<-factor(dat2[,pt2])
ng1<-nlevels(g1)
nt11<-nlevels(t11)
nt22<-nlevels(t22)
a1$Conf.Int<-data.frame(a2,Time1=sort(rep(levels(t11),nt22)),Time2=rep(levels(t22),nt11))
samples<-split(a1$Conf.Int,a1$Conf.Int$Group)
par(mfrow=c(1,ng1))
for (hh in 0:(ng1-1)){
datenhilf1<-samples[[hh+1]]
datenhilf<-split(datenhilf1,datenhilf1$Time1)
plot(1:nt22,rep(0.5,nt22),col="white",xlim=c(0,nt22+1), ylim=c(0,1),xaxt="n",ylab="",xlab=nadat2[pt2],cex.lab=1.5)
title(main=paste(nadat2[pg],levels(g1)[hh+1]))
for (jj in 0:(nt11-1)){
for (ss in 1:nt22){
points(rep(ss+jj/10,2), c(datenhilf[[jj+1]]$Lower[ss],datenhilf[[jj+1]]$Upper[ss]),pch=uu,cex=3,lwd=2,type="o",col=jj+1)
points((1:nt22)+jj/10,datenhilf[[jj+1]]$RTE,pch=15,lwd=2,type="o",col=jj+1)}}
axis(1,at=1:nt22, labels=levels(t22))
legend("top", col=c(1:nt11),paste(nadat2[pt1],levels(t11)),pch=c(rep(10,nt11)),lwd=c(rep(2,nt11)))}}
#------------------------End of the graphic implementation---------------------#
class(a1)<-"nparLD"
return(a1)}
if (model=="f2ldf1"){
a1<-f2.ld.f1(dat2[,1], dat2[,pt1], dat2[,pg1], dat2[,pg2], subject=dat2[,5], time.name=nadat2[pt1],
group1.name=nadat2[pg1], group2.name=nadat2[pg2], description=description,
time.order=time1.order, group1.order=group1.order, group2.order=group2.order,plot.RTE=FALSE,
show.covariance=show.covariance,order.warning=order.warning)
dat3<-data.frame(dat2,TG=ht2)
a2<-ld.ci(dat3[,1],dat3[,pt1],dat3[,5],group=dat3[,6],alpha=0.05,
time.name=nadat2[pt1],
group.name=nadat2[pg1], description=FALSE, time.order=levels(factor(dat3[,pt1])),
group.order=levels(factor(dat3[,6])),order.warning=FALSE)
a1$Conf.Int<-a2
a1$input<-input.list
if (plot.CI==TRUE){
uu<-"-"
g1<-factor(dat3[,pg1])
g2<-factor(dat3[,pg2])
t11<-factor(dat2[,pt1])
ng1<-nlevels(g1)
ng2<-nlevels(g2)
nt11<-nlevels(t11)
par(mfrow=c(1,ng1))
a1$Conf.Int<-data.frame(a2,Group1=sort(rep(levels(g1),nt11*ng2)),Group2=sort(rep(levels(g2),nt11)))
a1$input<-input.list
samples<-split(a1$Conf.Int,a1$Conf.Int$Group1)
for (hh in 0:(ng1-1)){
datenhilf1<-samples[[hh+1]]
datenhilf<-split(datenhilf1,datenhilf1$Group2)
plot(1:nt11,rep(0.5,nt11),col="white",xlim=c(0,nt11+1), ylim=c(0,1),xaxt="n",ylab="",xlab=nadat2[pt1],cex.lab=1.5)
title(main=paste(nadat2[pg1],levels(g1)[hh+1]))
for (jj in 0:(ng2-1)){
for (ss in 1:nt11){
points(rep(ss+jj/10,2), c(datenhilf[[jj+1]]$Lower[ss],datenhilf[[jj+1]]$Upper[ss]),pch=uu,cex=3,lwd=2,type="o",col=jj+1)
points((1:nt11)+jj/10,datenhilf[[jj+1]]$RTE,pch=15,lwd=2,type="o",col=jj+1)}}
axis(1,at=1:nt11, labels=levels(t11))
legend("top", col=c(1:ng2),paste(nadat2[pg2],levels(g2)),pch=c(rep(10,ng2)),lwd=c(rep(2,ng2)) )}}
#----------------------End of the graphic implementation-----------------------#
class(a1)<-"nparLD"
return(a1)}}
#-------------------------End-of-Function--------------------------------------#
}
