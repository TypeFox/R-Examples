GAIPE.RMSEA <-
function(rmsea,df,width=NULL,clevel=.95,N=c(100,1800,15),
PA_method=c("exact.fit","close.fit","not.close.fit"),H0rmsea=NULL,alpha=.05,
power=c(.8,.9,.95)){;n=N;N=seq(n[1],n[2],n[3]);par(mar=c(4,6,2,5))
plot(NA,ylim=c(min(N),max(N)+100),xlim=c(0,max(max(rmsea)+.1,.2)),
xlab=expression(RMSEA~"("~epsilon~")"),ylab='Sample Size (N)',axes=FALSE,main="",bg='white',xaxs='i',yaxs='i')
u = par("usr");pchwidth=c(15,16,17,18,3,4,8,7,9,10,12,13,18,19,20)
abline(v=c(.05,.08,.1),lty=1,col=colors()[205])
text(expression(epsilon~"="~.05,epsilon~"="~.08,epsilon~"="~.1),x=c(.05,.08,.1),y=rep(u[4]-40,3),
bg='white',col=colors()[205])
arrows(.05,u[4]-80,.05,u[4],col='white',length=0);arrows(.08,u[4]-80,.08,u[4],col='white',length=0)
arrows(.1,u[4]-80,.1,u[4],col='white',length=0)
points(rmsea,rep(u[3],length(rmsea)),pch=23,cex=1.5,lwd=1.5,bg='white')
arrows(u[1], u[3], u[2], u[3],xpd = TRUE,length=.0);arrows(u[1], u[3], u[1], u[4],xpd = TRUE,length=.0)
axis(2,at=seq(min(N),max(N),50),labels=F,tcl=-.3)
axis(2,at=seq(min(N),max(N),200),tcl=-.6)
rmsea.range=round(seq(0,max(max(rmsea)+.1,.2),.01),2)
axis(1,at=rmsea.range,labels=F,tcl=-.3);axis(1,at=seq(min(rmsea.range),max(rmsea.range),.02),tcl=-.6)
if(length(rmsea)==1){;ci=rep(0,2)
for(i in 1:length(N)){;CI=CI.RMSEA(rmsea,df,N[i],clevel)
ci=rbind(ci,c(CI$L,CI$U))};ci=ci[-1,]
for(i in 1:length(N)){lines(ci[i,],rep(N[i],2),lwd=1.8)}
if(!is.null(width)[1]){;for(k in 1:length(width)){
nn=AIPE.RMSEA(rmsea,df,width[k],clevel=clevel)
CI_=CI.RMSEA(rmsea,df,nn,clevel);CI_limit=c(CI_$L,CI_$U)
d_w=abs(CI_$U-CI_$L-width[k])
points(CI_limit,rep(nn,2),pch=pchwidth[k],cex=1,lwd=1.5)
lines(CI_limit,rep(nn,2),lwd=3)
text(CI_limit[2]+.01,nn,paste("N =",nn,""),cex=.8,lwd=1.8,lwd=2)}}}else{
for(j in 1:length(rmsea)){;ci=rep(0,2)
for(i in 1:length(N)){;CI=CI.RMSEA(rmsea[j],df,N[i],clevel)
ci=rbind(ci,c(CI$L,CI$U))};ci=ci[-1,]
lines(ci[,1],N,lty=j,lwd=2);lines(ci[,2],N,lty=j,lwd=2)
if(!is.null(width)[1]){;for(k in 1:length(width)){
nn=AIPE.RMSEA(rmsea[j],df,width[k],clevel=clevel)
CI_=CI.RMSEA(rmsea[j],df,nn,clevel);CI_limit=c(CI_$L,CI_$U)
points(CI_limit,rep(nn,2),pch=pchwidth[k],cex=1.25,lwd=2)
lines(CI_limit,rep(nn,2),lwd=2)
text(mean(CI_limit),nn-50,paste("N =",nn,""),cex=.8,lwd=1.8)}}
exp_rmsea_lab=c()
for(i in 1:length(rmsea))
exp_rmsea_lab=c(exp_rmsea_lab,substitute(Expected~epsilon==x,list(x=rmsea[i])))
legend("bottomright",legend=as.expression(exp_rmsea_lab),
lty=1:length(rmsea),lwd=2,cex=1,bg='white')}};if(!is.null(H0rmsea)){
if(!length(PA_method)==1|!is.element(PA_method,
c("exact.fit","close.fit","not.close.fit")))
stop("PA_method has to be correctly specified.")
HArmsea=rmsea[1];NN=seq(n[1],n[2],200)
NP=numeric(length(power));POWER=power;P=numeric(length(NN))
if(PA_method=="close.fit"|PA_method=="exact.fit"){;if(is.na(H0rmsea)|H0rmsea<0)
stop("H0rmsea has to be correctly specified.")
if(HArmsea<=H0rmsea&PA_method=="close.fit")
stop("For the test of close fit, H0rmsea has to be smaller than expected RMSEA.")
if(H0rmsea!=0&PA_method=="exact.fit")
stop("For the test of exact fit, H0rmsea has to be 0.")
for(k in 1:length(NN)){P[k]=1-pchisq(qchisq(1-alpha,df,H0rmsea^2*df*(NN[k]-1)),
df,HArmsea^2*df*(NN[k]-1))}
for(l in 1:length(power)){NP[l]=PA.RMSEA(df,"close.fit",H0rmsea,HArmsea,POWER[l],
alpha)}}
if(PA_method=="not.close.fit"){;if(is.na(H0rmsea)|H0rmsea<0)
stop("H0rmsea has to be correctly specified.");if(HArmsea>=H0rmsea)
stop("For the test of not-close fit, H0rmsea has to be larger than expected RMSEA.")
for(k in 1:length(NN)){P[k]=pchisq(qchisq(alpha,df,H0rmsea^2*df*(NN[k]-1)),df,
HArmsea^2*df*(NN[k]-1))}
for(l in 1:length(power)){NP[l]=PA.RMSEA(df,"not.close.fit",H0rmsea,HArmsea,
POWER[l],alpha)}}
arrows(u[2], u[3], u[2], u[4],xpd = TRUE,length=.0)
axis(4,at=seq(min(N),max(N),50),labels=F,tcl=-.3)
axis(4,at=seq(min(N),max(N),200),labels=round(P,2),tcl=-.6)
if(H0rmsea>=HArmsea)
mtext(substitute(Power~"("~H[0]*":"~epsilon>x~")",list(x=H0rmsea)),side=4,line=3)
if(H0rmsea<HArmsea&H0rmsea!=0)
mtext(substitute(Power~"("~H[0]*":"~epsilon<=x~")",list(x=H0rmsea)),side=4,line=3)
if(H0rmsea==0)
mtext(substitute(Power~"("~H[0]*":"~epsilon==x~")",list(x=H0rmsea)),side=4,line=3)
nr=ifelse(length(rmsea)==1,0,length(rmsea))
abline(h=NP,lty=1:length(power)+nr,lwd=2);if(!is.null(width)[1]){
legend("topright",legend=c(paste("Width =",width,""),paste("Power =",power,"")),
lty=c(rep(1,length(width)),1:length(power)+nr),cex=1,lwd=2,pch=c(pchwidth[1:length(width)],
rep(NA,3)),bg='white');NPP=sort(NP);for(a in 1:length(NPP))
text(x=max(max(rmsea)+.1,.2)-a*.02,y=NPP[a]+25,paste("N =",NPP[a],""),cex=.8)
};if(is.null(width)[1]){
legend("topright",legend=paste("Power =",power,""),lty=1:length(power),cex=1,
lwd=2,bg='white');NPP=sort(NP);for(a in 1:length(NPP))
text(x=max(max(rmsea)+.1,.2)-a*.02,y=NPP[a]+25,paste("N =",NPP[a],""),cex=.8)}
if(length(rmsea)==1){
title(substitute(list(df==x,Expected~epsilon==y, epsilon~"("~H[0]~")"~"="~z, epsilon~"("~H[a]~")"~"="~b),list(x=df,y=rmsea,z=H0rmsea,b=HArmsea)))}
if(length(rmsea)!=1){
title(substitute(list(df==x, epsilon~"("~H[0]~")"~"="~z, epsilon~"("~H[a]~")"~"="~b),list(x=df,z=H0rmsea,b=HArmsea)))}
}else{;if(!is.null(width)[1]){
legend("topright",legend=paste("Width =",width,""),lty=rep(1,length(width)),cex=1,lwd=2,pch=pchwidth[1:length(width)],bg='white')}
if(length(rmsea)==1){title(substitute(list(df==x,Expected~epsilon==y),list(x=df,y=rmsea)))}
if(length(rmsea)!=1){title(substitute(list(df==x),list(x=df)))}}}
