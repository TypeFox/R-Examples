FindLambda.MA <-
function(x,n,Z=1.96^2,correct=FALSE,a=c(-1,1),plot=FALSE){
p=x/n; li=sum(a[a<0]); ls=a%*%p

if (correct==TRUE){ N=exp( sum(log(n+1)) ); cpc= sum(abs(a))/ (2*(N-1))
ls=ls-cpc}

ls=ls-0.001*(ls-li) 
if (li==ls) {li=li-0.0001; ls=ls+0.0001} 

if (score.MA(x,n,Z,correct=FALSE,a=a,lambda=(li+ls)/2)==0) {raiz1=(li+ls)/2}
else { raiz1=uniroot(score.MA,c(li,ls),x=x,n=n,Z=Z,correct=correct,a=a)$root}

#if (plot){ par(mfcol=c(1,2))
#xgrid=seq(li,ls,length=100);fgrid=xgrid
#for (i in 1:100) fgrid[i]=score.MA(x,n,Z=Z,correct,a,lambda=xgrid[i])
#if (correct) {titulo="corrected score"} else {titulo="non corrected score"}
#plot(xgrid,fgrid,type="l",xlab="z",ylab="score",main=titulo,
#col="blue",lwd=2); abline(v=raiz1,col="red");abline(h=0)}

li=a%*%p; ls=sum(a[a>0]); if (correct==TRUE){li=li+cpc}

li=li+0.001*(ls-li) 

if (li==ls) {li=li-0.0001; ls=ls+0.0001}  

if (score.MA(x,n,Z,correct=FALSE,a=a,lambda=(li+ls)/2)==0) {raiz2=(li+ls)/2}
else {raiz2=uniroot(score.MA,c(li,ls),x=x,n=n,Z=Z,correct=correct,a=a)$root }

#if (plot){ xgrid=seq(li,ls,length=100);fgrid=xgrid
#for (i in 1:100) fgrid[i]=score.MA(x,n,Z=Z,correct,a,lambda=xgrid[i])
#if (correct) {titulo="corrected score"} else {titulo="non corrected score"}
#plot(xgrid,fgrid,type="l",xlab="z",ylab="score",main=titulo,
#col="blue",lwd=2); abline(v=raiz2,col="red");abline(h=0)}
c(raiz1,raiz2)}
