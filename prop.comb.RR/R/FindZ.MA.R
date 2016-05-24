FindZ.MA <-
function(x,n,correct=FALSE,a=c(-1,1),lambda=0,plot=FALSE){
p=x/n; resto=a%*%p-lambda

resto=round(resto,4)

if (resto==0) {raiz=0} 
else
{
# limits
li=4*(resto)^2/sum(a^2/n)
T=sum(x[a>0])+sum((n-x)[a<0])
if(resto>0) {ls=T*(resto)/(lambda-sum ( a[a<0]) )} else {
ls=(sum(n)-T)*resto/(lambda-sum (a[a>0]))}

if (li==ls) {li=li-0.0001; ls=ls+0.0001}

#correction case:
N=exp( sum(log(n+1)) ); cpc= sum(abs(a))/ (2*(N-1)); factor=resto/(abs(resto)-cpc)
if (correct==TRUE) {li=li/factor^2; ls=ls/factor^2}

li=round(li,4); ls=round(ls,4)
               
raiz=uniroot(score.MA,c(li,ls),x=x,n=n,correct=correct,a=a,lambda=lambda)$root
#if (plot){
#xgrid=seq(li,ls,length=100);fgrid=xgrid
#for (i in 1:100) fgrid[i]=score.MA(x,n,Z=xgrid[i],correct,a,lambda)
#if (correct) {titulo="corrected score"} else {titulo="non corrected score"}
#plot(xgrid,fgrid,type="l",xlab="z",ylab="score",main=titulo,
#col="blue",lwd=2); abline(v=raiz,col="red");abline(h=0)
#}
}
raiz}
