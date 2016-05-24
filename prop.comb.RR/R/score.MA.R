score.MA <-
function(x,n,Z,correct=FALSE,a=c(-1,1),lambda=0) {
b=1-2*x/n; p=x/n; L=a%*%p

rest=L-lambda  #S es cero hay que salir con solucion igual a cero
rest=round(rest,4)
if(rest==0) return(0)

if (correct==FALSE) {
r2=n^2*rest^2+ a^2*Z^2 +2*n*a*b*rest*Z
r2=round(r2,8)
sum(n)*rest+(sum(a)-2*lambda)*Z-sign(rest)*sum(sqrt(r2))}
else{
N=exp( sum(log(n+1)) ); cpc= sum(abs(a))/ (2*(N-1))
factor=rest/(abs(rest)-cpc)
r2=n^2*rest^2+a^2*Z^2*factor^4+2*n*a*b*rest*Z*factor^2
r2=round(r2,8)
sum(n)*rest+(sum(a)-2*lambda)*Z*factor^2-sign(rest)*sum(sqrt(r2))}}
