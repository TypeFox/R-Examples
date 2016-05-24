TwoSide.varyEffect <-
function(s1,s2,m,m1,delta,a1,r1,fdr){

h<-function(n)
{
a2=1-a1
alpha_star=r1*fdr/((m-m1)*(1-fdr))
beta_star=1-r1/m1
h=0
for(k in 1:length(delta))
{
h=h+(1-pnorm(qnorm(1-alpha_star/2)-abs(delta[k])*sqrt(n*a1*a2)))
}
h=h-r1
return(h)
}

if (h(s1) < 0 & h(s2)>0)
{
while(h(s1) < 0 & h(s2)>0 & (abs(s1-s2)>1 | abs(h(s2))>1)) 
{ 
s3=(s1+s2)/2; 
if(h(s3)>0) s2=s3 else s1=s3 
}
}

return(list("s1"=s1,"s2"=s2,"h(s1)"=h(s1),"h(s2)"=h(s2)))
}
