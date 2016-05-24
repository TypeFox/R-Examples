mycor_jgee1 <-
function(N,nr,nt,y,X,family,beta_new,corstr1,Mv,corstr2,maxclsz,R1,R2,scale.fix,scale.value) {

eta=X%*%beta_new
mu=family$linkinv(eta)
sd=sqrt(family$variance(mu))
res<-(as.vector(y)-mu)/sd 

if(scale.fix==0) {
fi<-sum(res^2)/(nr*sum(nt))
} else
if(scale.fix==1) {
fi<-scale.value
}

aindex=nr*cumsum(nt)
index=c(0,aindex[-length(aindex)])

if (corstr1=="independence") 
{alfa_hat<-0} else 
if (corstr1=="exchangeable") 
{
sum1<-0
sum3<-0
for ( i in  1:N)         {
for ( r in 1:nr)         {
for ( j in 1:nt[i])      {    
for ( jj in 1:nt[i])     {    
if  ( j!=jj)             {
#cat("i",i,"r",r,"j",j,"jj",jj,"\n")
sum2<-res[j+index[i]+nt[i]*(r-1)]*res[jj+index[i]+nt[i]*(r-1)]
sum1<-sum1+sum2
#cat("i",i,"r",r,"j",j,"jj",jj,"sum2",sum2,"sum1",sum1,"\n")
}
}
}
}
sum4<-nt[i]*(nt[i]-1)
sum3<-sum3+sum4
} #i
alfa_hat<-sum1/(nr*sum3*fi)
} else
if (corstr1=="AR-1") 
{ 
sum5<-0
sum6<-0
for ( i in  1:N)           {
for ( r in  1:nr)          {
for ( j in  1:nt[i])       {  
for ( jj in 1:nt[i])       {  
if( j>jj && abs(j-jj)==1)  {
#cat("i",i,"r",r,"j",j,"jj",jj,"\n")
sum7<-res[j+index[i]+nt[i]*(r-1)]*res[jj+index[i]+nt[i]*(r-1)]
sum5<-sum5+sum7           
#cat("i",i,"r",r,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
}
}
}
}
sum8<-nt[i]-1
sum6<-sum6+sum8
}
alfa_hat<-sum5/(nr*sum6*fi)
}  else
if (corstr1=="stat_M_dep") 
{  
alfa_hat=matrix(0,Mv,1)
for(m in 1:Mv) {
sum12<-0
sum14<-0
for ( i in  1:N)           {
for ( r in  1:nr)          {
for ( j in  1:nt[i])       {  
for ( jj in 1:nt[i])       {  
if( j>jj && abs(j-jj)==m)  {
#cat("m",m,"i",i,"r",r,"j",j,"jj",jj,"\n") 
sum11<-res[j+index[i]+nt[i]*(r-1)]*res[jj+index[i]+nt[i]*(r-1)]
sum12<-sum12+sum11 
#cat("m",m,"i",i,"r",r,"j",j,"jj",jj,"sum11",sum11,"sum12", sum12, "\n")          
} #if
}
}
}
sum13<-nt[i]-1
sum14<-sum14+sum13
}
alfa_hat[m]<-sum12/(nr*sum14*fi) 
} #m
} else
if (corstr1=="non_stat_M_dep") 
{  
alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
for( m in 1:Mv)            {
for ( j in  1:nt[1])       {  
for ( jj in 1:nt[1])       {  
if( j>jj && abs(j-jj)==m)  { 
sum16<-0                 
for ( i  in 1:N)     {
for ( r in 1:nr)     {
#cat("m",m,"j",j,"jj",jj,"i",i,"r",r,"\n") 
sum15<-res[j+index[1]+nt[i]*(r-1)]*res[jj+index[1]+nt[i]*(r-1)]
sum16<-sum15+sum16          
#cat("j",j,"jj",jj,"i",i,"r",r,"sum15",sum15,"sum16",sum16,"\n")
}
}
#cat("j",j,"jj",jj,"sum16",sum16,"\n")
alfa_hat[j,jj]<-sum16/(nr*N*fi) 
}
}
}
}
} else
if (corstr1=="unstructured") 
{  
alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
for ( j in 1:nt[1])  {  
for ( jj in 1:nt[1]) {
sum20<-0                
if (j > jj)          {
for ( i  in 1:N )    {
for ( r in 1:nr)     {
#cat("i",i,"r",r,"j",j,"jj",jj,"\n") 
sum21<-res[j+index[i]+nt[i]*(r-1)]*res[jj+index[i]+nt[i]*(r-1)]
sum20<-sum21+sum20           
}
}
#cat("j",j,"jj",jj,"sum20",sum20,"\n")
alfa_hat[j,jj]<-sum20/(nr*N*fi) 
}
}
}
} #end

if (corstr2=="independence") 
{gamma_hat<-0} else 
if (corstr2=="exchangeable") 
{
sum30<-0
for ( i in 1:N)    {
for ( j in 1:nt[i]){  
for ( r in 1:nr)   {    
for ( rr in 1:nr)  {
if  ( r != rr)     {
#cat("i",i,"j",j,"r",r,"rr",rr,"\n")
sum31<-res[j+index[i]+nt[i]*(r-1)]*res[j+index[i]+nt[i]*(rr-1)]
sum30<-sum30+sum31
#cat("i",i,"j",j,"r",r,"rr",rr,"sum31",sum31,"sum32",sum32,"\n")
}
}
}
}
}
gamma_hat<-sum30/(nr*(nr-1)*sum(nt)*fi)
} else
if (corstr2=="unstructured") 
{
gamma_hat<-matrix(0,nr,nr)
for ( r in 1:nr)   {  
for ( rr in 1:nr)  {
if  ( r > rr)      {
sum40<-0
for ( i in 1:N)       {
for ( j in 1:nt[1])   {  
sum41<-res[j+index[i]+nt[i]*(r-1)]*res[j+index[i]+nt[i]*(rr-1)]
#cat( "i",i,"j",j,"r",r,"rr",rr,"\n") 
sum40<-sum40+sum41
#cat("i",i,"j",j,"r",r,"rr",rr,"sum41",sum41,"sum40",sum40,"\n")
}
}
gamma_hat[r,rr]<-sum40/(sum(nt)*fi) 
}
}
}
} #


cor2<-matrix(0,nr,nr)
if (corstr2=="independence")                                        
{cor2<-diag(nr) } else
if (corstr2=="exchangeable")                                        
{ for (r1 in 1:nr) {
  for (r2 in 1:nr) {
   if (r1!=r2) 
   {cor2[r1,r2]<-gamma_hat} else 
   {cor2[r1,r2]<-1}
  }
  }
} else
if (corstr2=="unstructured")                                        
{ cor2<-gamma_hat+t(gamma_hat)
diag(cor2)=1
}else
if (corstr2=="fixed")
{cor2=R2
}


C1hat<-array(0,c(maxclsz,maxclsz,N))
Ehat<-array(0,c(maxclsz*nr,maxclsz*nr,N))

for(i in 1:N){

cor1<-matrix(0,nt[i],nt[i])

if (corstr1=="independence")                                        
{cor1<-diag(nt[i])} else
if (corstr1=="exchangeable")                                        
{ for (t1 in 1:nt[i]) {
  for (t2 in 1:nt[i]) {
   if (t1!=t2) 
   {cor1[t1,t2]<-alfa_hat} else 
   {cor1[t1,t2]<-1}
  }
  }
} else
if (corstr1=="AR-1")                                      
{ for (t1 in 1:nt[i]) {
  for (t2 in 1:nt[i]) {
  cor1[t1,t2]<-alfa_hat^abs(t1-t2)   
  }
  }
}  else
if (corstr1=="stat_M_dep")                                     
{ for (t1 in 1:nt[1]) {
  for (t2 in 1:nt[1]) {
     if (abs(t1-t2)==0)
     {cor1[t1,t2]<-1} else
     for(m in 1:Mv) {
     if (abs(t1-t2)==m)
     {cor1[t1,t2]<-alfa_hat[m]} 
     }
  }
  }
} else
if (corstr1=="non_stat_M_dep")                                     
{ 
 cor1=alfa_hat+t(alfa_hat)
 diag(cor1)=1
} else
if ( corstr1=="unstructured")                                        
{cor1=alfa_hat+t(alfa_hat)
 diag(cor1)=1
} else
if (corstr1=="fixed")
{cor1=R1
}

Ehat[1:(nr*nt[i]),1:(nr*nt[i]),i]<-kronecker(cor2,cor1) 
C1hat[1:nt[i],1:nt[i],i]<-cor1 
}

return(list("Ehat"=Ehat,"fi"=fi,"cor1"=C1hat,"cor2"=cor2))

}
