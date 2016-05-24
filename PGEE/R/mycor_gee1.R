mycor_gee1 <-
function(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R,scale.fix,scale.value) {

eta=X%*%beta_new
mu=family$linkinv(eta)
sd=sqrt(family$variance(mu))
res<-(as.vector(y)-mu)/sd 

if(scale.fix==0) {
fi<-sum(res^2)/(sum(nt))
} else
if(scale.fix==1) {
fi<-scale.value
}

aindex=cumsum(nt)
index=c(0,aindex[-length(aindex)])

if (corstr=="independence") 
{alfa_hat<-0} else 
if (corstr=="exchangeable") 
{
sum1<-0
sum3<-0
for ( i in  1:N)          {
for ( j in  1:nt[i])      {    
for ( jj in 1:nt[i])      {    
if  ( j!=jj)              {
#cat("i",i,"j",j,"jj",jj,"\n")
sum2<-res[j+index[i]]*res[jj+index[i]]
#cat("i",i,"j",j,"jj",jj,"\n")
sum1<-sum1+sum2
#cat("i",i,"j",j,"jj",jj,"sum2",sum2,"sum1",sum1,"\n")
}
}
}
sum4<-nt[i]*(nt[i]-1)
sum3<-sum3+sum4
} #i
alfa_hat<-sum1/(sum3*fi)
} else
if (corstr=="AR-1") 
{ 
sum5<-0
sum6<-0
for ( i in  1:N)           {
for ( j in  1:nt[i])       {  
for ( jj in 1:nt[i])       {  
if( j>jj && abs(j-jj)==1)  {
#cat("i",i,"j",j,"jj",jj,"\n")
sum7<-res[j+index[i]]*res[jj+index[i]]
sum5<-sum5+sum7           
#cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
}
}
}
sum8<-(nt[i]-1)
sum6<-sum6+sum8
} #i
alfa_hat<-sum5/(sum6*fi)
} else
if (corstr=="stat_M_dep") 
{  
alfa_hat=matrix(0,Mv,1)
for(m in 1:Mv) {
sum12<-0
sum14<-0
for ( i in  1:N)           {
for ( j in  1:nt[i])       {  
for ( jj in 1:nt[i])       {  
if( j>jj && abs(j-jj)==m)  {
#cat("m",m,"i",i,","j",j,"jj",jj,"\n") 
sum11<-res[j+index[i]]*res[jj+index[i]]
sum12<-sum12+sum11 
#cat("m",m,"i",i,"j",j,"jj",jj,"sum11",sum11,"sum12", sum12, "\n")          
} #if
}
}
sum13<-nt[i]-1
sum14<-sum14+sum13
} #i
alfa_hat[m]<-sum12/(sum14*fi) 
} #m
}  else
if (corstr=="non_stat_M_dep") 
{  
alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
for( m in 1:Mv)            {
for ( j in  1:nt[1])       {  
for ( jj in 1:nt[1])       {  
if( j>jj && abs(j-jj)==m)  { 
sum16<-0                 
for ( i  in 1:N)     {
#cat("m",m,"j",j,"jj",jj,"i",i"\n") 
sum15<-res[j+index[i]]*res[jj+index[i]]
sum16<-sum15+sum16          
#cat("j",j,"jj",jj,"i",i,"sum15",sum15,"sum16",sum16,"\n")
} #i
#cat("j",j,"jj",jj,"sum16",sum16,"\n")
alfa_hat[j,jj]<-sum16/(N*fi) 
}
}
}
}
} else
if (corstr=="unstructured") 
{  
alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
for ( j in 1:nt[1])  {  
for ( jj in 1:nt[1]) {
sum20<-0                
if (j > jj)          {
for ( i  in 1:N )    {
#cat("i",i,"j",j,"jj",jj,"\n") 
sum21<-res[j+index[i]]*res[jj+index[i]]
sum20<-sum21+sum20           
} #i
#cat("j",j,"jj",jj,"sum20",sum20,"\n")
alfa_hat[j,jj]<-sum20/(N*fi) 
}
}
}
} else
if (corstr=="fixed")
{alfa_hat=NULL
}

Ehat<-array(0,c(maxclsz,maxclsz,N))

for(i in 1:N){
cor1<-matrix(0,nt[i],nt[i])
if (corstr=="independence")                                        
{cor1<-diag(nt[i])} else
if (corstr=="exchangeable")                                        
{ for (t1 in 1:nt[i]) {
  for (t2 in 1:nt[i]) {
   if (t1!=t2) 
   {cor1[t1,t2]<-alfa_hat} else 
   {cor1[t1,t2]<-1}
  }
  }
} else
if (corstr=="AR-1")                                      
{ for (t1 in 1:nt[i]) {
  for (t2 in 1:nt[i]) {
  cor1[t1,t2]<-alfa_hat^abs(t1-t2)   
  }
  }
}  else
if (corstr=="stat_M_dep")                                     
{ for (t1 in 1:nt[i]) {
  for (t2 in 1:nt[i]) {
     if (abs(t1-t2)==0)
     {cor1[t1,t2]<-1} else
     for(m in 1:Mv) {
     if (abs(t1-t2)==m)
     {cor1[t1,t2]<-alfa_hat[m]} 
     }
  }
  }
} else
if (corstr=="non_stat_M_dep")                                     
{ 
 cor1=alfa_hat+t(alfa_hat)
 diag(cor1)=1
} else
if (corstr=="unstructured")                                        
{cor1=alfa_hat+t(alfa_hat)
 diag(cor1)=1
} else
if (corstr=="fixed")
{cor1=R
}

Ehat[1:nt[i],1:nt[i],i]<-cor1 
}

return(list("Ehat"=Ehat,"fi"=fi))

}
