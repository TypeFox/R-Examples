pssm.generate.data <-
function(theta1=.2,theta2=.2,phaz.progression=log(-log(.3)/4)*rep(1,5),
	phaz.survival=log(-log(.15)/4)*rep(1,15),accrual=3,followup=2,
  m=5,n=400,times=NULL,delta=.15,alloc=c(1,1),seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  #alloc placebo then rx
  if((accrual+followup)>m) {print('Accrual + Follow Up must less than M');return(NULL)}
  
 findint<-function(x,u)                   #returns the interval of x in the rows of u
t(apply(data.frame(inp=x,int=u),1,
		function(x){ 
			h=findInterval(x[1],x[-1]) 
		   return(c(x[-1][h],x[-1][h+1]))}))
		
cs=cumsum(exp(phaz.progression)) #cummulative hazard of progression over m intervals
fr=alloc[1]/sum(alloc)
n1=round(n*(fr))
n2=round(n*(1-fr))
n=n1+n2
p1=log(runif(n))/(-exp(c(rep(0,n1),rep(theta1,n2))))
s1=log(runif(n))/(-exp(c(rep(0,n1),rep(theta2,n2))))
cens=runif(n)*accrual+followup
if(is.null(times)){
st=matrix(0,n,5)
for (i in (1:3)) st[,i+1]=runif(n,(2*i-1)*m/8,(2*(i+1)-1)*m/8) #sample times for progression 
st[,5]=rep(m,n)}
else {
	st=matrix(0,n,length(times)+2)
	for (i in (1:length(times))) st[,i+1]=runif(n,times[i]-delta,times[i]+delta)
st[,length(times)+2]=m
		}
tprog=approx(c(0,cs),0:m,p1,rule=2)$y
tpc=pmin(tprog,cens)
cni=tprog<=cens
tpr01=findint(tpc,st)
tprog0=tpr01[,1]
tprog1=ifelse(cni,tpr01[,2],NA)
sp=((1:m)-1)*(2*(m+1)-(1:m))/2 #0  5  9 12 14
bc=matrix(0,m,m)
jprog=(1:n)[cni] #indicies of places where progression occured
tdeath=rep(m,n)
cdeath=rep(1,n)
cdeath[setdiff(1:n,jprog)]<-0

for (i in 1:m)  bc[i,i:m]=cumsum(exp(phaz.survival[(sp[i]+1):(sp[i]+1+m-i)]))
for (i in jprog){
	pr=tprog[i]
	itprog=floor(pr)+1
	d1=itprog-pr
	if (itprog<m){ 
		tm1=c(0,bc[itprog,itprog]*d1,bc[itprog,(itprog+1):m])
		tm2=c(0,(d1+0:(m-itprog)))}
	else {tm1=c(0,bc[itprog,itprog]*d1)
		tm2=c(0,d1)}
    tdth=approx(tm1,tm2,s1[i],rule=2)$y+pr
	if(tdth<cens[i]) tdeath[i]=tdth else {tdeath[i]=cens[i];cdeath[i]=0}
    if (tdeath[i]<tprog1[i]) tprog1[i]=tdeath[i] 
	
	}
return(data.frame(tprog0,tprog1,cdeath,tdeath,rx=c(rep(0,n1),rep(1,n2))))
}
