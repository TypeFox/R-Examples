BBRecap.custom.part=function(
data,
last.column.count=FALSE,
partition,
neval=1000,
by.incr=1,
prior.N=c("Rissanen","Uniform","one.over.N","one.over.N2"),
output=c("base","complete")){

 prior.N=match.arg(prior.N)
 output=match.arg(output)

 prior=switch(prior.N,
 Rissanen="Rissanen.prior",
 Uniform="Uniform.prior",
 one.over.N="1overN.prior",
 one.over.N2="1overN^2.prior"
 )	
 	
  if(!(class(data)=="data.frame" | class(data)=="matrix" | class(data)=="array" | class(data)=="table")){
       stop("input data must be a data.frame or a matrix object or an array")
   }
   
   if(class(data)=="table" | class(data)=="array"){
   	
   	n.occ=length(dim(data))
mm=matrix(ncol=n.occ,nrow=2^n.occ)

for(i in 1:2^n.occ){
mm[i,]=as.numeric(intToBits(i-1)>0)[1:n.occ]
}

data.vec=c(data)
data[1]=0
dd=c()
for(i in 1:length(data)){
		
	dd=c(dd,rep(mm[i,],data.vec[i]))
	
		}

data.matrix=matrix(dd,ncol=length(dim(data)),byrow=T)
  	
   }
   
 
   if(class(data.matrix)!="matrix"){
       data.matrix=as.matrix(data) 
       }
   
  if(last.column.count){

if (any(data[,ncol(data)]<0)){
    stop("Last column must contain non negative frequencies/counts")
}
      
  	data=as.matrix(data)
  	data.matrix=matrix(ncol=(ncol(data)-1))

for(i in 1:nrow(data)){

d=rep(data[i,1:(ncol(data)-1)],(data[i,ncol(data)]))
dd=matrix(d,ncol=(ncol(data)-1),byrow=T)
data.matrix=rbind(data.matrix,dd)

}

        
data.matrix=data.matrix[-1,]
  	  	
  }

  if(any(data.matrix!=0 & data.matrix!=1)) stop("data must be a binary matrix")
  
  if(sum(apply(data.matrix,1,sum)==0)){
 	warning("input data argument contains rows with all zeros: these rows will be removed and ignored")
  
        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

 t=ncol(data.matrix)
 M=nrow(data.matrix)

	p.p=c()

	for(r in 1:length(partition)){
	pp=unique(partition[[r]])
	p.p=c(p.p,pp)
	}

 	p.p=sort(unique(p.p))
 
 	cond1=!all(sort(unlist(list.historylabels(ncol(data.matrix))))==p.p)
 
 	cond2=max(table(unlist(partition)))>1
 
 	if(cond1 | cond2)
 		stop("The input argument 'partition' does not represent a partition of the set of all partial capture histories")

 	
 	prior.distr.N=function(x){
 	
 	if(prior=="Rissanen.prior") {out=(rissanen(x))}
 	if(prior=="Uniform.prior") {out=(1/(x^0))}
 	if(prior=="1overN.prior") {out=(1/(x^1))}
 	if(prior=="1overN^2.prior") {out=(1/(x^2))}
 	
 	return(out)
 	}
 
	partial=pch(data.matrix)
	
	mm1=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	mm0=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	cc=c()
	n.obs.1=c()
	n.obs.0=c()
	n.unobs=c()
	
	log.post.N=c()
 	post.N=c()
 	vv=c()
  	prior.inv.const=0
	
	for(k in 1:length(partition)){
	
	for(i in 1:nrow(data.matrix)){
	for(j in 1:ncol(data.matrix)){	
	mm1[i,j]=any(partition[[k]]==partial[i,j]) & data.matrix[i,j]==1
	mm0[i,j]=any(partition[[k]]==partial[i,j]) & data.matrix[i,j]==0
	}
	}
	n.obs.0[k]=sum(mm0)
	n.obs.1[k]=sum(mm1)
	}
	
	for(k in 1:length(partition)){
	
	for(j in 1:ncol(data.matrix)){	
	cc[j]=any(partition[[k]]==pch(matrix(rep(0,ncol(data.matrix)),nrow=1))[j])
	}
	n.unobs[k]=sum(cc)	
	}	

  for(l in 1:neval){
  	
	val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)
  	nn.0=n.obs.0+(n.unobs*(l-1)*by.incr)
  	nn.1=n.obs.1
 
 prior.inv.const=prior.inv.const+ prior.distr.N((M+l-1))
 
 log.post.N[l]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+sum(lbeta((nn.1+1),(nn.0+1)))+log(prior.distr.N((M+l-1)))

	}

	l.max=max(log.post.N)

	for(k in 1:neval){
	vv[k]<-exp(log.post.N[k]-l.max)
	}
	
	ss=sum(vv)
	post.N=vv/ss
	
#	return(post.N)

ord=order(post.N,decreasing=T)
p.max=ord[1]
mode.N=val[p.max]

mean.N<-round(sum(post.N*val))

median.N=M

g<-c()				
for(k in 1:neval)	{
g=c(g,post.N[k])
if (sum(g)<=0.5) median.N=median.N+1
			}

funzioneRMSE=function(x){
	sum((((x/val)-1)^2)*post.N)
	}		
	
estimate.N=round(optimize(funzioneRMSE,c(min(val),max(val)))$minimum)	

### Credible Set ###

alpha<-0.05	

g=0
d=0

aa=c()

ordine<-order(post.N,decreasing=T)

w<-val
w<-w[ordine]
p<-post.N
p<-p[ordine]

for(k in 1:neval)		{
if (g<(1-alpha)) 	{g=g+p[k]  
  		 	 d=d+1}
				}
aa<-w[1:d]
inf.lim<-min(aa)
sup.lim<-max(aa)

log.marg.likelihood=log(sum(exp(log.post.N-max(log.post.N)-log(prior.inv.const))))+max(log.post.N)

out=switch(output,
base=
list(Prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N.=c(inf.lim,sup.lim),log.marginal.likelihood=log.marg.likelihood),
complete=
list(Prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),log.marginal.likelihood=log.marg.likelihood,
N.range=val,posterior.N=post.N,Partition=partition)) 

  return(out)

}
