LBRecap.custom.part=function(
data,
last.column.count=FALSE,
partition,
neval=1000,
by.incr=1,
output=c("base","complete")){

 output=match.arg(output)
	
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
		
	partial=pch(data.matrix)
	
	mm1=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	mm0=matrix(NA,ncol=ncol(data.matrix),nrow=nrow(data.matrix))
	cc=c()
	n.obs.1=c()
	n.obs.0=c()
	n.unobs=c()
 	p.1=matrix(ncol = length(partition), nrow = neval)

 	log.lik=c()
	
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

	failure.conditions=F
	
	if(sum(as.logical(n.unobs))==1){failure.conditions=T}

  for(l in 1:neval){
  	
	val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)
  	nn.0=n.obs.0+(n.unobs*(l-1)*by.incr)
  	nn.1=n.obs.1

p.1[l,]=nn.1/(nn.0+nn.1)

n.00=nn.0[nn.0>0]
n.10=nn.1[nn.0>0]

n.01=nn.0[nn.1>0]
n.11=nn.1[nn.1>0]

log.lik[l]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+(sum(n.11*log(n.11/(n.11+n.01)))+(sum(n.00*log(n.00/(n.10+n.00)))))

	}

 AIC=-2*max(log.lik)+2*(length(partition)+1)

 u=which(log.lik==max(log.lik))-1
 
 N.hat=M+(u*by.incr)

 pH.hat=p.1[which(log.lik==max(log.lik)),]
 names(pH.hat)=paste("pH",1:length(pH.hat),sep="")

  ## CONFIDENCE INTERVAL (normal approximation)
 
	z2=(qnorm(0.975))^2
 
 Nplus=N.hat
 
 	if(N.hat<(M+neval*(by.incr))){
 i=1
 		
 while((2*(max(log.lik)-log.lik[(which(log.lik==(max(log.lik)))+i)])<z2 & (N.hat+i*(by.incr))< M+neval*(by.incr)))
 {
 	i=i+1 
 	}
 
 	 Nplus=N.hat+i*by.incr
 	}
 
 Nminus=N.hat
 	
 	if(N.hat>M){
 i=1	
 	
 while(2*(max(log.lik)-log.lik[(which(log.lik==(max(log.lik)))-i)])<z2 &(N.hat-i*(by.incr))>M){
 	i=i+1
 	}
 	 Nminus=N.hat-i*by.incr
	}
	

 out=switch(output,
 base=list(N.hat=N.hat,CI=c(Nminus,Nplus),AIC=AIC,L.Failure=FALSE),
 complete=list(N.hat=N.hat,CI=c(Nminus,Nplus),pH.hat=pH.hat,AIC=AIC,L.Failure=FALSE,log.lik=log.lik,N=seq(M,neval*by.incr+M,by=by.incr)[-neval+1],Partition=partition))
	
  if(failure.conditions==T){
  	n1=n.obs.1[which(n.unobs!=0)]
  	n0=n.obs.0[which(n.unobs!=0)]
  if((M*(t-1)-n0<=((M-1)*(t-1))/2-1) & n1==M) { 
  out=switch(output,
  base=list(N.hat=NA,CI=c(NA,NA),AIC=-Inf,L.Failure=TRUE),
  complete=list(N.hat=NA,CI=c(NA,NA),AIC=-Inf,L.Failure=TRUE,log.lik=log.lik,N=seq(M,neval*by.incr+M,by=by.incr)[-neval+1],Partition=partition))
		 }
	 }	
 
  return(out)
}  
