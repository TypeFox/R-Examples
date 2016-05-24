LBRecap.all=function(data, last.column.count=FALSE, neval=1000, by.incr=1,which.mod=c("all","standard"),sort=c("default","AIC")){

    ss=sort[1]
    mm=which.mod[1] 
 
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
   
  if(last.column.count==TRUE){
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
        #### BETTER USING data.frame and subset(...)
        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }
  	
 t=ncol(data.matrix)
 M=nrow(data.matrix)		
	
	mod.M0=LBRecap(data.matrix,mod="M0",neval=neval,by.incr=by.incr)
	mod.Mb=LBRecap(data.matrix,mod="Mb",neval=neval,by.incr=by.incr)
	mod.Mlogistic=LBRecap(data.matrix,neval=neval,by.incr=by.incr)	
	mod.Mc1=LBRecap(data.matrix,mod="Mc",markov.ord=1,neval=neval,by.incr=by.incr)
	mod.Mc2=LBRecap(data.matrix,mod="Mc",markov.ord=2,neval=neval,by.incr=by.incr)
	mod.Mc1b=LBRecap(data.matrix,mod="Mcb",markov.ord=1,neval=neval,by.incr=by.incr)
	mod.Mc2b=LBRecap(data.matrix,mod="Mcb",markov.ord=2,neval=neval,by.incr=by.incr)
        mod.Mt=LBRecap(data.matrix,mod="Mt",neval=neval)
		
output=data.frame()
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt")
npar=c(2,3,3,3,5,4,6,t+1)
AIC=c(mod.M0$AIC,mod.Mb$AIC,mod.Mlogistic$AIC,mod.Mc1$AIC,mod.Mc2$AIC,mod.Mc1b$AIC,mod.Mc2b$AIC,mod.Mt$AIC)
Nhat=c(mod.M0$N.hat,mod.Mb$N.hat,mod.Mlogistic$N.hat,mod.Mc1$N.hat,mod.Mc2$N.hat,mod.Mc1b$N.hat,mod.Mc2b$N.hat,mod.Mt$N.hat)
Ninf=c(mod.M0$CI[1],mod.Mb$CI[1],mod.Mlogistic$CI[1],mod.Mc1$CI[1],mod.Mc2$CI[1],mod.Mc1b$CI[1],mod.Mc2b$CI[1],mod.Mt$CI[1])
Nsup=c(mod.M0$CI[2],mod.Mb$CI[2],mod.Mlogistic$CI[2],mod.Mc1$CI[2],mod.Mc2$CI[2],mod.Mc1b$CI[2],mod.Mc2b$CI[2],mod.Mt$CI[2])	

         output=data.frame(Model=model,npar=npar,AIC=AIC,Nhat=Nhat,Ninf=Ninf,Nsup=Nsup)

if(which.mod=="all"){
	
	mod.Mlogistic.integer=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="integer")	
	mod.Mlogistic.counts=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts")	
	mod.Mlogistic.counts.integer=LBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts.integer")
	
	
	
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt","Mlogistic.integer","Mlogistic.counts","Mlogistic.counts.integer")
npar=c(2,3,3,3,5,4,6,t+1,rep(3,3))
AIC=c(mod.M0$AIC,mod.Mb$AIC,mod.Mlogistic$AIC,mod.Mc1$AIC,mod.Mc2$AIC,mod.Mc1b$AIC,mod.Mc2b$AIC,mod.Mt$AIC,mod.Mlogistic.integer$AIC,mod.Mlogistic.counts$AIC,mod.Mlogistic.counts.integer$AIC)
Nhat=c(mod.M0$N.hat,mod.Mb$N.hat,mod.Mlogistic$N.hat,mod.Mc1$N.hat,mod.Mc2$N.hat,mod.Mc1b$N.hat,mod.Mc2b$N.hat,mod.Mt$N.hat,mod.Mlogistic.integer$N.hat,mod.Mlogistic.counts$N.hat,mod.Mlogistic.counts.integer$N.hat)
Ninf=c(mod.M0$CI[1],mod.Mb$CI[1],mod.Mlogistic$CI[1],mod.Mc1$CI[1],mod.Mc2$CI[1],mod.Mc1b$CI[1],mod.Mc2b$CI[1],mod.Mt$CI[1],mod.Mlogistic.integer$CI[1],mod.Mlogistic.counts$CI[1],mod.Mlogistic.counts.integer$CI[1])	
Nsup=c(mod.M0$CI[2],mod.Mb$CI[2],mod.Mlogistic$CI[2],mod.Mc1$CI[2],mod.Mc2$CI[2],mod.Mc1b$CI[2],mod.Mc2b$CI[2],mod.Mt$CI[2],mod.Mlogistic.integer$CI[2],mod.Mlogistic.counts$CI[2],mod.Mlogistic.counts.integer$CI[2])	
		
output=data.frame(Model=model,npar=npar,AIC=AIC,Nhat=Nhat,Ninf=Ninf,Nsup=Nsup)

	
}


#if(ss=="npar"){
#    output=output[order(npar),]
#}
    

if(ss=="AIC"){
    output=output[order(AIC),]
}

print(paste("neval =",neval,sep=" "))
print(paste("M =",M,sep=" "))
return(output)
	
}
