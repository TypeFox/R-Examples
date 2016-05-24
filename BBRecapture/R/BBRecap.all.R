BBRecap.all=function(data, last.column.count=FALSE, neval=1000, by.incr=1,nsim=10000,
  burnin=round(nsim/10),nsim.ML=500,burnin.ML=round(nsim.ML/10), num.t = 50, 
  prior.N = c("Rissanen","Uniform","one.over.N","one.over.N2"),       
  which.mod=c("all","standard"), sort=c("default","log.ML")){
           
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
	
	mod.M0=BBRecap(data.matrix,mod="M0",neval=neval,by.incr=by.incr,output="complete")
	mod.Mb=BBRecap(data.matrix,mod="Mb",neval=neval,by.incr=by.incr,output="complete")
	mod.Mlogistic=BBRecap(data.matrix,neval=neval,by.incr=by.incr,output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)	
	mod.Mc1=BBRecap(data.matrix,mod="Mc",markov.ord=1,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc2=BBRecap(data.matrix,mod="Mc",markov.ord=2,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc1b=BBRecap(data.matrix,mod="Mcb",markov.ord=1,neval=neval,by.incr=by.incr,output="complete")
	mod.Mc2b=BBRecap(data.matrix,mod="Mcb",markov.ord=2,neval=neval,by.incr=by.incr,output="complete")
    mod.Mt=BBRecap(data.matrix,mod="Mt",neval=neval,output="complete")
		
output=data.frame()
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt")
npar=c(2,3,3,3,5,4,6,t+1)
log.marginal.likelihood=c(mod.M0$log.marginal.likelihood,mod.Mb$log.marginal.likelihood,mod.Mlogistic$log.marginal.likelihood,mod.Mc1$log.marginal.likelihood,mod.Mc2$log.marginal.likelihood,mod.Mc1b$log.marginal.likelihood,mod.Mc2b$log.marginal.likelihood,mod.Mt$log.marginal.likelihood)
Nhat.RMSE=c(mod.M0$N.hat.RMSE,mod.Mb$N.hat.RMSE,mod.Mlogistic$N.hat.RMSE,mod.Mc1$N.hat.RMSE,mod.Mc2$N.hat.RMSE,mod.Mc1b$N.hat.RMSE,mod.Mc2b$N.hat.RMSE,mod.Mt$N.hat.RMSE)
Nhat.mean=c(mod.M0$N.hat.mean,mod.Mb$N.hat.mean,mod.Mlogistic$N.hat.mean,mod.Mc1$N.hat.mean,mod.Mc2$N.hat.mean,mod.Mc1b$N.hat.mean,mod.Mc2b$N.hat.mean,mod.Mt$N.hat.mean)
Nhat.median=c(mod.M0$N.hat.median,mod.Mb$N.hat.median,mod.Mlogistic$N.hat.median,mod.Mc1$N.hat.median,mod.Mc2$N.hat.median,mod.Mc1b$N.hat.median,mod.Mc2b$N.hat.median,mod.Mt$N.hat.median)
Nhat.mode=c(mod.M0$N.hat.mode,mod.Mb$N.hat.mode,mod.Mlogistic$N.hat.mode,mod.Mc1$N.hat.mode,mod.Mc2$N.hat.mode,mod.Mc1b$N.hat.mode,mod.Mc2b$N.hat.mode,mod.Mt$N.hat.mode)
Ninf=c(mod.M0$HPD.N[1],mod.Mb$HPD.N[1],mod.Mlogistic$HPD.N[1],mod.Mc1$HPD.N[1],mod.Mc2$HPD.N[1],mod.Mc1b$HPD.N[1],mod.Mc2b$HPD.N[1],mod.Mt$HPD.N[1])
Nsup=c(mod.M0$HPD.N[2],mod.Mb$HPD.N[2],mod.Mlogistic$HPD.N[2],mod.Mc1$HPD.N[2],mod.Mc2$HPD.N[2],mod.Mc1b$HPD.N[2],mod.Mc2b$HPD.N[2],mod.Mt$HPD.N[2])	

         output=data.frame(Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Nhat.RMSE=Nhat.RMSE,Nhat.mean=Nhat.mean,Nhat.median=Nhat.median,Nhat.mode=Nhat.mode,Ninf=Ninf,Nsup=Nsup)

if(which.mod=="all"){
	
	mod.Mlogistic.integer=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="integer",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)	
	mod.Mlogistic.counts=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML,num.t=num.t,prior.N=prior.N)	
	mod.Mlogistic.counts.integer=BBRecap(data.matrix,neval=neval,by.incr=by.incr,mbc.function="counts.integer",output="complete.ML",nsim=nsim,nsim.ML=nsim.ML,burnin=burnin,burnin.ML=burnin.ML)
	
	
	
model=c("M0","Mb","Mlogistic","Mc1","Mc2","Mc1b","Mc2b","Mt","Mlogistic.integer","Mlogistic.counts","Mlogistic.counts.integer")
npar=c(2,3,3,3,5,4,6,t+1,rep(3,3))
log.marginal.likelihood=c(mod.M0$log.marginal.likelihood,mod.Mb$log.marginal.likelihood,mod.Mlogistic$log.marginal.likelihood,mod.Mc1$log.marginal.likelihood,mod.Mc2$log.marginal.likelihood,mod.Mc1b$log.marginal.likelihood,mod.Mc2b$log.marginal.likelihood,mod.Mt$log.marginal.likelihood,mod.Mlogistic.integer$log.marginal.likelihood,mod.Mlogistic.counts$log.marginal.likelihood,mod.Mlogistic.counts.integer$log.marginal.likelihood)
Nhat.RMSE=c(mod.M0$N.hat.RMSE,mod.Mb$N.hat.RMSE,mod.Mlogistic$N.hat.RMSE,mod.Mc1$N.hat.RMSE,mod.Mc2$N.hat.RMSE,mod.Mc1b$N.hat.RMSE,mod.Mc2b$N.hat.RMSE,mod.Mt$N.hat.RMSE,mod.Mlogistic.integer$N.hat.RMSE,mod.Mlogistic.counts$N.hat.RMSE,mod.Mlogistic.counts.integer$N.hat.RMSE)
Nhat.mean=c(mod.M0$N.hat.mean,mod.Mb$N.hat.mean,mod.Mlogistic$N.hat.mean,mod.Mc1$N.hat.mean,mod.Mc2$N.hat.mean,mod.Mc1b$N.hat.mean,mod.Mc2b$N.hat.mean,mod.Mt$N.hat.mean,mod.Mlogistic.integer$N.hat.mean,mod.Mlogistic.counts$N.hat.mean,mod.Mlogistic.counts.integer$N.hat.mean)
Nhat.median=c(mod.M0$N.hat.median,mod.Mb$N.hat.median,mod.Mlogistic$N.hat.median,mod.Mc1$N.hat.median,mod.Mc2$N.hat.median,mod.Mc1b$N.hat.median,mod.Mc2b$N.hat.median,mod.Mt$N.hat.median,mod.Mlogistic.integer$N.hat.median,mod.Mlogistic.counts$N.hat.median,mod.Mlogistic.counts.integer$N.hat.median)
Nhat.mode=c(mod.M0$N.hat.mode,mod.Mb$N.hat.mode,mod.Mlogistic$N.hat.mode,mod.Mc1$N.hat.mode,mod.Mc2$N.hat.mode,mod.Mc1b$N.hat.mode,mod.Mc2b$N.hat.mode,mod.Mt$N.hat.mode,mod.Mlogistic.integer$N.hat.mode,mod.Mlogistic.counts$N.hat.mode,mod.Mlogistic.counts.integer$N.hat.mode)
Ninf=c(mod.M0$HPD.N[1],mod.Mb$HPD.N[1],mod.Mlogistic$HPD.N[1],mod.Mc1$HPD.N[1],mod.Mc2$HPD.N[1],mod.Mc1b$HPD.N[1],mod.Mc2b$HPD.N[1],mod.Mt$HPD.N[1],mod.Mlogistic.integer$HPD.N[1],mod.Mlogistic.counts$HPD.N[1],mod.Mlogistic.counts.integer$HPD.N[1])	
Nsup=c(mod.M0$HPD.N[2],mod.Mb$HPD.N[2],mod.Mlogistic$HPD.N[2],mod.Mc1$HPD.N[2],mod.Mc2$HPD.N[2],mod.Mc1b$HPD.N[2],mod.Mc2b$HPD.N[2],mod.Mt$HPD.N[2],mod.Mlogistic.integer$HPD.N[2],mod.Mlogistic.counts$HPD.N[2],mod.Mlogistic.counts.integer$HPD.N[2])	
		
output=data.frame(Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Model=model,npar=npar,log.marginal.likelihood=log.marginal.likelihood,Nhat.RMSE=Nhat.RMSE,Nhat.mean=Nhat.mean,Nhat.median=Nhat.median,Nhat.mode=Nhat.mode,Ninf=Ninf,Nsup=Nsup)

	
}


#if(ss=="npar"){
#    output=output[order(npar),]
#}
    

if(ss=="log.ML"){
    output=output[order(log.marginal.likelihood),]
}

print(paste("neval =",neval,sep=" "))
print(paste("M =",M,sep=" "))
return(output)
	
}
