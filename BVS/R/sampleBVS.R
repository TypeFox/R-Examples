sampleBVS = function(data,forced=NULL,inform=FALSE,cov=NULL,rare=FALSE,mult.regions=FALSE,regions=NULL,hap=FALSE,
                     iter=10000,save.iter=0,outfile=NULL,status.file=NULL,old.results=NULL){
	library(MASS)
	library(msm)
   
    	
	##Set up sampled model matrix 
	M = dim(data)[2]-1
	c = dim(cov)[2]
	r = dim(cov)[1]
	if(length(cov)==0){
		c=1
		r=M}
	if(rare==FALSE){
		which.ind = 1:M} 
	if(rare==TRUE & mult.regions==FALSE){
	    which.ind = 1}
	if(rare==TRUE & mult.regions==TRUE){
		regions = as.factor(regions)
		which.ind = 1:length(levels(regions))}
	which = matrix(NA,nrow=iter+1,ncol=length(which.ind)+2+c)
	which.char = NULL
	if(c>0){
	  colnames(which) = c(paste("Coef.",which.ind,sep=""),paste("Alpha",1:c,sep=""),"Fitness","logPrM")}
	old.iter=0
	if(length(old.results)>0){
		old.iter = length(old.results$fitness) 
		which.char = c(which.char,old.results$which)
		which[1:old.iter,] = cbind(old.results$coef,old.results$alpha,old.results$fitness,old.results$logPrM)
	}
	
	
	
	##Initialize parameters
	a0 = qnorm((1-2^(-1/r)))
	##If old results are not present initialize randomly
	if(old.iter==0){
	  Z.current = rep(0,M)
	  Z.current[sample(1:M,5)]=1
	  a1.current = rep(0,c)}
	##If old results are present use the last model sampled
	if(old.iter>0){
	  Z.current= as.numeric(unlist(strsplit(old.results$which[old.iter],split="")))
	  if(is.matrix(old.results$alpha)==TRUE){
	    a1.current = old.results$alpha[old.iter,]}
	  if(is.matrix(old.results$alpha)==FALSE){
	  	a1.current = old.results$alpha[old.iter]}}
	  	
	##Calculate current fitness of initial parameters  
	lower.bound = rep(0,r)
	upper.bound = rep(Inf,r)
    lower.bound[Z.current==0] = -Inf
    upper.bound[Z.current==0] = 0
	t.current = rtnorm(r,lower=lower.bound,upper=upper.bound)
	fit.current = fitBVS(Z.current,data=data,forced=forced,cov=cov,a1=a1.current,rare=rare,mult.regions=mult.regions,regions=regions,
	                     hap=hap,inform=inform,which=which,which.char=which.char)
	if(old.iter==0){                     
	   which.char = cbind(which.char,paste(Z.current,collapse=""))
	   which[1,] = c(fit.current[which.ind],a1.current,fit.current[-which.ind][1:2])}
	
	
	
		
	##Run mutation for iter
	run.ind = 1:iter
	if(old.iter>0){
		run.ind = old.iter:iter}
	for(i in run.ind){
		if(length(status.file)>0){
		   cat("Processing iter",i,"\n",file=status.file,append=TRUE)}
		if(length(status.file)==0){
		   cat("Processing iter",i,"\n")}   
		   
		
		##Sample new a1 if we are using the informative prior.  If not keep a1=0
		if(inform==TRUE){
		  ##Sample a1 | current model and current t
	      V.hat = solve(diag(1,c,c)+t(cov)%*%cov)
	      alpha.hat = V.hat%*%(t(cov)%*%(t.current-a0))
	      a1.current = mvrnorm(1,mu=alpha.hat,Sigma=V.hat)
	    
	      ##Sample latent variable | current a1 and current model 
	      lower.bound = rep(0,r)
		  upper.bound = rep(Inf,r)
          lower.bound[Z.current==0] = -Inf
          upper.bound[Z.current==0] = 0
          mean = a0+as.matrix(cov)%*%as.matrix(a1.current)
	      t.current = rtnorm(r,mean=mean,lower=lower.bound,upper=upper.bound)
	      
	      ##Calculate new current fitness | current a1
	      ll = fit.current[-which.ind][1] + fit.current[-which.ind][2]
	        ##Calculate new prior
	        eta = mean
            lprob.inc = pnorm(0,mean=eta,lower.tail=FALSE,log.p=TRUE)
	        lprob.ninc = pnorm(0,mean=eta,lower.tail=TRUE,log.p=TRUE)
	        logPrM = t(Z.current)%*%lprob.inc + t(1-Z.current)%*%lprob.ninc
	      fit.current[-which.ind][1] = ll - logPrM
	      fit.current[-which.ind][2] = logPrM}
	    
	      ##Sample Model
	      model.type=0:1
		  mutate.ind = sample(1:M,1)
		  Z.new = Z.current
		  Z.new[mutate.ind] = sample(model.type[!model.type==Z.current[mutate.ind]],1)
		  
		
		  ##Calculate acceptance ratio 
		  fit.new = fitBVS(Z.new,data=data,forced=forced,cov=cov,a1=a1.current,rare=rare,mult.regions=mult.regions,regions=regions,
		                   hap=hap,inform=inform,which=which,which.char=which.char)
		  a = min(1,(exp(as.numeric(fit.new[length(which.ind)+1])-fit.current[length(which.ind)+1]))^(-1))
		  rand = runif(1)
		  ##If rand is less than or equal to a accept the new model and replace Z.current and fit.current
		  if(rand<=a){
		  	Z.current = Z.new
		  	fit.current = fit.new}	
		
		  	
		##Place fitness results of iteration in which matrix
		which.char = c(which.char,paste(Z.current,collapse=""))
		which[c(1:dim(which)[1])[is.na(which[,1])][1],] = c(fit.current[which.ind],a1.current,fit.current[-which.ind][1:2])
		
		  
		##Save results after every save.iter iteration if save.iter>0
		  
		if(i %% save.iter ==0 & save.iter>0){
			results = list(as.numeric(which[is.na(which[,1])==FALSE,(1+length(which.ind)+c)]),as.numeric(which[is.na(which[,1])==FALSE,(2+length(which.ind)+c)]),
	                       which.char,which[is.na(which[,1])==FALSE,(which.ind)],which[is.na(which[,1])==FALSE,(length(which.ind)+1):(length(which.ind)+c)])
	        names(results) = c("fitness","logPrM","which","coef","alpha")
	        save(results,file=outfile)
			}}
	
	##Return final results		
	which = which[is.na(which[,1])==FALSE,]
	results = list(as.numeric(which[,(1+length(which.ind)+c)]),as.numeric(which[,(2+length(which.ind)+c)]),which.char,
	               which[,(which.ind)],which[,(length(which.ind)+1):(length(which.ind)+c)])
	names(results) = c("fitness","logPrM","which","coef","alpha")
	return(results)
	}