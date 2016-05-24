summaryBVS = function(BVS.out,data=data,forced=NULL,cov=NULL,burnin=1000,regions=NULL,rare=FALSE,mult.regions=FALSE,inform=FALSE){
	if(burnin>0){
	  which = BVS.out$which[-c(1:burnin)]
	  if(rare==FALSE || mult.regions==TRUE){
	    coef = BVS.out$coef[-c(1:burnin),]}
	  if(rare==TRUE & mult.regions==FALSE){
        coef = BVS.out$coef[-c(1:burnin)]}
	  fitness = BVS.out$fitness[-c(1:burnin)]
	  logPrM = BVS.out$logPrM[-c(1:burnin)]
	  if(is.matrix(BVS.out$alpha)==FALSE){
	  	a1 = BVS.out$alpha[-c(1:burnin)]}
	  if(is.matrix(BVS.out$alpha)==TRUE){
	  	a1 = BVS.out$alpha[-c(1:burnin),]}
	  }
	if(burnin==0){
	  which = BVS.out$which
	  coef = BVS.out$coef
	  fitness = BVS.out$fitness
	  logPrM = BVS.out$logPrM
	  a1 = BVS.out$alpha}
	  
	  
	if(is.matrix(a1)==FALSE){
		a1 = as.matrix(a1)}
		
	
	##Parameters
	p = dim(data)[2] - 1
	snps = colnames(data)[-1]
	c = dim(cov)[2]
	r = dim(cov)[1]
	if(length(cov)==0){
	   c=1
	   r=p}
	a0 = qnorm((1-2^(-1/r)))
	
	##coef.ind
	if(rare==FALSE){coef.ind=1:p}
	if(rare==TRUE&mult.regions==FALSE){coef.ind=1}
	if(rare==TRUE&mult.regions==TRUE){
		coef.ind=1:length(unique(regions))}
	   
	##Make sure Null model is in the results
	null.fit = fitBVS(rep(0,p),data=data,forced=NULL,cov=cov,a1=rep(0,c),rare=rare,mult.regions=mult.regions,regions=regions,inform=inform,which=NULL)
	Z.null = paste(rep(0,p),collapse="")
	if(sum(which==Z.null)==0){
		which = c(Z.null,which)
		coef = rbind(rep(0,length(coef.ind)),coef)
		fitness = c(null.fit[length(coef.ind)+1],fitness)
		logPrM = c(null.fit[length(coef.ind)+2],logPrM)}
	which =as.matrix(which)	
	rownames(which) = c(1:dim(which)[1])
	
	##Post expectation of alpha
	if(is.matrix(BVS.out$alpha)==FALSE){
	  	post.alpha = mean(a1)}
	if(is.matrix(BVS.out$alpha)==TRUE){
	  	post.alpha = apply(a1,2,mean)}
		
	##If inform==TRUE average prior inclusion probability across multiple values of alpha!
	if(inform==TRUE){
	  eta = t(a0+as.matrix(cov)%*%t(a1))
	  inc.prob = function(x){mean(pnorm(0,mean=x,lower.tail=FALSE))}
	  Prior.marg = apply(eta,2,inc.prob)
	  lprob.inc = log(Prior.marg)
	  lprob.ninc = log(1-Prior.marg)
	  Nulllprob.ninc = pnorm(0,mean=a0,lower.tail=TRUE,log.p=TRUE)
	  NullPrior.marg = 1 - exp(Nulllprob.ninc)
	  ##Global Priors
	  Prior.null = exp(sum(lprob.ninc))	    
      Prior.alt = 1-Prior.null}
      
	if(inform==FALSE){
	  NullPrior.marg = (1/(p+1))
      Prior.null = .5
      Prior.alt =.5	
	}
	
	 
	##get unique models and recalc fitness with avg. prior inc probability if inform == TRUE
	u.which = unique(which)
	u.ind = as.numeric(rownames(u.which))
	u.which.char = u.which
	if(p>1){
	  u.which = t(apply(u.which,1,function(x){as.numeric(unlist(strsplit(x,split="")))}))}
	if(p==1){
	  u.which = as.matrix(apply(u.which,1,function(x){as.numeric(unlist(strsplit(x,split="")))}))}
	if(rare==FALSE || mult.regions==TRUE){
	   u.coef = coef[u.ind,]}
	if(rare==TRUE&mult.regions==FALSE){
	   u.coef = coef[u.ind]}
	u.fitness = fitness[u.ind]
	new.lPrM = logPrM[u.ind]
	if(inform==TRUE){   
	  new.lPrM = apply(u.which,1,function(x){t(x)%*%lprob.inc + t(1-x)%*%lprob.ninc})
	  u.fitness = fitness[u.ind] + logPrM[u.ind] - new.lPrM}
	
	##Calculate Posterior model probabilities
    PrMgivenD = exp(-u.fitness+min(u.fitness))/sum(exp(-u.fitness+min(u.fitness)))
	 
	
    ##Global Posterior Prob & BF
	Post.null <- PrMgivenD[u.which.char==Z.null]
	Post.alt <- 1 - Post.null
	Global.BF = (Post.alt/Post.null)
	
	
	##Marginal Inclusion probabilities
	Post.marg = t(u.which!=0)%*%as.matrix(PrMgivenD)
	Post.marg.n = t(u.which==0)%*%as.matrix(PrMgivenD)
	Global.marg.BF <- log(Post.marg) - log(Post.marg.n) - log(NullPrior.marg) +  log(1-NullPrior.marg)
	which = cbind(u.which,exp(new.lPrM),PrMgivenD)
	colnames(which) = c(colnames(data)[-1],"Prior","PostProb")
	
	##Gene Inclusion probabilities
	which.r = NULL
	Region.marg.BF = NULL
	Region.post.inc = NULL
	if(length(regions)>0){
		regions = as.factor(regions)
		u.regions = levels(regions)
		which.regions = model.matrix(~regions-1)
		which.r <- u.which %*% which.regions
        colnames(which.r) = u.regions
        probne0.r <- as.numeric((PrMgivenD * 10000) %*% (which.r > 0)) / 10000
        probe0.r <- as.numeric((PrMgivenD * 10000) %*% (which.r == 0)) / 10000
        NullPrior.marg.r <- (1 - NullPrior.marg)^apply(which.regions, 2, sum)
        Region.post.inc = probne0.r
        Region.marg.BF <- log(probne0.r) - log(probe0.r) - log(1 - NullPrior.marg.r) + log(NullPrior.marg.r)
        Region.marg.BF = exp(Region.marg.BF)
		}
	
	##Posterior estimates of coefficients and sd
	post.coef=NULL
	if(is.matrix(coef)==TRUE){
	  ##Posterior estimate of OR	
	  post.coef = t(u.coef)%*%as.matrix(PrMgivenD)
	  }
    

	##Save Results
	results=list(Global.BF,exp(Global.marg.BF),Region.marg.BF,post.alpha,post.coef,which,which.r,u.coef)
	names(results)=c("Global","MargBF","Marg.RBF","PostAlpha","PostCoef","Which","Which.r","Coef")  
    return(results)
	}