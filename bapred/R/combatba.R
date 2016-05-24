combatba <-
function(x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

    mod=NULL; numCovs = NULL; par.prior=TRUE; prior.plots=FALSE

    if(length(unique(batch))==1) {
      params <- list(xadj=x, meanoverall=colMeans(x), var.pooled=matrix(nrow=ncol(x), ncol=1, data=apply(x, 2, var)*(1 - 1/nrow(x))))
      params$batch <- batch
	  params$nbatches <- length(unique(batch))
	
      class(params) <- "combat"
	
      return(params)    
	}
	  
	  
	mod = cbind(mod,batch)
	
	# check for intercept, and drop if present
	check = apply(mod, 2, function(x) all(x == 1))
	mod = as.matrix(mod[,!check])
	
	colnames(mod)[ncol(mod)] = "Batch"
	
	if(sum(check) > 0 & !is.null(numCovs)) numCovs = numCovs-1
	
	design <- design.mat(mod,numCov = numCovs)	

	batches <- list.batch(mod)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(x))
	if(NAs){stop(paste("Found", sum(is.na(x)), "Missing data Values"))}
        #print(x[1:2,])
	##standardize data across genes
	# cat('Standardizing data across genes\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%as.matrix(x)}
    ## if (NAs){B.hat=apply(t(x),1,Beta.NA,design)} #Standarization Model
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((t(x)-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}
      if (NAs){var.pooled <- apply(t(x)-t(design%*%B.hat),1,var,na.rm=T)}

	meanoverall <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;meanoverall <- meanoverall+t(tmp%*%B.hat)}	
	s.data <- (t(x)-meanoverall)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	# cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){
		gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
	} 
    ##  if (NAs){
	##	gamma.hat=apply(s.data,1,Beta.NA,batch.design)
	##}
	
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
		}

	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		# cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],
				delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}
	##else{
	##	# cat("Finding nonparametric adjustments\n")
	##	for (i in 1:n.batch){
	##		temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
	##		gamma.star <- rbind(gamma.star,temp[1,])
	##		delta.star <- rbind(delta.star,temp[2,])
	##		}
	##	}


	### Normalize the data ###
	# cat("Adjusting the data\n")

	bayesdata <- s.data
	j <- 1
	for (i in 1:length(batches)){
            id = batches[[i]]
		bayesdata[,id] <- (bayesdata[,id]-t(batch.design[id,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+meanoverall
	
    params <- list(xadj=t(bayesdata), meanoverall=meanoverall[,1], var.pooled=var.pooled)
	params$batch <- batch
	params$nbatches <- length(unique(batch))
	
    class(params) <- "combat"
	
    return(params)

}
