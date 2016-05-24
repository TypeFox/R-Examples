combatbaaddon <-
function(params, x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

   if(class(params) != "combat")
     stop("Input parameter 'params' has to be of class 'combat'.")
	 
   if(ncol(params$xadj) != ncol(x))
     stop("Number of variables in test data matrix different to that of training data matrix.")	 
	 
	mod = matrix(nrow=length(batch), ncol=1, data=batch)
	colnames(mod)[ncol(mod)] = "Batch"

	design <- design.mat(mod,numCov = NULL)	

	batches <- list.batch(mod)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	##standardize Data across genes

      B.hat <- solve(t(design)%*%design)%*%t(design)%*%as.matrix(x)

	meanoverall <- matrix(nrow=ncol(x), ncol=n.array, data=params$meanoverall)
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;meanoverall <- meanoverall+t(tmp%*%B.hat)}	
	s.data <- (t(x)-meanoverall)/(sqrt(params$var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	# cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch, drop=FALSE]		
		

	gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))

	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
	}


	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)


	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL

	# cat("Finding parametric adjustments\n")
	for (i in 1:n.batch){
		temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],
			delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
		gamma.star <- rbind(gamma.star,temp[1,])
		delta.star <- rbind(delta.star,temp[2,])
	}


	### Normalize the Data ###
	# cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in 1:length(batches)){
            id = batches[[i]]
		bayesdata[,id] <- (bayesdata[,id]-t(batch.design[id,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(params$var.pooled)%*%t(rep(1,n.array))))+meanoverall

      xadj <- t(bayesdata)

      return(xadj)

}
