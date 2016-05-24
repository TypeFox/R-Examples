posterior.predictive.sample <-
function(MCMC.output,posterior.predictive.sample.size,output.file,prefix=''){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
			mcmc.generations <- length(which(a0 != 0))
			n.pops <- nrow(last.params$counts)
			n.loci <- ncol(last.params$counts)
				sampled.generations <- sample(c(1:mcmc.generations),posterior.predictive.sample.size,replace=TRUE)
			posterior.sample.a0 <- a0[sampled.generations]
			posterior.sample.aD <- aD[sampled.generations]
			posterior.sample.aE <- matrix(aE[,sampled.generations],nrow=length(last.params$E),ncol=posterior.predictive.sample.size)
			posterior.sample.a2 <- a2[sampled.generations]
			posterior.sample.beta <- beta[sampled.generations]
				if(exists("phi_mat")){
					posterior.sample.phi_mat <- phi_mat[,sampled.generations]
				}
			posterior.sample.covariance <- matrix(0,nrow=n.pops,ncol=n.pops)
			posterior.sample.thetas <- matrix(0,nrow=n.pops,ncol=n.loci) 
			posterior.sample.mu <- matrix(0,nrow=n.pops,ncol=n.loci)  
			posterior.sample.allele.frequencies <- matrix(0,nrow=n.pops,ncol=n.loci)  
			posterior.sample.allele.counts <- matrix(0,nrow=n.pops,ncol=n.loci)  
			posterior.sample.Fst <- array(dim=c(n.pops,n.pops,posterior.predictive.sample.size))
			
			progress <- txtProgressBar(min=0,posterior.predictive.sample.size,char="|",style=3)
				for(i in 1:posterior.predictive.sample.size){
					posterior.sample.covariance <- Covariance(posterior.sample.a0[i],
															posterior.sample.aD[i],
															posterior.sample.aE[,i],
															posterior.sample.a2[i],
															last.params$D,
															last.params$E,
															last.params$delta)
					posterior.sample.thetas <- t(mvrnorm(last.params$loci,numeric(n.pops),posterior.sample.covariance))
					posterior.sample.mu <- matrix(rnorm(n.loci,0,sd=sqrt(1/posterior.sample.beta[i])),nrow=n.pops,ncol=n.loci,byrow=TRUE)
					posterior.allele.frequencies <- transform_frequencies(posterior.sample.thetas,posterior.sample.mu)
						if(!exists("phi_mat")){
							posterior.sample.allele.counts <- simulate_allele_count_data(posterior.allele.frequencies,last.params$sample_sizes)
						}
						if(exists("phi_mat")){
							posterior.sample.allele.counts <- simulate_allele_count_data(posterior.allele.frequencies,last.params$sample_sizes,posterior.sample.phi_mat[,i])							
						}
					posterior.sample.Fst[,,i] <- calculate.all.pairwise.Fst(posterior.sample.allele.counts,last.params$sample_sizes)
						setTxtProgressBar(progress,i)
				}
			
			observed.Fst <- calculate.all.pairwise.Fst(last.params$counts,last.params$sample_sizes)
				D <- last.params$D
				E <- last.params$E
			save(observed.Fst,posterior.sample.Fst,D,E,posterior.predictive.sample.size,file=paste(prefix,paste(output.file,".Robj",sep=''),sep=''))
		})
	}
