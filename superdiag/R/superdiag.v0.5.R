
## CREATE THE DIAGNOSTIC FUNCTION 10312011, TSUNG-HAN TSAI AND JEFF GILL
superdiag <- function(mcmcoutput, burnin=10000, confidence.gr=0.95, frac1.gw=0.1, frac2.gw=0.5, eps.hw=0.1, 
                       pvalue.hw=0.05, q.rl=0.025, r.rl=0.005, s.rl=0.95, eps.rl=0.001) {
	
	# mcmcoutput: input chains from jags, bugs, etc.
	# burnin: the number of burn-in iterations for the sampler
	# confidence.gr: 1-alpha for the Gelman and Rubin test
	# frac1.gw: frac1 for the Geweke test
	# frac2.gw: frac2 for the Geweke test
	# eps.hw: epsilon for the Heidelberger and Welch test
	# pvalue.hw: p-value for the Heidelberger and Welch test
	# q.rl: q-parameter for the Raftery and Lewis test
	# r.rl: r-parameter for the Raftery and Lewis test
	# s.rl: s-parameter for the Raftery and Lewis test
	# eps.rl: convergence epsilon for the Raftery and Lewis test
	
	# CREATE A FUNCTION TO DISCARD THE BURN-IN PERIOD
	burn <- function(input.matrix, burnin) {
		out <- input.matrix[-(1:burnin),];
		return(out);
	}
	
	# THE INPUTS SHOULD BE SAVED AS "mcmc" OR "mcmc.list" CLASS
	if (class(mcmcoutput) != "mcmc" & class(mcmcoutput) != "mcmc.list" & class(mcmcoutput) != "list") 
		stop("The inputs have to be mcmc, mcmc.list, or list objects.");
	# CONVERT "mcmc" INTO "mcmc.list"
	if (class(mcmcoutput) == "mcmc") {
		mcmcoutput <- as.mcmc.list(mcmcoutput);
	}
	
	para.names <- dimnames(mcmcoutput[[1]])[[2]];  # PARAMETERS NAMES
	n.chains <- length(mcmcoutput);  # THE NUMBER OF CHAINS
	dim.chain <- sapply(mcmcoutput, dim);  # THE DIMENSION OF EACH CHAIN
	t.iter <- dim.chain[1];	  # THE NUMBER OF ITERATIONS BEFORE DELETING THE BURN-IN PERIOD
	diff.dim <- dim.chain - dim.chain;
	if (sum(diff.dim != 0) != 0) stop("The number of iterations or variables is not equal for all chains.");
	
	# DISCARD THE BURN-IN PERIOD
	if (burnin <= 0) {
		warning("The burn-in period is negative or zero");
		mcmcburnin <- mcmcoutput;
	}
	else {
		if (burnin >= t.iter) {
			stop("The number of iterations is less than the burn-in period.");
		}
		else {
			mcmcburnin <- lapply(mcmcoutput, burn, burnin=burnin);
		}
	}
	
	# SAVE THE SAMPLES AS A MCMC LIST AFTER DISCARDING THE BURN-IN PERIOD
	mcmcburnin.list <- vector("list", n.chains);
	for (i in 1:n.chains) {
		mcmcburnin.list[[i]] <- as.mcmc(mcmcburnin[[i]]);		
	}	
	mcmcburnin.mcmclist <- as.mcmc.list(mcmcburnin.list);	
	
	# THE TOTAL SAMPLES FOR ALL CHAINS
	t.samples <- as.matrix(mcmcburnin.mcmclist);

	# REARRANGE THE PRINTED RESULTS	
	geweke.chains <- matrix(NA, nrow=n.chains, ncol=dim.chain[2]);
	heidel.list <- vector("list", n.chains);
	raftery.list <- vector("list", n.chains);
	
	# SETUP DIFFERENT WINDOW SPECIFICATIONS FOR GEWEKE
	geweke.windows <- matrix(c(frac1.gw,frac2.gw),ncol=2)
	for (i in 2:n.chains) {
		win1 <- runif(1,0,0.99); win2 <- 1-runif(1,win1,1)
		geweke.windows <- rbind(geweke.windows,c(win1,win2))
	}

	# SETUP DIFFERENT PARAMETER SPECIFICATIONS FOR HEIDELBERGER AND WELCH
	heidel.params <- matrix(c(eps.hw,pvalue.hw),ncol=2)
	pvals <- c(0.1,0.05,0.025,0.01,0.005)
	for (i in 2:n.chains) {
		param1 <- runif(1,0.01,0.2); param2 <- sample(x=pvals,size=1)
		heidel.params <- rbind(heidel.params,c(param1,param2))
	}

	# SETUP DIFFERENT PARAMETER SPECIFICATIONS FOR RAFTERY AND LEWIS
	raft.params <- matrix(c(q.rl, r.rl, s.rl, eps.rl),ncol=4)
	qvals <-c(0.25,0.1,0.05,0.01,0.001)
	rvals <- c(0.001,0.0025,0.0005,0.001,0.005)
	svals <- c(0.9,0.95,0.975,0.99,0.999)
	evals <- c(0.005,0.0025,0.001,0.0005,0.0002)
	for (i in 2:n.chains) {
		param1 <- sample(x=qvals,size=1); param2 <- sample(x=rvals,size=1)
		param3 <- sample(x=svals,size=1); param4 <- sample(x=evals,size=1)
		raft.params <- rbind(raft.params,c(param1,param2,param3,param4))
	}

	# RUN DIAGNOSTICS BY CHAIN IN ORDER TO AVOID ONE CHAIN RUINS ALL CHAINS
	for (i in 1:n.chains) {
		geweke <- suppressWarnings(try(geweke.diag(mcmcburnin.mcmclist[[i]], geweke.windows[i,1], geweke.windows[i,2]), silent=TRUE));
		if (class(geweke) == "geweke.diag")  geweke.chains[i,] <- t(geweke[1]$z);
		heidel.list[[i]] <- suppressWarnings(try(heidel.diag(mcmcburnin.mcmclist[[i]], heidel.params[i,1], heidel.params[i,2]), silent=TRUE));
		raftery.list[[i]] <- suppressWarnings(try(raftery.diag(mcmcburnin.mcmclist[[i]], q=raft.params[i,1],r=raft.params[i,2],
					s=raft.params[i,3],converge.eps=raft.params[i,4]), silent=TRUE));		
	}	
	colnames(geweke.chains) <- para.names;
	
	# PROVIDE THE BASIC INFORMATION OF MCMC SAMPLES
	cat(paste("Number of chains =", n.chains, "\n"));
	cat(paste("Number of iterations =", t.iter, "per chain before discarding the burn-in period\n"));
	cat(paste("The burn-in period =", burnin, "per chain\n"))
	cat(paste("Sample size in total =", dim(t.samples)[1], "\n"))			
	cat("\n");	

	# REPORT RESULTS OF DIAGNOSTICS OF CONVERGENCE
	if (n.chains < 2) {
		chain.name <- "chain1";		
		rownames(geweke.chains) <- chain.name;		
		cat("The Gelman-Rubin diagnostic is not reported since the number of chains is less than 2.\n");
		cat("\n");
		cat("The Geweke diagnostic:\n");
		cat(paste("Fraction in 1st window =", frac1.gw, "\n"));
		cat(paste("Fraction in 2nd window =", frac2.gw, "\n"));
		cat(paste("Z-scores:\n"));
		print(t(geweke.chains));
		cat("\n");
		cat("The Heidelberger-Welch diagnostic:\n");
		cat("\n");
		print(heidel.list);
		cat("\n");
		cat("The Raftery-Lewis diagnostic:\n");
		cat("\n");
		print(raftery.list);
		cat("\n");
	}
	else {
		chain.name <- "chain1";
		for (i in 2:n.chains) {
			namei <- paste("chain ", i, sep="");
			chain.name <- c(chain.name, namei);
		}
		rownames(geweke.chains) <- chain.name;

		cat("********** The Geweke diagnostic: **********\n");
		dimnames(geweke.windows)[[2]] <- c("Window From Start","Window From Stop")
		cat(paste("Z-scores:\n"));
		print(rbind( t(geweke.chains), t(round(geweke.windows,5))));
		cat("\n");
		cat("********** The Gelman-Rubin diagnostic: **********\n");
		print(gelman.diag(mcmcburnin.mcmclist, confidence.gr));
		cat("\n");
		cat("********** The Heidelberger-Welch diagnostic: **********\n");
		cat("\n");
		for (i in 1:n.chains)  {
			cat(paste("Chain ",i,", epsilon=",round(heidel.params[i,1],3),", alpha=",round(heidel.params[i,2],3),sep=""))
			print(heidel.list[[i]]);
			cat("\n");
		}
		cat("********** The Raftery-Lewis diagnostic: **********\n");
		cat("\n");
		for (i in 1:n.chains)  {
			cat(paste("Chain ",i,", converge.eps = ",round(raft.params[i,4],4),sep=""))
			print(raftery.list[[i]]);
			cat("\n");
		}
		cat("\n");	
	}
	# RETURN THE MCMC SAMPLES WITH THE BURN-IN DISCARDED
	superdiag.chains <- list(mcmc.samples = mcmcburnin.mcmclist);
}
