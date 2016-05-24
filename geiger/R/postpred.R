## Exposed functions


pp.mcmc <- function(phy, d, Ngens = 1000000, sampleFreq = 1000, printFreq = 1000, prop.width = 1, model = "BM", eb.type = "exponential", clade = NULL, rlm.maxit = 20){
	
	## in the functional call, most things should be self explanatory. the only one that isn't is eb.type. this just specifies whether the node.height test is testing for a linear or exponential decline in rates. options are "exponential" or "linear".
	
	## initial housekeeping and data extraction ##
	td <- treedata(phy, d);
	phy <- td$phy;
	d <- td$data[,1]; names(d) <- rownames(td$data);
	ltt <- sort(branching.times(phy), decreasing = TRUE)
    ltt <- c(0, (max(ltt) - ltt)/max(ltt))
	bt <- branching.times(phy);
	bt <- abs(bt  - max(bt));
	disparity.data <- dtt(phy, d, plot=FALSE); 
	
	if(model == "BM") {
		output <- matrix(data = NA, nrow =(Ngens/sampleFreq)+1, ncol = 6) ;
		colnames(output) <- c("generation", "logLk", 'Sigma',"node.height.slope.lm", "node.height.slope.rlm", "MDI" );
		ephy <- phy; # tree used for likelihood computation later
		start.sig <- exp(runif(1,-20, 5))
		start.lk <- .reml.lk(phy, d, start.sig);
		sim.d <- sim.char(phy, par = start.sig)[,,1];
		
	}
			
	if(model == "EB") {
		acdc.prior <- .pp.acdc.prior(phy, increase.max = 0);

		start.a <- runif(1, acdc.prior$min, 0);
		start.sig <- exp(runif(1, 0, 5))
		output <- matrix(data = NA, nrow =(Ngens/sampleFreq)+1, ncol = 7) ;
		colnames(output) <- c("generation", "logLk", 'Sigma', "a", "node.height.slope.lm", "node.height.slope.rlm", "MDI");
		ephy <- rescale.phylo(phy, model="EB", start.a);
		start.lk <- .reml.lk(ephy, d, start.sig)
		sim.d <- sim.char(ephy, par = start.sig)[,,1];
			
	} 

	if(model == "edge.shift") {
		if(is.null(clade)) stop("you must specify a clade in which the putative shift occured to use this option");
		start.edge.scalar <- runif(1, 0, 5);
		start.sig <- exp(runif(1, 0, 5));
		output <- matrix(data = NA, nrow =(Ngens/sampleFreq)+1, ncol = 7) ;
		colnames(output) <- c("generation", "logLk", 'Sigma', "edge.scalar", "node.height.slope.lm", "node.height.slope.rlm", "MDI");
		ephy <- .scale.edge.tree(phy, clade, start.edge.scalar);
		start.lk <- .reml.lk(ephy, d, start.sig);
		sim.d <- sim.char(ephy, par = start.sig)[,,1];
			
	} 

		
	## compute the dtt and node height statistics	
	disparity.sim <- dtt(phy, sim.d, plot=FALSE);	
	a.b.c <- as.numeric(.area.between.curves(ltt, disparity.sim$dtt, disparity.data$dtt))
	pics <- abs(pic(sim.d, phy));
	if(eb.type == "exponential") {
		pics <- log(pics);
	}	
	post.lm.slope <- summary(lm(pics~bt))$coefficients[2,1];
	while(1) {
		post.rlm.slope <- suppressWarnings(rlm(pics~bt, maxit = rlm.maxit));
		if(post.rlm.slope$converged  == F) {
			rlm.maxit <- rlm.maxit + 10;
			post.rlm.slope <- suppressWarnings(rlm(pics~bt, maxit = rlm.maxit));
		} else if(post.rlm.slope$converged  == T)
		break 
	}
	post.rlm.slope <- as.numeric(post.rlm.slope$coefficients[2]);
	
	
	## do ML simulation ##
	
	if(model == "BM") {
		output[1, ] <- c(0, start.lk, start.sig, post.lm.slope, post.rlm.slope, a.b.c);
		curr <- c(start.lk, start.sig);
	} else if(model == "EB") {
		output[1, ] <- c(0, start.lk, start.sig, start.a, post.lm.slope, post.rlm.slope, a.b.c);
		curr <- c(start.lk, start.sig, start.a);		
	} else if(model == "edge.shift") {
		output[1, ] <- c(0, start.lk, start.sig, start.edge.scalar, post.lm.slope, post.rlm.slope, a.b.c);
		curr <- c(start.lk, start.sig, start.edge.scalar);		
	}

	
	gen.counter <- 1;
	acc.counter <- 0;
	print(c("Gen", "Accept" ));
	
	for(ngen in 1:Ngens) {
		
		prop.sig <- exp(.pp.getproposal(log(curr[2]), prop.width));
		if(model == "EB") {
			prop.a <- .pp.getproposal(curr[3], prop.width/10, bound = T, max = acdc.prior[[2]], min = acdc.prior[[1]]);
			ephy <- rescale.phylo(phy, model="EB", prop.a);
		}
		
		if(model == "edge.shift") {
			prop.edge.scalar <- .pp.getproposal(curr[3], prop.width/10, bound = T, min = 0, max = 50);
			ephy <- .scale.edge.tree(phy, clade, prop.edge.scalar);		
			}
		
		prop.lk <- .reml.lk(ephy, d, prop.sig);
		
		p <- runif(1)
		
		if(exp(prop.lk - curr[1]) >=p) {
			if(model == "BM") {
				curr <- c(prop.lk, prop.sig);
			} else if(model == "EB") {
				curr <- c(prop.lk, prop.sig, prop.a);
			} else if(model == "edge.shift") {
				curr <- c(prop.lk, prop.sig, prop.edge.scalar);
			}
			acc.counter <- acc.counter + 1;	
		} 		
		if(ngen %% sampleFreq == 0) {
			gen.counter <- gen.counter + 1;
			sim.d <- sim.char(ephy, par = curr[2])[,,1];
			disparity.sim <- dtt(phy, sim.d, plot=FALSE);
			a.b.c <- as.numeric(.area.between.curves(ltt, disparity.sim$dtt, disparity.data$dtt))
			pics <- abs(pic(sim.d, phy));
			if(eb.type == "exponential") {
			pics <- log(pics);
			}	
			post.lm.slope <- summary(lm(pics~bt))$coefficients[2,1];
			while(1) {
				post.rlm.slope <- suppressWarnings(rlm(pics~bt, maxit = rlm.maxit));
				if(post.rlm.slope$converged  == F) {
					rlm.maxit <- rlm.maxit + 10;
					post.rlm.slope <- suppressWarnings(rlm(pics~bt, maxit = rlm.maxit));
				} else if(post.rlm.slope$converged  == T)
				break 
			}
			post.rlm.slope <- as.numeric(post.rlm.slope$coefficients[2]);
					
			if(model == "BM") {
				output[gen.counter, ] <- c(ngen, curr, post.lm.slope, post.rlm.slope, a.b.c);
			} else if(model == "EB") {
				output[gen.counter, ] <- c(ngen, curr, post.lm.slope, post.rlm.slope, a.b.c);
			} else if(model == "edge.shift") {
				output[gen.counter, ] <- c(ngen, curr, post.lm.slope, post.rlm.slope, a.b.c);
			}

		}
		
		if(ngen %% printFreq == 0) {
			print(c(round(ngen), round(acc.counter/ngen, 2) ));
		}
	}
	return(as.data.frame(output));
}









nh.test <- function(phy, d, regression.type, log = TRUE, rlm.maxit = 20 , show.plot = TRUE, ...) {

	## computes the node height test of freckleton and harvey ##

	bt <- branching.times(phy);
	bt <- abs(bt - max(bt));
	pics <- abs(pic(d, phy));
	if(show.plot == T) {
		if(log==T) {
			plot(bt, log(pics), main = "node height test", xlab = "branching times", ylab = "log(absolute value of contrasts)", ...);
		} else {
			plot(bt, pics, main = "node height test", xlab = "branching times", ylab = "absolute value of contrasts", ...);			
		}
	}
	if(log == T) {
		if(length(which(pics==0))>0) {
			pics[which(pics==0)] <- 10^-5
		}
		if(regression.type == "lm") return(summary(lm(log(pics)~bt)));
		if(regression.type == "rlm") return(summary(rlm(log(pics)~bt, maxit = rlm.maxit)));	
	} else {
		if(regression.type == "lm") return(summary(lm(pics~bt)));
		if(regression.type == "rlm") return(summary(rlm(pics~bt, maxit = rlm.maxit)));	
		
	}
}







## Internal functions
.pp.getproposal<-function(currentState, psi, bound = FALSE, min = NA, max = NA, prop.type = "sw") {
    ## internal function for mcmc for generating proposals ##
    
    # This function returns a proposal using a sliding window
    # Arguments;
	
	# currentState state <- the currentState state of the parameter;
	# psi <- a tuning parameter - here, some proportion of the standard deviation from the calibration step;
	# min <- the lower bound for the prior distribution;
	# max <- the upper bound for the prior distribution;
	
    # Value;
	# prop <- the proposal;
	if(prop.type == "sw") {
        prop <- currentState + ((runif(1) - 0.5) * psi);
	} else if(prop.type == "mlt") {
		u <- runif(1);
		prop <- currentState * exp((u-0.5)*psi)
		
    }
	
	
	if(bound == TRUE) {
		if(is.na(min) & is.na(max)) stop("you must specify limits if using a bounded proposal!")
		if(!is.na(min) & is.na(max)) max <- Inf;
		if(is.na(min) & !is.na(max)) min <- -Inf;
		if (prop < min) { prop <- min + (min - prop) } # if lower than lower bound, bounce back;
		if (prop > max) { prop <- max - (prop - max) } # if higher than upper bound, bounce back;
	}
    
	if(prop.type == "sw") return(prop);
	if(prop.type == "mlt") return(list(prop = prop, hastings = u))
}

.scale.edge.tree <- function(phy, clade, edge.scalar) {
	node <- .slater.mrca(phy, clade[1], clade[2]);
	edge <- which(phy$edge[,2] == node);
	phy$edge.length[edge] <- phy$edge.length[edge] * edge.scalar
	return(phy);
} 






.reml.lk <- function(phy, d, sigma) {
	
	 ## function to compute the BM likelihood using REML ##

	contr <- pic(d, phy, scaled = F, var.contrasts = T);
	sum.cont <- sum(apply(contr, 1, .cont.foo, sigma =sigma));
	
	n <- length(d)-1;
	
	lk <- -0.5 * ((n * log(2*pi*sigma)) + sum.cont);
	return(lk);
	
}	






.cont.foo <- function(x, sigma) {
	## internal function of reml.lk ##
	return(as.numeric(log(x[2]) + ((as.numeric(x[1])^2) / (sigma*as.numeric(x[2])))));

	}

.pp.acdc.prior <- function (phy, decrease.max = 1e-05, increase.max = 1e+05) {
	# ## This useful function helps set a prior for the EB parameter (described in manuscript text and already used as internal function of fitContinuousMCMC)
 	
    max.bt <- max(heights(phy))
    
    prior.min <- log(decrease.max)/max.bt
    prior.max <- log(increase.max)/max.bt
    if(prior.max == -Inf) prior.max <- 0;
    return(list(min = prior.min, max = prior.max))
}



