
# Arguments:
#	ephy: a bammdata object
# rate: either "speciation", "extinction" or "net diversification". default to "speciation".
#	traits: a vector of trait data, with names corresponding to tips in the ephy object
#	reps: number of permutations to do
#	return.full: include the permutated and observed correlations in the retunred object?
#	method: either spearman,  pearson, mann-whitney, kruskal
#	logrates:  log-transform the rates before analysis? This can only matter for the pearson correlation
#	two.tailed: perform a two tailed statistical test? (in which case it will double the p-value)
# traitorder: for one tail test, specify the direction of correlation ("positive" or "negative") for countinues trait or a string 
#             indicating states of binary trait with increasing speciation rate, separted by comma (e.g., 'A, B'). Currently, this function only
#             perform two-tailed test for categorical data with more than two states.
# Returns:
#	estimate: the mean observed correlation coefficient
#	p.value: the probability that the observed correlation is less than or equal to
#				a sample from the null distribution. If you are doing a one-tailed
#				test on a continuous trait, a small pvalue means that your observed correlation is larger
#				than the null distribution. P-values approaching 1 mean that the observed
#				correlation is more negative than expected.
#	gen: if return.full=TRUE, the vector of generations sampled from the bammdata objected for permutation
#	obs.corr: if return.full=TRUE, a vector of observed correlation coefficients between traits and tip speciation rate for every sampled generation
#	null: if return.full=TRUE, a vector of permuted correlation coefficients for everty sample generation

# ephy<-ed
# traits<-bitrait
# reps<-100
# return.full=T
# method="m"
# logrates=T
# two.tailed=T
# nthreads=6
# traitorder='0,1';
traitDependentBAMM <- function(ephy, traits, reps, rate = 'speciation', return.full = FALSE, method = 'spearman', logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1) {

	if (ephy$type != 'diversification'){
		stop("Function currently supports only bammdata objects from speciationextinction analyses\n");
	}
  
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}
	ratetype.option <- c("speciation", "extinction", "net diversification");
	ratetype <- ratetype.option[grep(paste("^", rate, sep = ''), ratetype.option, ignore.case = TRUE, perl = TRUE)];
  if (length(ratetype) == 0) {
    stop("Rate must be one of 'speciation', 'extinction', or 'net diversification', only the initial letter is needed\n")
  }
  
	if (ratetype == "net diversification" & logrates == TRUE) {
		cat("WARNING: Net diversification might be negative and logged rates would then produce NaNs.\n");
	}
	
	# check if species in ephy and traits match
	if (!identical(sort(names(traits)), sort(ephy$tip.label))) {
		stop('Species names in the bamm-data object and trait data are not identical. You may want to run subtreeBAMM and subset the trait data in order to reduce the two datasets to a common taxon set.')
	}

	method.option <- c("spearman",  "pearson", "mann-whitney", "kruskal");
	method <- method.option[grep(paste("^", method, sep = ''), method.option, ignore.case = TRUE)];
	if (length(method) == 0) {
		stop("method must be one of 'spearman', 'pearson', 'mann-whitney', or 'kruskal', only the initial letter is needed");
	}

	# check if the trait is right class

	if (method == 'spearman' | method == "pearson") {
		if (! is.numeric(traits)){
			cat(paste("selected ", method, ", but the trait is not numeric, converted the trait into a numeric vector\n",sep=''));
			traits <- as.numeric(traits);
		}
	} else if (method == "mann-whitney"| method == "kruskal") {
      
		if (length(unique(traits[! is.na(traits)])) == 1) {
			stop(paste("selected ", method, ", but the trait only has one level\n", sep = ''));
		}
		if (method == "mann-whitney") {
			if (length(unique(traits[! is.na(traits)])) > 2) {
				stop(paste("selected ", method, ", but the trait has more than two levels\n", sep = ''));
			}
		}
	}

	#check if the traitorder is specified
	trait.state <- NA;
	if (two.tailed == FALSE) {
		if (is.na(traitorder)) {
			stop("selected one-tail test, but traitorder is not specified\n");
		}
		if ( method == "kruskal") {
			stop(" currently one-tail test is only available for continous or binary trait");
		}
		if (method == 'spearman' | method == "pearson") {
			direction.option <- c("positive", "negative");
			direction <- direction.option[grep(paste("^",traitorder,sep=''), direction.option, ignore.case = TRUE)];
			if (length(direction) == 0) {
				stop(" for one-tail test with continous trait, traitorder must be either 'positive' or 'negative', only the initial letter is needed");
			} else {
				cat(paste("select one-tailed ", method, " test\nAlternative hypothesis: the trait is ", direction, "ly correlated with speciation rate\n", sep = ''));
			}
		} else {
			traitorder <- gsub(" ", "", traitorder);
			trait.state <- as.character(unlist(strsplit(x = traitorder, split = ",")));
			if (length(trait.state) != 2) {
				stop("please specify the traitorder for binary trait:\nTwo states separated by comma, and the state that is expected to have lower speciation rate first\n");
			} else {
				cat(paste("selected one-tail ", method, " test\nAlternative hypothesis: species with trait ", trait.state[2], " has higher speciation rate than those with trait ", trait.state[1], "\n", sep = ''));
			}
			for (i in trait.state) {
				if (sum(traits == i) == 0) {
					stop(paste("no species with state ", i," \n", sep = ''));
				}
			}
		}
	}
	if (ratetype == 'speciation') {
		tiprates <- ephy$tipLambda;
	} else if (ratetype == "extinction") {
		tiprates <- ephy$tipMu;
	} else {
		tiprates <- lapply(1:length(ephy$tipLambda), function(i) {ephy$tipLambda[[i]] - ephy$tipMu[[i]]});
	}
	tipstates <- ephy$tipStates;
	#tiprates <- tiprates[ephy$tip.label];
	traits <- traits[ephy$tip.label];
	stat.mu <- 0;
	if (method == "mann-whitney") {
		trait.stat.count <- table(traits);
		trait.stat.count <- trait.stat.count[! is.na(names(trait.stat.count))];
		stat.mu <- prod(trait.stat.count) / 2;
	}

	if (logrates) {
		tiprates <- lapply(1:length(tiprates), function(x){ log(tiprates[[x]]) });
	}

	#randomly sample generations from BAMM posterior
	gen <- sample(1:length(tiprates), size = reps, replace = TRUE);

	gen.tiprates <- list();
	for (l in 1:length(gen)) {
		gen.tiprates[[l]] <- data.frame(rates = tiprates[[gen[l]]], states = tipstates[[gen[l]]], stringsAsFactors = FALSE);
	}
  
	rm("tiprates", "tipstates");
	permute_tiprates <- function(m) {
		tt <- m$states;
		tlam <- m$rates;
		index <- unique(tt);
		lvec <- numeric(length(index));
		for (k in 1:length(index)) {
			lvec[k] <- tlam[tt == index[k]][1];
		}
		new_index <- sample(index, size = length(index), replace = FALSE);
		x <- rep(0,length(tt));
		for (xx in 1:length(index)) {
			x[which(tt == index[xx])] <- lvec[which(index == new_index[xx])];
		}
		x;   
	}
	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads);
		p.gen.tiprates <- parallel::parLapply(cl, gen.tiprates, permute_tiprates);
		parallel::stopCluster(cl);
	} else {
		p.gen.tiprates <- lapply(gen.tiprates, permute_tiprates);
	}
	xgen.tiprates <- list();
	for (l in 1:length(gen)) {
		xgen.tiprates[[l]] <- gen.tiprates[[l]]$rates;
	}
	gen.tiprates <- xgen.tiprates; rm("xgen.tiprates");

	cortest <- function(rates, traits, method) {
	    if (sd(rates, na.rm = TRUE) == 0) {
	    	return(0);
	    } else {
			return(cor.test(rates, traits, method = method, exact = FALSE)$estimate);  
	    }
	}
	manntest <- function(rates, traits, two.tailed, trait.state) {
		if (two.tailed) {
			return(wilcox.test(rates ~ traits, exact = FALSE)$statistic);
		} else {
			return(wilcox.test(rates[which(traits == trait.state[2])], rates[which(traits == trait.state[1])], exact = FALSE)$statistic);
		}
	}

	kruskaltest <- function(rates, traits) {
		return(kruskal.test(rates ~ traits)$statistic);
	}
 

	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads);
		if (method == 'spearman' | method == "pearson") {
			obs <- parallel::parLapply(cl, gen.tiprates,cortest, traits, method);
			permu <- parallel::parLapply(cl, p.gen.tiprates, cortest, traits, method);
		} else if (method == "mann-whitney") {
			obs <- parallel::parLapply(cl, gen.tiprates, manntest, traits, two.tailed, trait.state);
			permu <- parallel::parLapply(cl, p.gen.tiprates, manntest, traits, two.tailed, trait.state);
		} else {
			obs <- parallel::parLapply(cl, gen.tiprates, kruskaltest, traits);
			permu <- parallel::parLapply(cl, p.gen.tiprates, kruskaltest, traits);
		}
		parallel::stopCluster(cl);
	} else {
		if (method == 'spearman' | method == "pearson") {
			obs <- lapply(gen.tiprates, cortest, traits, method);
			permu <- lapply(p.gen.tiprates, cortest, traits, method);
		} else if (method == "mann-whitney") {      
			obs <- lapply(gen.tiprates, manntest, traits, two.tailed, trait.state);
			permu <- lapply(p.gen.tiprates, manntest, traits, two.tailed, trait.state);
		} else {    
			obs <- lapply(gen.tiprates, kruskaltest, traits);
			permu <- lapply(p.gen.tiprates, kruskaltest, traits);    
		}
	}
	obs <- unlist(obs);
	permu <- unlist(permu);
	obs <- obs - stat.mu;
	permu <- permu - stat.mu;

  

	if (two.tailed) {
		pval <- sum(abs(obs) <= abs(permu)) / length(permu);
	} else {
		if (method == "spearman" | method == "pearson") {
			if (direction == 'positive') {
				pval <- sum(obs <= permu) / length(permu);
			} else {
				pval <- sum(obs >= permu) / length(permu);
			}
		} else {
			pval <- sum(obs <= permu) / length(permu);
		}
	}
	if (method == "spearman" | method == "pearson") {
		obj <- list(estimate = mean(as.numeric(obs)), p.value = pval, method = method, two.tailed = two.tailed);
	} else {
		if (ratetype == 'speciation') {
			ave.tiprate <- getTipRates(ephy)$lambda.avg;
		} else if (ratetype == 'extinction') {
			ave.tiprate <- getTipRates(ephy)$mu.avg;
		} else {
			ave.tiprate <- getTipRates(ephy)$lambda.avg - getTipRates(ephy)$mu.avg;
		}
		l <- lapply(unique(traits[! is.na(traits)]), function(x) {
			median(ave.tiprate[which(traits == x)], na.rm = TRUE);
		});
		names(l) <- as.character(unique(traits[! is.na(traits)]));
		obj <- list(estimate = l, p.value = pval, method = method, two.tailed = two.tailed);
	}
	obj$rate <- ratetype;
	if (return.full) {
		obj$obs.corr <- as.numeric(obs);
		obj$gen <- gen;
		obj$null <- as.numeric(permu);
	}

	return(obj);
}

