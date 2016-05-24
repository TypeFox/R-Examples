.laser.gsr<-
#gsr<-
function (Source, Search, Replace, char = FALSE)
{
    if (length(Search) != length(Replace))
        stop("Search and Replace Must Have Equal Number of Items\n")
    Changed <- as.character(Source)
    if (char == FALSE) {
        for (i in 1:length(Search)) {
            Changed <- replace(Changed, Changed == Search[i],
                Replace[i])
        }
    }
    else if (char == TRUE) {
        for (i in 1:length(Search)) {
            Changed <- replace(Changed, Changed == Search[i],
                paste(Replace[i]))
        }
    }
    Changed
}


.laser.gettipdata<-
#getTipdata<-
function (tipdata, phy)
{
    if (is.data.frame(tipdata)) {
        x <- as.vector(tipdata[, 1])
        names(x) <- row.names(tipdata)
    }
    else {
        x <- tipdata
    }
    if (class(phy) != "phylo")
        stop("object \"phy\" is not of class \"phylo\"")
    if (is.null(phy$edge.length))
        stop("your tree has no branch lengths: invalid input")
    tmp <- phy$edge
    nb.tip <- length(phy$tip.label)
    if (phy$Nnode != nb.tip - 1)
        stop("\"phy\" is not fully dichotomous")
    if (length(x) != nb.tip)
        stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x)))
        stop("method can't be used with missing data...")
    phenotype <- as.numeric(rep(NA, nb.tip + phy$Nnode))
    names(phenotype) <- 1:max(phy$edge)
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    }
    else {
        if (!any(is.na(match(names(x), phy$tip.label))))
            phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels\ndid not match\n")
        }
    }
    phy$phenotype <- rep(NA, (nb.tip + phy$Nnode - 1))
    for (i in 1:length(phenotype)) {
        phy$phenotype[as.numeric(phy$edge[, 2]) == names(phenotype[i])] <- phenotype[i]
    }
    return(phy)
}

.laser.splitedgematrix<-
#splitEdgeMatrix<-
function (phy, node)
{
    x <- branching.times(phy)
    rootnode <- length(phy$tip.label) + 1
    phy$tag <- rep(1, nrow(phy$edge))
    if (node >= rootnode) {
        node.desc <- node
        pos <- 1
        phy$tag[phy$edge[, 1] == node.desc[1]] <- 2
        while (pos != (length(node.desc) + 1)) {
            temp <- .get.desc.of.node(node.desc[pos], phy)
            temp <- temp[temp > rootnode]
            for (k in 1:length(temp)) {
                phy$tag[phy$edge[, 1] == temp[k]] <- 2
            }
            node.desc <- c(node.desc, temp)
            pos <- pos + 1
        }
    }
    else if (node > 0)
        phy$tag[phy$edge[, 2] == node] <- 2
    z <- cbind(phy$edge, .laser.gsr(phy$edge[, 1], names(x), x), phy$edge.length,
        phy$phenotype, phy$tag)
    z <- matrix(as.numeric(z), dim(z))
    z <- as.data.frame(z)
    return(z)
}

.laser.getlambda<-
#getLambda.internal<-
function (zmat, rootnode, rbounds, para = 0.01, eps, combined = TRUE)
{
    int <- zmat[zmat[, 2] > rootnode, ]
    term <- zmat[zmat[, 2] < rootnode, ]
    nint <- nrow(int)
    nterm <- nrow(term)
    betaF <- function(r, t1) {
        xf <- (exp(r * t1) - 1)/(exp(r * t1) - eps)
        xf
    }
    Lfunc_tax <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm,
            5] - 1) * log(betaF(r, term[1:nterm, 4]))))
    }
    Lfunc_phy <- function(p) {
        r <- p
        (nint * log(r) - r * sum(int[1:nint, 4]) - sum(log(1 -
            (eps * exp(-r * int[1:nint, 3])))))
    }
    Lfunc_comb <- function(p) {
        r <- p
        (sum(log(1 - betaF(r, term[1:nterm, 4]))) + sum((term[1:nterm,
            5] - 1) * log(betaF(r, term[1:nterm, 4]))) + nint *
            log(r) - r * sum(int[1:nint, 4]) - sum(log(1 - (eps *
            exp(-r * int[1:nint, 3])))))
    }
    res <- list()
    if (combined == TRUE) {
        if (nrow(int) == 0)
            tempres <- optimize(Lfunc_tax, interval = rbounds,
                maximum = TRUE)
        else tempres <- optimize(Lfunc_comb, interval = rbounds,
            maximum = TRUE)
    }
    else {
        tempres <- optimize(Lfunc_tax, interval = rbounds, maximum = TRUE)
    }
    res$LH <- tempres$objective
    res$lambda <- tempres$maximum/(1 - eps)
    res$r <- tempres$maximum
    res$eps <- eps
    res <- as.data.frame(res)
    return(res)
}

.laser.fitNDR_1rate<-
#fitNDR_1rate<-
function (phy, eps = 0, rbounds = c(1e-04, 0.5), combined = TRUE)
{
    z <- .laser.splitedgematrix(phy, phy$Nnode)
    r1 <- .laser.getlambda(z, rootnode = (length(phy$tip.label) +
        1), rbounds = rbounds, para = 0.01, eps, combined = combined)
    res <- list()
    res$LH <- r1$LH
    res$aic <- (-2 * r1$LH) + 2
    res$r <- r1$r
    res$lambda <- r1$lambda
    res$eps <- eps
    res <- as.data.frame(res)
    return(res)
}

sim.mecca <-
#MECCAsim <-
function(phy, richness, cladeAges, model, prop, makeNewTipTrees = TRUE, mytiptrees = NULL, hotBranches = NULL) {

	root.age = max(branching.times(phy));

	if(makeNewTipTrees == FALSE & is.null(mytiptrees)) stop("there are no tip trees provided");
	cladeAncStates <- .mecca.nodesim(phy, model, prop, hotBranches);

	Svar <- numeric(length(richness));
	Smean <- numeric(length(richness));

	if(makeNewTipTrees == TRUE) {
		treelist<-list();
		treelist[which(richness==1)] <- cladeAges[richness==1];
		for(i in which(richness>1)) {
			treelist[[i]]<- .mecca.fasttreesim(n=richness[i], lambda = prop$birth, mu=prop$death, rho = 1, origin = cladeAges[i]); ## generate the tree  under prior drawn lambda
		}
	} else {
		treelist <- mytiptrees;
	}


	for(i in 1:length(treelist)) {

		if(!is.numeric(treelist[[i]])) {
			if(!model == "twoRate") dat <- .mecca.tipsim(treelist[[i]], model, prop, cladeAncStates[i], root.age);
			if(model == "twoRate") {

				if (i %in% hotBranches) { Branch = "hot" } else { Branch <- "not" }
				dat <- .mecca.tipsim(phy = treelist[[i]], model = model, prop = prop, cladeAncStates[i], BranchState = Branch);
			}

 			Svar[i]<-var(dat); # here is the variance of the tip data for a character evolved under the prior for sigma. we want to know if it provides a good approximation to our data, so...
			Smean[i]<-mean(dat);

		} else {
			if(model == "twoRate") {
				if (i %in% hotBranches) { Branch = "hot" } else { Branch <- "not" }

					dat <- .mecca.singletiptrait(model, prop, cladeAges[i], cladeAncStates[i], root.age, BranchState = Branch);
				}
				else {
				dat <- .mecca.singletiptrait(model, prop, cladeAges[i], cladeAncStates[i], root.age);
			}


			Smean[i] <- dat;
			Svar[i] <- 0;
		}

	}
	return(list("trees" = treelist, "Smean" = Smean, "Svar" = Svar));
	}

.mecca.makecaloutput<-
#MakeCalOutput <-
function(model, Ncalibrations, Nclades) {


	if(model == "BM") k <- 5;
	if(model == "Trend" | model == "twoRate") k <- 6;

	output <- matrix(NA, nrow = Ncalibrations, ncol = k + 2*Nclades);


	if(model == "BM") colnames(output) <- c("birth", "death","logLk", "sigmasq","root", paste("var", seq(1, Nclades), sep = ""), paste("Mean", seq(1, Nclades), sep = ""));

		if(model == "Trend") colnames(output) <- c("birth", "death","logLk", "sigmasq","root","mu", paste("var", seq(1, Nclades), sep = ""), paste("Mean", seq(1, Nclades), sep = ""));

	if(model == "twoRate") colnames(output) <- c("birth", "death","logLk", "sigmasq1","sigmasq2","root", paste("var", seq(1, Nclades), sep = ""), paste("Mean", seq(1, Nclades), sep = ""));


	return(output);

}

.mecca.proposal<-
#Mecca.Proposal <-
function(model, trait.params, trait.widths, prior.list, SigmaBounds, scale, ngen) {


	if(model == "BM") {

		if (ngen %% 2 == 1) {

			trait.params$sigmasq <- .mecca.getproposal(trait.params$sigmasq, trait.widths$sigmasq * scale, min = SigmaBounds[1], max = SigmaBounds[2])

		} else {

			trait.params$root <- .mecca.getproposal(trait.params$root, trait.widths$root * scale, min = prior.list$priorMean[1], max = prior.list$priorMean[2])
		}
	} # end BM


	if(model == "twoRate") {

		if (ngen %% 2 == 1) {

			trait.params$sigmasq1 <- .mecca.getproposal(trait.params$sigmasq1, trait.widths$sigmasq1 * scale, min = SigmaBounds[1], max = SigmaBounds[2])
			trait.params$sigmasq2 <- .mecca.getproposal(trait.params$sigmasq2, trait.widths$sigmasq2 * scale, min = SigmaBounds[1], max = SigmaBounds[2])

		} else {

		trait.params$root <- .mecca.getproposal(trait.params$root, trait.widths$root * scale, min = prior.list$priorMean[1], max = prior.list$priorMean[2])
		}

	}

	if(model == "Trend"){

		if (ngen %% 2 == 1) {

			trait.params$sigmasq <- .mecca.getproposal(as.numeric(trait.params$sigmasq), as.numeric(trait.widths$sigmasq) * scale, min = SigmaBounds[1], max = SigmaBounds[2])
			trait.params$mu <- .mecca.getproposal(trait.params$mu, trait.widths$mu * scale, min = prior.list$priorMu[1], max = prior.list$priorMu[2])

		} else {

		trait.params$root <- .mecca.getproposal(trait.params$root, trait.widths$root * scale, min = prior.list$priorMean[1], max = prior.list$priorMean[2])
		}

	}


	return(trait.params)
}

.mecca.calproposal <-
#MeccaCalProposal <-
function(model, prior.list, sigmaPriorType, SigmaBounds, rootPriorType, thetaB, thetaD, propWidth) {

		while(1) {
			birth = .mecca.getproposal(thetaB,propWidth, 0, Inf); # prior for lambda
			death = .mecca.getproposal(thetaD, propWidth, 0, Inf); # prior for mu
			if(birth > death) (break);

		}
		sigma = .mecca.priordrawsig(prior.list$priorSigma[1], prior.list$priorSigma[2], sigmaPriorType, SigmaBounds)   # Prior for sigma^2
		root = .mecca.priordrawmean(prior.list$priorMean[1], prior.list$priorMean[2], rootPriorType); # prior for the root state
		if(model == "BM") return(prop = list(birth = birth, death = death, sigmasq = sigma, root = root));

		if(model == "Trend") return(prop = list(birth = birth, death = death, sigmasq = sigma, root = root, mu = runif(1, min = prior.list$priorMu[1],max = prior.list$priorMu[2])));

		if(model == "twoRate") return(prop = list(birth = birth, death = death, sigmasq1 = sigma, sigmasq2 =.mecca.priordrawsig(prior.list$priorSigma[1], prior.list$priorSigma[2], sigmaPriorType, SigmaBounds), root = root));

		### do bounded here later

}

startingpt.mecca <-
#MeccaStartValues <-
function(calibrationOutput, phy, cladeMean, cladeVariance, tolerance = 0.01, plsComponents, BoxCox = TRUE) {

	if(BoxCox == TRUE & length(calibrationOutput) < 6) stop("there are no parameters in the calibration output to do boxcox transformation. Try rerunning calibration with BoxCox = TRUE");

	m <- match(phy$tip.label, names(cladeMean));
	cladeMean <- cladeMean[m];
	m <- match(phy$tip.label, names(cladeVariance));
	cladeVariance <- cladeVariance[m];

	if( length(which(cladeVariance == 0))>0) {
		cladeVariance <- cladeVariance[-which(cladeVariance==0)];
	}


	obs <- c(cladeVariance, cladeMean);
	names(obs) <- c(paste(names(cladeVariance), "var", sep = "_"), paste(names(cladeMean), "mean", sep = "_"));
	stdobs <- obs; # observed stats for standardization and pls

	bdcal <- as.matrix(calibrationOutput$diversification);
	bmcal <- as.matrix(calibrationOutput$trait);

	startingBirth <- as.numeric(bdcal[1,1]);
	startingDeath <- as.numeric(bdcal[1,2]);


	if(BoxCox == TRUE) {
		stdz <- calibrationOutput$stdz;
		lambda <- calibrationOutput$lambda;
		GM <- calibrationOutput$GM;
		boxcox <- calibrationOutput$BoxCox;
	}


	K <- ncol(plsComponents); # gives the number of PLS components
	Nsims <- nrow(bdcal)

	if(BoxCox == TRUE) {

		for(i in 1:length(stdobs)) {
			stdobs[i] <- 1 + (stdobs[i] - stdz$min[i]) / (stdz$max[i] - stdz$min[i]);
			stdobs[i] <-(stdobs[i]^lambda[i] - 1) / (lambda[i] * GM[i]^(lambda[i]-1));
			stdobs[i] <- (stdobs[i] - boxcox$BCmeans[i])/boxcox$BCstd[i];

			}
		}

	obsPLS <- numeric(K);

	for(i in 1:length(obsPLS)) {
		obsPLS[i] <- sum(stdobs * plsComponents[, i]);
	} # creates the pls asjusted oberved data

	plsSim <- matrix(NA, nrow = Nsims, ncol = K);
	params <- bmcal[ , 1: K];
	stats <- bmcal[ , - seq(1, K)];

	for(i in 1:Nsims) {

		plsSim[i, ] <- .mecca.extractpls(stats[i,], plsComponents, K); ## this needs adjusting for different numbers of pls components

	} ## now this holds the pls adjusted distances

	dmat <- numeric(Nsims);

	for (i in 1:Nsims) {

		dmat[i] <- dist(rbind(obsPLS, plsSim[i, ]));
		} # computes the euclidean distances for the simulations

	pdmat <- cbind(params, dmat);
	bmretained <- pdmat[order(as.data.frame(pdmat)$dmat),][seq(1, Nsims * tolerance) , ];

	## and now we compute the starting values (mean) and proposal width (sd) for the MCMC ##

	tuning <- matrix(data = NA, nrow = K, ncol = 2);
	rownames(tuning) <- colnames(bmretained)[1:K];
	colnames(tuning) <- c("starting", "width");

	for(i in 1:K) {

	tuning[i, 1] <- as.numeric(bmretained[1,i]);
	tuning[i, 2] <- sd(bmretained[, i]);

	}

    if(BoxCox == TRUE) {
		return(list("tuning" = tuning, "startingBirth" = startingBirth, "startingDeath" = startingDeath, "dcrit" = max(bmretained[,K+1]), "obsTraits" = obs, "plsObserved" = obsPLS, "plsLoadings" = plsComponents, "stdz" = stdz, "lambda" = lambda, "GM" = GM, "BoxCox" = boxcox));
		}

	if(BoxCox == FALSE) {
		return(list("tuning" = tuning, "startingBirth" = startingBirth, "startingDeath" = startingDeath, "dcrit" = max(bmretained[,K+1]), "obsTraits" = obs, "plsObserved" = obsPLS, "plsLoadings" = plsComponents));
		}
	}

.mecca.nodestatestworate <-
#NodeStatesTwoRate <-
function(phy, prop, hotBranches) {


	cbind(phy$edge, phy$edge.length) -> edge.mat
 	edge.mat <- edge.mat[- which(edge.mat[ ,2] <= length(phy$tip.label)), ];
	rootNode <- length(phy$tip.label) + 1;
	rootState <- prop$root;

	currentNode <- rootNode;
	nodeStates <- numeric(length(phy$tip.label) - 1);
	names(nodeStates) <- seq(length(phy$tip.label) + 1, (length(phy$tip.label) * 2) - 1, 1);
	nodeStates[1] <- rootState

	for (i in 2:length(nodeStates)) {
		currentNode <- currentNode + 1;
		parent <- edge.mat[which(edge.mat[ ,2] == currentNode), 1];
		parentValue <- nodeStates[which(names(nodeStates) == parent)];
		if(currentNode %in% hotBranches) {
			sigma <- prop$sigmasq2;
			} else {
			sigma <- prop$sigmasq1
			}

		nodeStates[i] <- rnorm(1, parentValue, sqrt(sigma * edge.mat[which(edge.mat[ ,2] == currentNode), 3]));
		}

	return(nodeStates);

	}

.mecca.normalpriorratio <-
#NormalpriorRatio <-
#UNUSED
function(currentState, proposedState, mean, sd) {

	h <- dnorm(proposedState, mean, sd) /  dnorm(currentState, mean, sd); # compute the prior density ratio

	if(h > 1) { currentState <- proposedState } # if the prior ratio is greater than one, accept the proposal

	if (h < 1) {

		rnum <- runif(1);

		if (h < rnum) {
			currentState <- proposedState;
			}else{
				currentState <- currentState; # if the prior ratio is less than one, accept proposal with probability h
				}
	 # return the new state
		}
		return(currentState);
	}

.mecca.tipsim<-
#TipSim <-
function(phy, model, prop, node.state, root.age, BranchState = "not") {


	if(model == "BM") {

		dat <- .mecca.fastbm(phy, node.state, prop$sigmasq, mu = 0)$tipStates;

	}


	if(model == "Trend") {

		dat <- .mecca.fastbm(phy, node.state, prop$sigmasq, prop$mu)$tipStates;


	}

	if(model == "twoRate") {
		if(BranchState == "not") {
			dat <- .mecca.fastbm(phy, node.state, prop$sigmasq1, mu = 0)$tipStates;
		} else if(BranchState == "hot") {
			dat <- .mecca.fastbm(phy, node.state, prop$sigmasq2, mu = 0)$tipStates;

		}
			}

	return(dat);
}

.mecca.acceptance <-
#acceptance <-
function(trait.params, prop.params, distribution, params) {

	if(distribution == "uniform") {
		return("accept");
	}

	if(distribution == "normal") {
		curr.rates <- trait.params[grep("sigma", colnames(trait.params))];
		prop.rates <- prop.params[grep("sigma", colnames(prop.params))];

		k <- length(curr.rates);
		priorden <- numeric(k)
		for(i in 1:k) {
		priorden[i] <- dnorm(as.numeric(prop.rates[i]), mean = params[1], sd = params[2]) / dnorm(as.numeric(curr.rates[i]), mean = params[1], sd = params[2]);
		}
		priorden <- prod(priorden);
		p <- runif(1);

		if(priorden >= p) {

			return("accept")
			} else {

				return("reject")
				}

		}

	}

.mecca.allbranches <-
#allBranches <-
function(phy, node) {
	node <- as.numeric(node);
	n <- length(phy$tip.label);
	if(node <= n) return(node);


	l <- .get.desc.of.node(node, phy);
	for(j in l) {
		if(j > n) l<-c(l, .mecca.allbranches(phy, j));
	}
	return(l);
}

.mecca.boxcox <-
#boxcoxTransform <-
function(summaries, stdz, lambda, GM, boxcox) {

	for(i in 1:length(summaries)) {
	summaries[i] <- 1 + (summaries[i] - stdz$min[i]) / (stdz$max[i] - stdz$min[i]);
	summaries[i] <-(summaries[i]^lambda[i] - 1) / (lambda[i] * GM[i]^(lambda[i]-1));
	summaries[i] <- (summaries[i] - boxcox$BCmeans[i])/boxcox$BCstd[i];

			}

	return(summaries);

	}

#mecca.dat<-mecca.lm<-NA
calibrate.mecca <-
#calibrateMECCA <-
function(phy, richness, model = c("BM", "Trend", "twoRate"), prior.list = list(priorSigma = c(-4.961845, 4.247066), priorMean = c(-10, 10)), Ncalibrations = 10000, sigmaPriorType = "uniform", rootPriorType = "uniform", divSampleFreq = 0, initPars = "ML", propWidth = 0.1, SigmaBounds = c(-4.961845, 4.247066), hotclade = NULL, BOXCOX = TRUE, optimRange =c(-1000, 10)) {

## add requirement for pls
    ## require(pls)


##Arguments
	# phy <- phylogenetic tree
	# richness <- a vector of species richness values corresponding to tips in the tree
	# cladeMean <- a vector of clade mean trait values corresponding to tips in the tree
	# cladeVariance <- a vector of clade variance trait values corresponding to tips in the tree
	# Ncalibrations <- number of calibration steps - by default 10000 following Wegmann et al. 2009
	# priors <- concatenated vectors giving the min and max for the priors on the four parameters. NB if the root state prior is normal the values are assumed to be mean and sd instead
	# rootPriorType <- the root prior distribution - the default is uniform but a normal can be specified at present. Other prior types should probably be added
	# propWidth - only used for the b/d parameters here

## Values
	# a list containing calibration output for MECCA
	## these first three vectors will hold the distances between observed and simulated data ##


	## house keeping block - check everything is in the correct order, get rid of node labels and generate a vector of clade age values ##

	name.check(phy, richness)->nc;
	if(!nc == "OK") { stop("names in tree and data do not match") }
	if(!sum(is.na(richness) == 0)) { stop("richness contains missing values") }
	Nclades <- length(richness);
	phy$node.label <- NULL;
	hotBranches <- NULL;

	match(phy$tip.label, names(richness)) -> m;
	richness[m] -> richness;
	cladeAges <-  .mecca.cladeage(phy);

	if(model =="BM") k <- 2; if(model == "ACDC"| model == "Trend"| model == "twoRate" | model == "SSP") k <- 3; if(model == "OU") k <- 4;


	if(model=="twoRate") {

		if(is.null(hotclade)) { stop("if using a two rate model you need to specify the hot clade") }

		hotBranches <- .mecca.gethotbranches(phy, hotclade[1], hotclade[2]);

	}

	if(model == "Trend") {
		if (is.null(prior.list$priorMu)) {prior.list$priorMu = c(-0.5,0.5); warning("No mu prior specified for Trend model. Using default prior of -0.5:0.5") }
		}

	output <- .mecca.makecaloutput(model, Ncalibrations, Nclades);


	## end house keeping block ##

	if (length(initPars) > 1) {
		thetaB <- initPars[1]; thetaD <- initPars[2];
		} else if (initPars == "ML") {
			  ## generate starting states
				mleBD <- .mecca.getdivMLE(phy, richness); thetaB <- mleBD$lamda; thetaD <- mleBD$mu
		}

	Lk <- .mecca.logLP(phy, thetaB, thetaD) + .mecca.logLT(phy, richness, thetaB, thetaD);

	pb <- txtProgressBar(min = 0, max = Ncalibrations, char = "-", style = 3); #cat("Calibrating...\n\n")

	for(ncal in 1:Ncalibrations) {     ## begin the calibration ##

		cvar <- numeric(length(richness)); # will hold simulated calibration variances
		cmeans <- numeric(length(richness)); # will hold simulated clade mean values

		prop <- .mecca.calproposal(model, prior.list, sigmaPriorType, SigmaBounds, rootPriorType, thetaB, thetaD, propWidth)

		## we can now compute the taxonomic/phylogenetic likelihood of the tree data under the b/d priors

		propLk <- .mecca.logLP(phy, prop$birth, prop$death) + .mecca.logLT(phy, richness, prop$birth, prop$death);

		if(propLk== - Inf) {
			LKratio <- 0 ;

			} else {

			LKratio <- exp(propLk - Lk);
			}

			if(LKratio >= 1) {

				thetaB <- prop$birth; thetaD <- prop$death; Lk <- propLk; ## accept and update

			}

			if (LKratio < 1) {

			p <- runif(1);
			if (p <= LKratio) {
				thetaB <- prop$birth; thetaD <- prop$death; Lk <- propLk;## accept and update

			} else {

				prop$birth <- thetaB; prop$death <- thetaD; ## reject and modify b/d proposals for simulation

			}

		}

		##
		if(ncal  == 1 | divSampleFreq == 0 |ncal  %% divSampleFreq == 0) {
			calSim <- sim.mecca(phy = phy, richness = richness, cladeAges = cladeAges, model = model, prop = prop, makeNewTipTrees = TRUE, mytiptrees = NULL, hotBranches = hotBranches);
			tipTrees <- calSim$trees;

			} else {

			calSim <- .mecca.generatetipsummaries(phy, model, tipTrees, prop = prop, hotBranches = hotBranches, richness=richness);

		}

  		output[ncal, ] <- .mecca.paramsample(model, prop, Lk, calSim);
  		setTxtProgressBar(pb, ncal); # move the progress bar
		} ## end the calibration loop

	close(pb); # close the progress bar

  	div.cal <- output[, c(1,2,3)]; ## Birth death Parameters and stats
	trait.cal <- output[, -c(1,2,3)]; ## BM parameters and stats
	if(length(which(apply(trait.cal, 2, function(x) sum(x) == 0) == TRUE))>0) {
		trait.cal <- trait.cal[, -which(apply(trait.cal, 2, function(x) sum(x) == 0) == TRUE)];
	}

	if(BOXCOX == TRUE) {
		trait.cal <- as.data.frame(trait.cal)
		stat <- trait.cal[, -(1:k)];
		param <- trait.cal[ ,1:k];

		myMax<-array(0, dim = length(stat));
		myMin<-array(0, dim = length(stat));
		lambda<-array(0, dim = length(stat));
		myGM<-array(0, dim = length(stat));
		myBCMeans <- array(0, dim = length(stat));
		myBCSDs <- array(0, dim = length(stat));

		## preliminary transform
		for(i in 1:length(stat)){
			myMax[i]<-max(stat[,i]);
			myMin[i]<-min(stat[,i]);
			stat[ ,i] <- 1+(stat[ ,i] - myMin[i])/(myMax[i]-myMin[ i]);

			}

        mecca.lm <- mecca.env <- NULL

		lmreturn=function(stati, param, ...){
		  ## FIXME: unable to figure out how to make .mecca.maxboxcox work without mecca.lm and mecca.env exported to global environment (prompting a NOTE from R CMD check -- JME)
		  ## fixed by GJS following RF, but check...

		  mecca.dat=as.data.frame(cbind(y=stati, param))
		  mecca.lm<-as.formula(mecca.dat)
		  mecca.env<-environment(mecca.lm)
		  #optimize(.mecca.maxboxcox, ..., linmod = lm(mecca.lm, mecca.env$object))$maximum
		  optimize(.mecca.maxboxcox, ..., env=mecca.env)$maximum
		}
		for(i in 1:length(stat)) {
		  #JME d <<- cbind(stat[,i], param);
		  #JME d <<- as.data.frame(d)

		  lambda[i] = lmreturn(stat[,i], param, interval = optimRange, maximum = TRUE)
		  #JME mylm <<- lm(as.formula(d), data=d);
		  print(paste(colnames(stat)[i], lambda[i]))
		}

			for(i in 1:length(stat)) {
			myGM[i]<-exp(mean(log(stat[,i])));

			stat[,i]<-(stat[ ,i]^lambda[i] - 1) / (lambda[i] * myGM[i]^(lambda[i]-1));
			myBCSDs[i] <- sd(stat[,i]);
			myBCMeans[i] <- mean(stat[,i]);
			stat[,i]<-(stat[,i] -myBCMeans[i])/myBCSDs[i];
			}




			res=list("diversification" = div.cal, "trait" = cbind(param, stat), "stdz" = list("min" = myMin, "max" = myMax), "lambda" = lambda, "GM"=myGM, "BoxCox" = list("BCmeans" = myBCMeans, "BCstd" = myBCSDs));

		}


	if(BOXCOX == FALSE) {
		res=list("diversification" = div.cal, "trait" = trait.cal);

    }
    class(res)=c("mecca", "calibration", class(res))
    return(res)


} # end


.mecca.maxboxcox <-
  #maxboxcox <-
  function(lambda, env) {
    myboxcox<-boxcox(env$rval, lambda=lambda, data=env$object, interp=T, eps=1/50, plotit = FALSE);
    return(myboxcox$y[1]);
  }

.mecca.cladeage <-
#cladeAge <-
function(phy) {
# a function to pull tip lineage ages from a phylogenetic tree

# Arguments
	# a phylogenetic tree

#Values
	# a named vector of ages for terminal lineages

			which(phy$edge[ ,2] <= length(phy$tip.label)) -> tips;
			cladeAges <- phy$edge.length[tips];
			tipNumbers <- phy$edge[tips, 2];
			names(cladeAges) <- phy$tip.label[tipNumbers];
			match(phy$tip.label, names(cladeAges)) -> m;
			cladeAges[m]->cladeAges;
			return(cladeAges)

}

.mecca.cladestates <-
#cladeStates <-
function(nodeStates, phy) {

	cladeAncStates <- numeric(length(phy$tip.label));

	for(k in 1:length(cladeAncStates)) {
		parent <- .mecca.parentnode(phy, k);
		cladeAncStates[k] <- nodeStates[which(names(nodeStates) == parent)];

				}
	return(cladeAncStates);

}

.mecca.extractpls <-
#extractpls <-
function(p, plsdat, Ncomp) {
	x <- numeric(Ncomp);
	for(i in 1:Ncomp) {
		x[i] <- sum(p * plsdat[, i])
		}
	return(x)
	}

.mecca.fastbm <-
#fastBM2 <-
function(tree, alpha, sig2, mu) {

	n <- length(tree$edge.length);
	x <- rnorm(n=n, mean = mu * tree$edge.length, sd=sqrt(sig2*tree$edge.length))
	y <- array(alpha, dim=n+1);

	for(i in 1:n){
		y[tree$edge[i,2]]<-y[tree$edge[i,1]]+x[i];
	}
	nodeStates<-y[(length(tree$tip.label)+1):(n+1)];
	names(nodeStates)<-as.character((length(tree$tip.label)+1):(n+1));
	tipStates<-y[1:length(tree$tip.label)];
	names(tipStates)<-tree$tip.label;

	return(list(nodeStates=nodeStates, tipStates=tipStates));
}

.mecca.fasttreesim <-
#fastTreeSim <-
function (n, lambda, mu, rho, origin) {

## TreeSim Function by Tanja Standler, modified by Daniel Wegmann for speed
## n  <- number of desired species in the tree
## lambda <- the speciation, or birth rate
## mu <- the extinction, or death rate
## rho <- each tip is included in the final tree with probability rho. should be set to 1.
## origin <- the time of origin, or age of the clade

    edge<-matrix(data=NA, nrow=2*n-1, ncol=2);
    edge[1,]<- c(-1, -2);
    leaves <- c(-2)
    timecreation <- array(data=0, dim=2*n);
    time <- 0
    maxspecies <- -2
    edge.length <- array(data=0, dim=2*n-1);

    specevents = array(dim=n)
    specevents[1]<-origin;
    r<-runif(n-1, 0,1);

    if (lambda > mu) {
      lamb1 <- rho * lambda
      mu1 <- mu - lambda * (1 - rho)
      AA<-exp((-lamb1 + mu1) * origin);
      YY<-lamb1 - mu1 * AA;
      XX<-(1 - AA) * r;
      specevents[2:n] <- 1/(lamb1 - mu1) * log((YY - mu1 * XX)/(YY - lamb1 * XX));
    } else {
      specevents[2:n] <- (origin * r)/(1 + lambda * rho * origin * (1 + r));
    }

    specevents <- sort(specevents, decreasing = TRUE)
    for (index in 2:n) {
        time = time + (specevents[index - 1] - specevents[index]);
        species <- sample(leaves, 1)
        del <- which(leaves == species)

        edge.length[which(edge[,2] == species)] <- time - timecreation[-species]

	edge[2*index-2,]<-c(species, maxspecies - 1);
	edge[2*index-1,]<-c(species, maxspecies - 2);

	leaves <- c(leaves[-del], maxspecies - 1, maxspecies - 2)
        maxspecies <- maxspecies - 2
        timecreation[2*index-1]<-time;
	timecreation[2*index]<-time;
    }

    m<-match(leaves, edge[,2]);
    edge.length[m]<-specevents[1]-timecreation[-edge[m,2]];


    nodes <- (length(leaves)) * 2
    leaf = 1
    interior = length(leaves) + 1
    if (nodes != 2) {
         for (j in 1:nodes) {
            if (sum(match(leaves, -j, 0)) == 0) {
                posvalue <- interior;
                interior <- interior + 1;
            }
            else {
                posvalue <- leaf
                leaf <- leaf + 1
            }
            edge[which(edge == -j)] <- posvalue;
        }
    }

    phy <- list(edge = edge)
    phy$tip.label <- paste("t", sample(length(leaves)), sep = "")
    phy$edge.length <- edge.length;
    phy$Nnode <- length(leaves)
    class(phy) <- "phylo"
    return(phy)
}

.mecca.generatetipsummaries <-
#generateTipSummaries <-
function(phy, model, cladeAges, treelist, prop, hotBranches = NULL, richness=richness) {

	Smean	<- numeric(length(treelist))
	Svar <- numeric(length(treelist))

    root.age <- max(branching.times(phy));
	cladeAncStates <- .mecca.nodesim(phy, model, prop, hotBranches);

	for(i in 1:length(treelist)) {

		if(!is.numeric(treelist[[i]])) {
			if(!model == "twoRate") dat <- .mecca.tipsim(treelist[[i]], model, prop, cladeAncStates[i], root.age);
			if(model == "twoRate") {
				tipNum <- which(phy$tip.label == names(richness)[i]);
				if (tipNum %in% hotBranches) { Branch = "hot" } else { Branch <- "not" }
				dat <- .mecca.tipsim(phy = treelist[[i]], model = model, prop = prop, cladeAncStates[i], BranchState = Branch);
			}

 			Svar[i]<-var(dat); # here is the variance of the tip data for a character evolved under the prior for sigma. we want to know if it provides a good approximation to our data, so...
			Smean[i]<-mean(dat);

		} else {
				if(model == "twoRate") {
                    if (i %in% hotBranches) {
                        Branch = "hot"
                    } else {
                        Branch = "not"
                    }

                    dat <- .mecca.singletiptrait(model, prop, cladeAges[i], cladeAncStates[i], root.age, BranchState = Branch);
				} else {
                    dat <- .mecca.singletiptrait(model, prop, cladeAges[i], cladeAncStates[i], root.age);
                }

                Smean[i] <- dat;
                Svar[i] <- 0;
		}

	}

return(list("Smean" = Smean, "Svar" = Svar));
	}

.mecca.gethotbranches <-
#get.hot.branches <-
function(phy, tip1, tip2) {


	if(is.null(tip2)) {

		if(is.numeric(tip1)) {

			return(tip1);
		} else {

			return(which(phy$tip.label==tip1));
		}

	} else {

		mrca <- .mecca.getmrca(phy, tip1, tip2);
		return(.mecca.allbranches(phy, mrca))
	}
}

.mecca.getdivMLE <-
#getDivMLE <-
function(phy, richness) {
	#require(laser)
# a function to compute the maximum likelihood estimators of birth and death from a richness tree. uses the function fitNDR_1rate from Laser (Rabosky 2006)


# Arguments
	# phy <- a phylogenetic tree
	# richness <- a vector of species richness values for the tips in phy

# Values

	# lik <- the log Likelihood of the MLE
	# div <- the  MLE of diversification rate
	# lambda <- the  MLE of birth rate
	# mu <- the  MLE of extinction rate
	phy2 <- .laser.gettipdata(richness, phy);
	ext <- seq(0, 0.99, 0.005); # gives extinction fractions from 0 to 99%

	lik<-numeric(length(ext));
	div<-numeric(length(ext));
	lam<-numeric(length(ext));

	for(i in 1:length(ext)) {

		.laser.fitNDR_1rate(phy2, eps = ext[i], rbounds = c(0.0001, .5), combined = TRUE) -> fit;
		lik[i] <- fit$LH;
		div[i] <- fit$r;
		lam[i] <- fit$lambda;

		}
names(lik) <- ext;
names(div) <- ext;
names(lam) <- ext;


return(list("lik" = max(lik), "div" = div[which(lik==max(lik))],"lamda" = lam[which(lik==max(lik))], "mu" = ext[which(lik==max(lik))] * lam[which(lik==max(lik))]
));
}

.mecca.getproposal <-
#getProposal <-
function(currentState, psi, min, max) {

# This function returns a proposal using a sliding window
# Arguments;

	# current state <- the current state of the parameter;
	# psi <- a tuning parameter - here, some proportion of the standard deviation from the calibration step;
	# min <- the lower bound for the prior distribution;
	# max <- the upper bound for the prior distribution;

# Value;
	# prop <- the proposal;

	prop <- currentState + ((runif(1) - 0.5) * psi);

	if (prop < min) { prop <- min + (min - prop) } # if lower than lower bound, bounce back;
	if (prop > max) { prop <- max - (prop - max) } # if higher than upper bound, bounce back;

	return(prop);

	}

.mecca.getmrca <-
#getmrca <-
function(phy, tip1, tip2) {

	if(is.null(tip2)) {

		bb<-which(phy$tip.label==tip1);
		mrca<-.mecca.parentnode(phy, bb);
	} else {
		nn <- phy$Nnode;
		nt <- length(phy$tip.label);
		minLeaves <- nt;
		mrca <- NULL;
		for(i in 1:nn) {
			leaves <- node.leaves(phy, i + nt);
			if(tip1 %in% leaves & tip2 %in% leaves) {
				ll <- length(leaves);
				if(ll < minLeaves) { mrca <- i + nt; minLeaves <- ll }
			}
		}
	}
	return(mrca);
	}

.mecca.logLP <-
#logLP <-
function(phy, b, d) {
# function that computes the phylogenetic likelihood for a backbone tree following Rabosky et al 2007

#Arguments
	# phy <- a phylogenetic tree
	# b <- the proposed birth rate
	# d <- the proposed death rate

	r <- b - d; # computes the diversification rate
	a <- d / b; # computes the extinction fraction
	N <- phy$Nnode - 1; # computes the number of internal branches

	Xi <- cbind(phy$edge, phy$edge.length);
	Xi <- Xi[-which(Xi[ ,2] <= length(phy$tip.label)), ];
	Xi <- cbind(Xi, numeric(length(Xi[ ,3])));
	bt <- branching.times(phy);
	for (i in 1: length(Xi[ ,4])) {
		Xi[i, 4] <- bt[which(names(bt)==Xi[i, 1])];
		}



	part3 <- numeric(length(Xi[ ,4]));

	for (p in 1:length(part3)) {

		part3[p]<- log(1 - a * exp(-r*Xi[p, 4]));

		}

	logLp <- (N * log(r) - r * sum(Xi[ ,3]) - sum(part3));

	return(logLp);

	}

.mecca.logLT <-
#logLT <-
function(phy, richness, b, d) {
# function that computes the taxonomic likelihood for terminals with species richness in a backbone tree following Rabosky et al 2007

#Arguments
	# phy <- a phylogenetic tree
	# richness <- a vector of species richness values
	# b <- the proposed birth rate
	# d <- the proposed death rate


	r <- b-d;
	a <- d/b;
	t <- .mecca.cladeage(phy);
	n <- richness[match(phy$tip.label, names(richness))];
	beta <- (exp(r * t) - 1) / (exp(r * t) - a);

	logLT <- sum(log(1-beta)) + sum((n-1)*log(beta));
	return(logLT);

	}

.mecca.nodesim <-
#nodeSim <-
function(phy, model, prop, hotBranches = NULL) {


	if(model == "BM") {

		cladeAncStates <- .mecca.cladestates(.mecca.fastbm(phy, prop$root, prop$sigmasq, mu = 0)$nodeStates, phy);

	}


	if(model == "Trend") {

		cladeAncStates <- .mecca.cladestates(.mecca.fastbm(phy, prop$root, prop$sigma, prop$mu)$nodeStates, phy);


	}

	if(model == "twoRate") {

		cladeAncStates <- .mecca.cladestates(.mecca.nodestatestworate(phy, prop, hotBranches), phy);


	}

	return(cladeAncStates)

}

.mecca.paramsample <-
#paramSample <-
function(model, prop, Lk, calSim) {



	if(model == "BM") x <- c(prop$birth, prop$death, Lk, log(prop$sigmasq),prop$root,   calSim$Svar, calSim$Smean);


	if(model == "Trend") x <- c(prop$birth, prop$death, Lk,log(prop$sigmasq), prop$root, prop$mu,  calSim$Svar, calSim$Smean);

	if(model == "twoRate") x <- c(prop$birth, prop$death, Lk, log(prop$sigmasq1),log(prop$sigmasq2), prop$root, calSim$Svar, calSim$Smean);

	return(x);

}

.mecca.parentnode <-
#parent.node <-
function(phy, tip){

## a function to return the parent node of a given tip

# Arguments
	# phy <- a phylogenetic tree
	# tip <- a tip name

# Value
	# a node number

	return(phy$edge[which(phy$edge[ ,2] == tip), 1]);

	}

.mecca.priordrawbirth <-
#priorDrawBirth <-
#UNUSED
function(min, max) {
## function to generate draws from the birth rate prior for calibration
	return(runif(1, min, max));
	}

.mecca.priordrawdeath <-
#priorDrawDeath <-
#UNUSED
function(min, max ) {
## function to generate draws from the death rate prior for calibration

	return(runif(1, min, max));
	}

.mecca.priordrawmean <-
#priorDrawMean <-
function(rootA, rootB, rootPriorType = "uniform") {
## function to generate draws from the root state prior for calibration - note that this function requires and extra argument be specified in the calibration call - whether the prior is normal or uniform. Normal priors might be used if fossil info is available. Other prior types could easily be added though.


	if (rootPriorType == "uniform") {
		return(runif(1, rootA, rootB));
		}

	if(rootPriorType == "normal") {
		return(rnorm(1, mean = rootA, sd = rootB));
		}

	}

.mecca.priordrawsig <-
#priorDrawSig <-
function(min , max, sigmaPriorType, SigmaBounds) {
## function to generate draws from the Brownian rate prior for calibration
	if(sigmaPriorType == "uniform") return(exp(runif(1, min, max)));
	if(sigmaPriorType == "normal") {
		if(is.null(SigmaBounds)) {
			return(exp(rnorm(1, mean = min, sd = max)))
		} else {
			while(1) {
				s <- rnorm(1, mean = min, sd = max)
				if(s > SigmaBounds[1] & s < SigmaBounds[2]) (break)
			}
			return(exp(s))
			}
		}
	}

mecca <-
#runMECCA <-
function(phy, richness, cladeMean, cladeVariance, model = c("BM", "Trend", "twoRate"), prior.list = list(priorSigma = c(-4.961845, 4.247066), priorMean = c(-10, 10)), start = start, Ngens = 10000, printFreq = 100, sigmaPriorType = "uniform", rootPriorType = "uniform", SigmaBounds = c(-4.961845, 4.247066),hotclade = NULL, divPropWidth = 0.1, scale = 1, divSampleFreq = 0, BoxCox = TRUE, outputName ="mecca") {  ##

    ## require pls package
    ## require(pls)

	## To run MECCA, the following arguments are required
	# phy - an incompletely resolved phylogenetic tree
	# richness - a named vector of species richness values for the tips in the phylogenetic tree
	# cladeMean - a named vector of mean trait values for the tips in the phylogenetic tree
	# cladeVariance - a named vector of trait variances for the tips in the phylogenetic tree
	# cal - the output of calibrateMECCA
	# Ngens - the number of generations that MECCA will be run for


	## house keeping block - check everything is in the correct order, get rid of node labels and generate a vector of clade age values ##
	name.check(phy, richness) -> nc;
	if(!nc == "OK") { stop("names in tree and data do not match") }
	if(!sum(is.na(richness) == 0)) { stop("richness contains missing values") }
	if(!sum(is.na(cladeMean) == 0)) { stop("cladeMeans contains missing values") }
	if(!sum(is.na(cladeVariance) == 0)) { stop("cladeVariance contains missing values") }
	phy$node.label <- NULL;
	hotBranches <- NULL;

	match(phy$tip.label, names(richness)) -> m;
	richness[m] -> richness;
	cladeMean[m] -> cladeMean;
	cladeVariance[m] -> cladeVariance;
	cladeVariance <- cladeVariance[-which(cladeVariance==0)];

	if(sigmaPriorType == "uniform") SigmaBounds <- prior.list$priorSigma;

	cladeAges <- .mecca.cladeage(phy);

	dcrit <- start$dcrit;
	obs <- start$obsTraits;
	plsobs <- start$plsObserved
	pls <- start$plsLoadings;
	tuning <- as.data.frame(t(start$tuning));
	ncomp <- ncol(pls);

	if(BoxCox == TRUE) {
		stdz <- start$stdz;
		lambda <- start$lambda;
		GM <- start$GM;
		boxcox <- start$BoxCox;
	}



	## end of name checking and data ordering ##

	if(model =="BM") k <- 2; if(model == "Trend"|model == "twoRate") k <- 3;

	if(model=="twoRate") {

		if(is.null(hotclade)) { stop("if using a two rate model you need to specify the hot clade") }

		hotBranches <- .mecca.gethotbranches(phy, hotclade[1], hotclade[2]);

	}


	if(model == "Trend") {
		if (is.null(prior.list$priorMu)) {prior.list$priorMu = c(-0.5,0.5); warning("No mu prior specified for Trend model. Using default prior of -0.5:0.5") }
		}


	## we generate a matrix to hold our sampled parameters and their associated summary statistics ##

	## will hold the pls transformed summary data

	distSim <- file(paste(outputName, "distSimFile.txt", sep = "_"), "w");
	cat(colnames(tuning), paste("pls", 1:length(plsobs), sep = ""), file = distSim, sep = "\t","\n");

	## will hold the untransformed summary data
	bmSim <- file(paste(outputName, "bmSimFile.txt", sep = "_"), "w");
	cat(colnames(tuning), names(obs), file = bmSim, sep = "\t","\n");

	bdSim <- file(paste(outputName, "bdSimFile.txt", sep = "_"), "w");
	cat(paste("Lambda", "Mu",  "lkl"), file = bdSim, sep = "\t","\n");



	## ABCtoolkit requires a matrix of observed data in the same order as the posterior sample. MECCA will provide this for you as yourMECCAoutputName$observed

	## pls transformed observed
	names(plsobs) <- c(paste("pls", 1:length(plsobs), sep = ""));

	write.table(t(plsobs), paste(outputName, "distObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F);

	write.table(t(obs), paste(outputName, "ObsFile.txt", sep = "_"), row.names = F, sep = "\t", quote = F);

	# TRACK PROPOSALS ETC

	#proposalMat <- matrix(data = NA, nrow = Ngens, ncol = 3 + length(obs));
	#colnames(proposalMat) <- c("rate", "root", "distance", names(obs));

	## end of matrix creation ##



	nAcceptBD = 0;
	nAcceptSigma = 0;
	nAcceptRoot = 0;

	## we initialize MECCA with the mean of the within-tolerance parameters sampled during calibration

	thetaB <- start$startingBirth;
	thetaD <- start$startingDeath;

	trait.params <- tuning[1,];
	trait.widths <- tuning[2,];

	Lk <- .mecca.logLP(phy, thetaB, thetaD) + .mecca.logLT(phy, richness, thetaB, thetaD);

	## and now we can run MECCA ##


	for(ngen in 1:Ngens) {

 	while(1) {
 		## use a sliding window proposal range with width derived from the calibration step ##
		thetaBprop<-.mecca.getproposal(thetaB, divPropWidth, min = 0, max = Inf);
		thetaDprop<-.mecca.getproposal(thetaD, divPropWidth, min = 0, max = Inf);

		if(thetaBprop > thetaDprop) ## we will only proceed with the proposed values if they are greater than 0 and b > d

		break;
	}

	propLk <- .mecca.logLP(phy, thetaBprop, thetaDprop) + .mecca.logLT(phy, richness, thetaBprop, thetaDprop); # compute the joint taxonomic/phylogenetic likelihood of the tree/richness data based on propB and propD


	if(propLk== - Inf) {
		LKratio <- 0 ;
	} else {
		LKratio <- exp(propLk - Lk);
	}

	#### At this point we need to either accept or reject the proposed B and D values - we accept if LKratio >/= 1 and accept with probability p if not. If we accept, we simulate trees under these values and simulate trait evolution on them for abc - if we don't accept, we simulate data and trees using the current states and store the likelihoods and summaries and distances. We might still reject the proposals for sigmasq and root state if traits aren't close though
	if(LKratio >= 1 || runif(1) < LKratio) {
		thetaB <- thetaBprop; thetaD <- thetaDprop; Lk <- propLk;
		nAcceptBD <- nAcceptBD + 1;
	}

	cat(thetaB, thetaD, Lk, file = bdSim, sep = "\t","\n")

	##### now that we have updated our diversification parameters or rejected the proposed move we can move on to simulating trees and traits ####

	prop.params <- .mecca.proposal(model, trait.params, trait.widths, prior.list, SigmaBounds, scale, ngen);

	sim.props <- .mecca.simparamlist(thetaB, thetaD, prop.params);

	while(1) {
		if(ngen  == 1 || divSampleFreq == 0 || ngen  %% divSampleFreq == 0) {
			Sim <- sim.mecca(phy, richness, cladeAges, model, sim.props, makeNewTipTrees = TRUE, mytiptrees = NULL, hotBranches = hotBranches);
			tipTrees <- Sim$trees;
		} else {
			Sim <- .mecca.generatetipsummaries(phy, model, cladeAges, tipTrees, prop = sim.props, hotBranches = hotBranches, richness=richness); ## needs editing!!!!!!
		}

		if(sum(Sim$Svar==0)>0) Sim$Svar <- Sim$Svar[-which(Sim$Svar==0)];
		 #### compute the PLS transformed summaries

		 sims <- c(Sim$Svar, Sim$Smean);
		 names(sims) <- names(obs)
		 if(BoxCox == TRUE) {
			bxsims <- .mecca.boxcox(as.numeric(sims), stdz, lambda, GM, boxcox)
			}
		 if(sum(is.na(bxsims))==0) break;
 	 }

 	 plsSim <- .mecca.extractpls(as.numeric(bxsims), pls, ncomp);
 	 dSim <- abs(dist(rbind(plsobs, plsSim))[1]);
 	 #proposalMat[ngen, ] <- c(thetaSprop, thetaMprop, dSim, sims);

   ####### start block for odd generations

 	if(dSim <= dcrit) {

 		if(ngen %% 2 == 1) {

 			acc <- .mecca.acceptance(trait.params, prop.params, sigmaPriorType, prior.list$priorSigma);

 		} else acc <- "accept";

 	} else if(dSim > dcrit) {
 		acc <- "reject";

 	}

 	 if(acc == "accept") {
		trait.params <- prop.params; ## accept the proposals
		if(ngen %%2 == 1) nAcceptSigma = nAcceptSigma + 1;
		if(ngen %% 2 == 0) nAcceptRoot = nAcceptRoot + 1;
		cat(as.numeric(trait.params), plsSim, file = distSim, sep = "\t", "\n");
		cat(as.numeric(trait.params), sims, file = bmSim, sep = "\t", "\n");

	}

	if(acc == "reject") {

		sim.curr <- .mecca.simparamlist(thetaB, thetaD, trait.params);

	while(1) {
		Sim <- .mecca.generatetipsummaries(phy,model, cladeAges, tipTrees, sim.curr, hotBranches = hotBranches, richness=richness)
		if(sum(Sim$Svar==0)>0) Sim$Svar <- Sim$Svar[-which(Sim$Svar==0)];
			sims <- c(Sim$Svar, Sim$Smean);

			if(BoxCox == TRUE) {
 	 			bxsims <- .mecca.boxcox(sims, stdz, lambda, GM, boxcox)
 	 		}
 	 		if(sum(is.na(bxsims))==0) break;
 	 		}
 	 		plsSim <- .mecca.extractpls(bxsims, pls, ncomp);
				cat(as.numeric(trait.params), plsSim, file = distSim, sep = "\t", "\n");
				cat(as.numeric(trait.params), sims, file = bmSim, sep = "\t", "\n");
			}



	if(ngen %% printFreq == 0) {
			cat("Generation", ngen, "/", round(nAcceptBD/ngen, 2), round(nAcceptSigma/(ngen/2), 2), round(nAcceptRoot/(ngen/2), 2), "\n\n");
			}

	} ## END MCMC


	close(bmSim); close(bdSim); close(distSim);



}

.mecca.simparamlist <-
#sim.param.list <-
function(birth, death, trait) {

	rates <- grep("sigma", colnames(trait));
	trait[,rates] <- exp(trait[,rates]);

	trait <- as.list(trait);

	trait$birth <- birth;
	trait$death <- death;

	return(trait);


}

.mecca.singletiptrait <-
#singleTipTrait <-
function(model, prop, clade.age, clade.anc.state, root.age, BranchState="not") {
#JME: default to BranchState is 'not' (was previously Branch, with no global binding)


	if(model == "BM") {
		dat <- rnorm(1, mean = clade.anc.state, sd = sqrt(prop$sigmasq * clade.age));
	}

	if(model == "Trend") {
		dat <- clade.anc.state + rnorm(1, mean = clade.age*prop$mu, sd =sqrt(prop$sigma * clade.age))
	}

	if(model == "twoRate") {
		if(BranchState == "not") {
			dat <- rnorm(1, mean = clade.anc.state, sd = sqrt(prop$sigmasq1 * clade.age));
		} else if(BranchState == "hot") {
			dat <- rnorm(1, mean = clade.anc.state, sd = sqrt(prop$sigmasq2 * clade.age));
	}
}
	return(dat);

}
