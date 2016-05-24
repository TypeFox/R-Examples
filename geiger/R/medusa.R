medusaVersion = 1.41;

medusa <- function (phy, richness = NULL, criterion = c("aicc", "aic"), partitions = NA, threshold = NA, model = c("mixed", "bd", "yule"),
	cut = c("both", "stem", "node"), stepBack = TRUE, init = c(r = 0.05, epsilon = 0.5), ncores = NULL, verbose = FALSE, ...) {

	## CHECK ARGUMENTS
	initialE <- init[["epsilon"]];
	initialR <- init[["r"]];
	sp = c(initialR, initialE);

	## Determine whether multiple cores can be used (i.e. non-Windows, non-GUI).
	## This is apparently breaking the build. Fix!
	fx <- .get.parallel(ncores);

	criterion <- match.arg(criterion, choices = c("aicc", "aic"), several.ok = FALSE);
	model <- match.arg(model, choices = c("mixed", "bd", "yule"), several.ok = FALSE);
	shiftCut <- match.arg(cut, choices = c("both", "stem", "node"), several.ok = FALSE);

	## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
	if (!any(c("phylo", "multiPhylo") %in% class(phy))) {
		stop("'phy' must either be a phylo or multiPhylo object");
	}
	richness <- .check.richness(phy = phy, richness = richness);
	phyData <- .treedata.medusa(phy = phy, richness = richness, warnings = FALSE); ## modified prune.tree.merge.data for multiple trees (jme)

	## Determine correct AICc threshold from tree size (based on simulations), or use user-provided threshold
	threshold_N <- .threshold.medusa(threshold = threshold, phy = phyData$phy, criterion = criterion);

	## Limit on number of piecewise models fitted; based on tree size, aicc correction factor,
	## and flavour of model fitted (i.e. # parameters estimated; birth-death or pure-birth)
	model.limit <- .get.max.model.limit(richness = richness, partitions = partitions, model = model, verbose = verbose);

	# internal function for running a single tree and richness data.frame (jme)
	medusa_runner <- function (phy, richness) {

		## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
		pend.nodes <- seq_len(length(phy$tip.label)); # Calculate pendant splits just once, keep track through various models
		int.nodes <- unique(phy$edge[,1])[-1]; # Omit root node
		root.node <- length(phy$tip.label) + 1;
		all.nodes <- c(pend.nodes, root.node, int.nodes);

		## Store pertinent information: branch times, richness, descendants
		#cat("Preparing data for analysis... ");
		obj <- .make.cache.medusa(phy = phy, richness = richness, fx = fx, shiftCut = shiftCut);
		#cat("done.\n");
		desc <- list(stem = obj$desc.stem, node = obj$desc.node);

		# z.orig is no longer used, but helpful when debugging
		z <- z.orig <- obj$z;

		## Fit single model to entire tree
		fit <- .fit.base.medusa(z = z, sp = sp, model = model, criterion = criterion);

	# If only one model is desired (i.e. base model), don't bother with all of the precalculations.
		prefit <- NULL;
		if (model.limit != 1) {

		## Needed downstream; do not recalculate
		## Gives the number of tips associated with an internal node; determines whether a node is 'virgin' or not
			num.tips <- fx(all.nodes, function (x) length(obj$tips[[x]]));

		## Pre-fit pendant edges so these values need not be re(re(re))calculated; amounts to ~25% of all calculations
		## Will show particular performance gain for edges with many fossil observations
			#cat("Optimizing parameters for pendant edges... ");
			# Will always be shiftCut="stem"; if mixed model, keep only best fit and throw out other in .prefit.medusa
			tips <- fx(pend.nodes, .prefit.tip.medusa, z = z, sp = sp, model = model, criterion = criterion);
			#cat("done.\n");

		## Pre-fit virgin internal nodes; should deliver performance gain for early models, and especially for large trees
		## Remain useful until a split is accepted within the clade
			virgin.stem <- list();
			virgin.node <- list();

			if (shiftCut == "stem" || shiftCut == "both") {
				virgin.stem <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
					model = model, shiftCut = "stem", criterion = criterion);
			}
			if (shiftCut == "node" || shiftCut == "both") {
				virgin.node <- fx(int.nodes, .prefit.medusa, z = z, desc = desc, sp = sp,
					model = model, shiftCut = "node", criterion = criterion);
			}
			virgin.nodes <- list(stem = virgin.stem, node = virgin.node);
			#cat("done.\n\n");
			prefit <- list(tips = tips, virgin.nodes = virgin.nodes);
		}

		## this will stop in one of two ways:
		## 1) model doesn't improve AIC by 'threshold' (determined by MEDUSA, or set by user)
		## 2) hit the maximum allowed number of models (almost certainly only by user specified value of 'partitions')
		runStepwiseAIC <- function (fit) {
			optModel <- fit;
			i <- 1;
			done <- FALSE;

			cat("Step 1: lnLik=", round(optModel$lnLik, digits=7), "; ", criterion, "=",
				round(as.numeric(optModel[criterion]), digits=7), "; model=", optModel$model[1], "\n", sep="");

			while (!done & i < model.limit) {
				i <- i + 1;
				node.list <- all.nodes[-fit$split.at];
				res <- fx(node.list, .update.fit.medusa, z = z, desc = desc, fit = fit, prefit = prefit, num.tips = num.tips,
					root.node = root.node, model = model, criterion = criterion, shiftCut = shiftCut);

			# Select model with best score according to the specific criterion employed (default aicc)
				fit <- res[[which.min(unlist(lapply(res, "[[", criterion)))]];
				z.best <- .split.z.at.node.medusa(node = tail(fit$split.at, 1), z = z, desc = desc, shiftCut = tail(fit$cut.at, 1))$z;
				step <- rbind(optModel$step, c("add", tail(fit$split.at, 1)));

			# Consider parameter removal
				if (stepBack) {
					backFit <- .back.step.medusa(currentModel = fit, z = z.best, step = step, model = model, criterion = criterion);
					fit <- backFit$fit;
					z.best <- backFit$z;
					step <- backFit$step;
				}
				fit$step <- step;
			# Compare last accepted model to current best model
				if (as.numeric(optModel[criterion]) - as.numeric(fit[criterion]) < threshold_N) {
					cat("\nNo significant increase in ", criterion, " score. Disregarding subsequent piecewise models.\n\n", sep="");
					done <- TRUE;
					break;
				} else {

					if (stepBack && !is.null(backFit$remove)) {.print.removed.shifts(remove=backFit$remove);}
					optModel <- fit;
					z <- z.best;

					cat("Step ", i, ": lnLik=", round(fit$lnLik, digits=7), "; ", criterion, "=",
						round(as.numeric(optModel[criterion]), digits=7), "; shift at node ", tail(fit$split.at,1), "; model=",
						tail(fit$model,1), "; cut=", tail(fit$cut.at,1), "; # shifts=", length(fit$split.at) - 1, "\n", sep="");
				}
			}
			return(list(optModel=optModel, z=z));
		}

		res <- runStepwiseAIC(fit = fit);

		optModel <- res$optModel;
		z <- res$z;
		modelSummary <- .summary.modelfit.medusa(optModel);

		summarizeZ <- function (z) {
			opt.model <- data.frame(cut = optModel$cut.at, split = optModel$split.at, optModel$par, lnLp = optModel$lnLik.part, stringsAsFactors = FALSE);
			break.pts <- as.numeric(optModel$split.at);
			cut.at <- optModel$cut.at;
			rownames(z) <- phy$hash[z[, "dec"]]; # identify identical edges among trees

			# temp hack. old format used NA for split.at for root
			root <- min(z[,"anc"]);
			break.pts[which(break.pts == root)] <- NA;

			# collect shift times
			internal <- is.na(z[, "n.t"]);
			res <- numeric(nrow(z));
			res[] <- NA;
			check <- !is.na(break.pts);
			cut.at <- cut.at[check];
			break.pts <- break.pts[check];
			if (length(break.pts)) {
				for (i in 1:length(break.pts)) {
					idx <- which(z[, "dec"] == break.pts[i]);
					if (cut.at[i] == "node" && internal[idx]) {
						res[idx] <- z[idx, "t.1"];
					} else if (cut.at[i] == "stem" || !internal[idx]) { # pendant node
						res[idx] <- z[idx, "t.0"];
					} else {
						stop("'cut' should either be stem or node");
					}
				}
			}
			anc <- match(z[, "anc"], z[, "dec"]);
			z <- cbind(z, shift = match(z[, "dec"], optModel$split), t.shift = res);
			z <- cbind(z, r = optModel$par[z[, "partition"],"r"], epsilon = optModel$par[z[, "partition"],"epsilon"]);
			z <- cbind(z, ancestral.r = z[anc, "r"], ancestral.epsilon = z[anc, "epsilon"]);
			return(z);
		}

		zSummary <- summarizeZ(z);

		control <- list(threshold = structure(threshold_N, names = criterion), partitions = partitions);
		results <- list(control = control, cache = list(desc = desc, phy = phy), model = optModel,
			summary = modelSummary, zSummary = zSummary, medusaVersion = medusaVersion);
		class(results) <- c("medusa", class(results));
		return(results);
	}

	## I don't like this format. Why store N (possibly 10,000) copies of richness?
	# seems like hash is ONLY useful with mulitple trees.
	if ("multiPhylo" %in% class(phy)) { ## deal with multiple trees
		#res <- lapply(phyData, function (x) medusa_runner(x$phy, x$richness));
		res <- lapply(phyData$phy, function (x) medusa_runner(phy = x, richness = phyData$richness));
		names(res) <- names(phy);
		res$richness <- phyData$richness;
		class(res) <- c("multimedusa", class(res));

	} else {
		res <- medusa_runner(phyData$phy, phyData$richness);
		res$model$prof.par <- .get.profile.likelihoods(res); # add profile likelihoods on parameter values
		res$summary <- cbind(res$summary, res$model$prof.par);
		res$richness <- phyData$richness;

		print(res$summary);
	}

	invisible(res);
}

##############################################################################
##############################################################################


# make sure things are in the correct order and of correct format
.check.richness <- function (phy, richness = NULL) {
	if (is.null(richness)) {
		if ("multiPhylo" %in% class(phy)) {
			phy <- phy[[1]];
		}
		richness <- data.frame(taxon = phy$tip.label, n.taxa = 1);
	} else {
		richness = data.frame(richness, stringsAsFactors = FALSE);
		if (length(richness[1, ]) == 2) {
			if (colnames(richness)[1] != "taxon" || colnames(richness)[2] != "n.taxa") {
				if (class(richness[, 1]) == "factor" & class(richness[, 2]) == "integer") {
					colnames(richness) = c("taxon", "n.taxa");
				} else if (class(richness[, 1]) == "integer" & class(richness[, 2]) == "factor") {
					colnames(richness) = c("n.taxa", "taxon");
				} else {
					stop("'richness' data appear incorrectly formated: see medusa()");
				}
			}
		}
	}
	return(richness);
}

## Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
##   Perhaps relax on these column names, may cause too many problems
## May also include 'exemplar' column; in that case, rename relevant tip.label before pruning. this is currently broke.
#prune.tree.merge.data
.treedata.medusa <- function (phy, richness = NULL, ...) {

	## MODIFIED -- jme
	if ("multiPhylo" %in% class(phy)) { ## deal with multiple trees
		res = lapply(phy, .treedata.medusa, richness, ...);
		tips = c();
		for (i in 1:length(res)) {
			tips = union(tips, res[[i]]$phy$tip.label);
		}
		nn = sapply(res, function (x) Ntip(x$phy));
		if (!all(nn == length(tips))) {
			stop("'phy' appears to have trees with incompletely overlapping sets of tips");
		}
		trees = lapply(res, "[[", "phy");
		class(trees) = "multiPhylo";
		trees = hashes.phylo(phy = trees, tips = tips) ## FROM PHYLO package [returns trees with $hash object -- a 'key' for each split]
		#for (i in 1:length(res)) res[[i]]$phy = trees[[i]];
		#return(res);
		return(list(phy = trees, richness = richness));
	}
	## END -- jme

	# Rename exemplar taxa with taxon name in richness file. Doesn't work in current implementation
	if (!is.null(richness$exemplar)) {
		# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
		# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar);
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na]);
	}

	# checking for typo; if same size, nothing should be dropped
	check <- FALSE;
	if (length(phy$tip.label) == length(richness[, 1])) {
		check <- TRUE;
	}

	# Prune tree down to lineages with assigned richnesses
	temp <- richness[, "n.taxa"];
	names(temp) <- richness[, "taxon"];
	pruned <- treedata(phy, temp, ...); # geiger function calling ape (namecheck)
	if (check) {
		if (length(phy$tip.label) != length(pruned$phy$tip.label)) {
			stop("'richness' data and tip labels in 'phy' do not appear to match fully");
		}
	}
	phy <- pruned$phy;
	rr <- pruned$data;
	richness <- structure(data.frame(taxon = rownames(rr), n.taxa = rr[, 1], row.names = NULL));
	# Check the tree
	# plotNN(phy); # Node numbers (ape-style) plotted

	return(list(phy = phy, richness = richness));
}

## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion ("model.limit") e.g. when k = n-1 (i.e. when denominator of aicc correction is undefined).
## k <- (3*i-1) # when both birth and death are estimated, where i is the number of piecewise models
## This occurs when i = n/3
## If Yule, max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
## Alternatively use aicc threshold itself as a stopping criterion.
# AICc = AIC + 2*k*(k+1)/(n-k-1);
.get.max.model.limit <- function (richness, partitions, model, verbose) {
	samp.size <- (2 * nrow(richness) - 1);
	if (model == "bd" || model == "mixed") {
		max.model.limit <- as.integer(samp.size/3) - ((!(samp.size%%3)) * 1);
	} else {
		max.model.limit <- as.integer(samp.size/2) - ((!(samp.size%%2)) * 1);
	}

	if (!is.na(partitions)) {
		flag <- "'partitions' should either be NA or a positive integer specifying the maximum number of piecewise models to consider";
		if (!is.numeric(partitions) | partitions <= 0) {
			stop(flag);
		}
		if (partitions > max.model.limit) {
			model.limit <- max.model.limit;
			warning("Supplied 'partitions' is in excess of the maximal number that can be considered");
		} else {
			model.limit <- partitions;
		}
	} else {
		model.limit <- max.model.limit;
	}

	if (verbose) {
		cat("Limiting consideration to a maximum of ", model.limit, " piecewise", sep = "");
		if (model == "bd") {
			cat(" birth-death models");
		} else if (model == "mixed") {
			cat(" mixed models");
		} else {
			cat(" pure-birth (Yule) models");
		}
		cat(" (or until threshold is not satisfied).\n\n");
	}
	return(model.limit);
}

## Fitted curve from random b-d simulations
## Value corresponds to 95th percentile of AICc(split) - AICc(no-split) for no-split simulations
## x-shifted power function
## User may override this with 'threshold' argument
.threshold.medusa <- function (threshold, phy, criterion) {
	if ("multiPhylo" %in% class(phy)) {
		phy <- phy[[1]];
	}
	N <- Ntip(phy);
	a <- -35.9410523803326;
	b <- 6.7372587299747;
	c <- -0.100615083407549;
	Offset <- 27.5166786643334;
	y <- a * (N - b)^c + Offset;
	if (y < 0) {
		y <- 0;
	}

	# check uer-specified value
	if (!is.na(threshold)) {
		if (!is.numeric(threshold)) {
			stop("'threshold' should be numeric, specifying the improvement in AIC model fit that is deemed significant");
		} else if (threshold < 0) {
			stop("'threshold' should be a positive, specifying the improvement in AIC model fit that is deemed significant");
		} else {
			cat("Using user-specified ", criterion, "-threshold of ", threshold,
				" (rather than the default value for a tree of ", N, " tips: ", y, ").\n\n", sep="");
			return(threshold);
		}
	} else {
		cat("Appropriate  ", criterion, "-threshold for a tree of ", N, " tips is: ", y, ".\n\n", sep="");
		return(y);
	}
}

## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## In addition, every node's descendants are also calculated.  The element 'desc' is a list.
## $desc[i] contains the indices within $edge, $t.start, etc., of all descendants of node 'i'
## (in ape node numbering format).
## f: either mclapply or lapply -- see .get.parallel()
.make.cache.medusa <- function (phy, richness, fx, shiftCut) {
	n.tips <- length(phy$tip.label);
	n.int <- nrow(phy$edge) - n.tips;

	## Ape numbers the tips first
	i.int <- seq_len(n.int);
	interior <- phy$edge[, 2] %in% phy$edge[, 1];
	bt <- branching.times(phy);

	# Consider only internal edges first. may be zero if only 2 tips.
	if (n.int > 0) {
		edges.int <- matrix(phy$edge[interior,], nrow = n.int, ncol = 2);
		colnames(edges.int) <- c("anc", "dec");

		t.0 <- bt[match(edges.int[,1], (n.tips + 1):max(edges.int))];
		t.1 <- c(t.0[i.int] - phy$edge.length[interior]);

		z.internal <- cbind(edges.int, t.0, t.1, t.len = t.0 - t.1,
			n.0 = rep(1, n.int), n.t = rep(NA, n.int));
	}

	# Now, pendant edges;
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),];
	colnames(edges.pendant) <- c("anc", "dec");

	t.0 <- bt[match(edges.pendant[, 1], (n.tips + 1):max(edges.pendant))];
	t.1 <- rep(0, n.tips);
	# cannot assume richness ordering necessarily matches that of tip labels
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)];

	z.pendant <- cbind(edges.pendant, t.0, t.1, t.len = t.0 - t.1,
		n.0 = rep(1, n.tips), n.t = ext.richness);

	if (n.int > 0) {
		z <- rbind(z.internal, z.pendant);
	} else { # case with only 2 pendant edges.
		z <- z.pendant;
	}

	z <- cbind(z, partition = rep(1, length(z[, 1]))) # Stores piecewise model structure
	rownames(z) <- NULL;

	# Used for identifying descendant nodes below i.e. tracking breakpoints
	all.edges <- as.matrix(z[, c("anc", "dec")]);
	desc.stem <- list();
	desc.node <- list();

	# only calculate what is needed given value of shiftCut
	if (shiftCut == "both" || shiftCut == "stem") {
		desc.stem <- fx(seq_len(max(all.edges)), .descendants.cutAtStem.idx, all.edges = all.edges);
	}
	if (shiftCut == "both" || shiftCut == "node") {
		if (length(desc.stem) != 0) {
			root <- min(z[,"anc"]);
			desc.node <- fx(desc.stem, .strip.stem); # much faster
			desc.node[root] <- desc.stem[root];
		} else {
			desc.node <- fx(seq_len(max(all.edges)), .descendants.cutAtNode.idx, all.edges = all.edges);
			# for the case of tips (no descendants), add themselves as their descendant
			for (i in 1:length(desc.node)) {
				if (length(desc.node[[i]]) == 0) { # tips
					desc.node[[i]] <- .descendants.cutAtStem.idx(node.list = i, all.edges = all.edges);
				}
			}
		}
	}
	tips=NULL;
	#tips <- .cache.descendants(phy)$tips;
	res <- list(z = z, desc.stem = desc.stem, desc.node = desc.node, tips = tips);
	return(res);
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtStem <- function (node, all.edges) {
	ans <- numeric();
	ans <- node;
	repeat {
		node <- all.edges[all.edges[, 1] %in% node, 2];
		if (length(node) > 0) {
			ans <- c(ans, node);
		} else {
			break;
		}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtStem.idx <- function (node.list, all.edges) {
	which(all.edges[, 1] == node.list | all.edges[, 2] %in% .descendants.cutAtStem(node.list, all.edges));
}

## This generates the indices of all descendants of a node, using ape's edge matrix.
## Deals with row numbers of the edge matrix rather than node numbers of the tree.
.descendants.cutAtNode <- function (node, all.edges) {
	ans <- numeric();
	repeat {
		node <- all.edges[all.edges[, 1] %in% node, 2]
		if (length(node) > 0) {
			ans <- c(ans, node);
		} else {
			break;
		}
	}
	return(unlist(ans));
}

## The function 'descendants' returns the indices of all descendants within the edge matrix.
.descendants.cutAtNode.idx <- function (node.list, all.edges) {
	which(all.edges[, 1] == node.list | all.edges[, 2] %in% .descendants.cutAtNode(node.list, all.edges));
}

# Remove stem node from previously calculated set of stem descendants; about a billion times faster than descendants.cutAtNode
.strip.stem <- function (x) {
	if (length(x) == 1) { # if it is a tip, leave it alone
		return(x);
	}
	y <- unlist(x);
	return(y[-1]);
}


## Only used for base model
#.fit.base.medusa <- function (z, sp, model, fixPar, criterion) {
.fit.base.medusa <- function (z, sp, model, criterion) {
	#fit <- .get.optimal.model.flavour(z = z, sp = sp, model = model, fixPar = fixPar, criterion = criterion);
	fit <- .get.optimal.model.flavour(z = z, sp = sp, model = model, criterion = criterion);
	model.fit <- .calculate.modelfit.medusa(fit = fit, z = z);
	steps <- c("add", as.numeric(min(z[,"anc"])));

	return(list(par=matrix(fit$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=fit$lnLik, lnLik=fit$lnLik,
		split.at=min(z[,"anc"]), aic=round(model.fit$aic, digits=7), aicc=round(model.fit$aicc, digits=7), num.par=model.fit$k,
		cut.at="node", model=fit$model, step=matrix(steps, nrow=1, dimnames=list(NULL,c("step", "node")))));
}

## When model == mixed, fit both and find optimal flavour
#.get.optimal.model.flavour <- function (z, sp, model, fixPar, criterion) {
.get.optimal.model.flavour <- function (z, sp, model, criterion) {
	fit.bd <- NULL;
	fit.yule <- NULL;
	fit <- NULL;

	if (model == "yule" | model == "mixed") {
		fit.yule <- .fit.partition.medusa(z = z, sp = sp, model = "yule");
		fit.yule$model <- "yule";
	}
	if (model == "bd" | model == "mixed") {
		if (is.na(sp[2])) {sp[2] <- 0.5;}
		fit.bd <- .fit.partition.medusa(z = z, sp = sp, model = "bd");
		fit.bd$model <- "bd";
	}
	# if (model != "mixed" && model != "bd" && model != "yule") { # i.e. the constrained models
		# fit <- medusaMLFitPartition(z=z, sp=sp, model=model, fixPar=fixPar);
		# fit$model <- model;
		# return(fit);
	# }
## Figure out which model fits best
	if (is.null(fit.bd)) {
		fit <- fit.yule;
	} else if (is.null(fit.yule)) {
		fit <- fit.bd;
	} else {
## Considering both models
		fit <- .get.best.partial.model(fit1 = fit.yule, fit2 = fit.bd, z = z, criterion = criterion);
	}
	return(fit);
}

## Used for comparing models fit to the same partition
.get.best.partial.model <- function (fit1, fit2, z, criterion) {
	n <- (length(z[,1]) + 1); # the number of nodes involved
## Add '1' to parameters to account for break
	k1 <- 1 + sum(!is.na(fit1$par));
	k2 <- 1 + sum(!is.na(fit2$par));

	if (n - k1 <= 1 || n - k2 <= 1) { # deals with single edges, where AICc correction becomes undefined. use AIC.
		if (.get.AIC(fit1$lnLik, k1) < .get.AIC(fit2$lnLik, k2)) {
			return(fit1);
		} else {
			return(fit2);
		}
	} else {
		if (criterion == "aicc") {
			if (.get.AICc(fit1$lnLik, k1, n) < .get.AICc(fit2$lnLik, k2 ,n)) {
				return(fit1);
			} else {
				return(fit2);
			}
		} else {
			if (.get.AIC(fit1$lnLik, k1) < .get.AIC(fit2$lnLik, k2)) {
				return(fit1);
			} else {
				return(fit2);
			}
		}
	}
}





# replace these at some point with a general function
.get.AIC <- function (lnLik, k) {
	return(-2 * lnLik + 2*k);
}

.get.AICc <- function (lnLik, k, n) {
	return(-2 * lnLik + 2*k*n/(n-k-1));
}



#.prefit.tip.medusa <- function (node, z, sp, model, fixPar, criterion) {
.prefit.tip.medusa <- function (node, z, sp, model, criterion) {
	z.tip <- z[z[,"dec"] == node,,drop=FALSE];

	# tips are always better fit by yule. don't bother with BD.
	if (model == "yule" || model == "mixed") {
		if (z.tip[,"n.t"] == 1) { # single tip. nothing has happened along edge, so r = 0
			return(list(par=c(0, NA), lnLik=0, model="yule"));
		} else { # unresolved clade
			# no t.len information available, only depth
			# MLE birth rate is: log(n.t) / depth
			r <- as.numeric(log(z.tip[,"n.t"]) / z.tip[,"t.len"]);
			lik <- as.numeric(-z.tip[,"t.len"] * r + (z.tip[,"n.t"] - 1) * log(1 - exp(-z.tip[,"t.len"] * r)));
			return(list(par=c(r, NA), lnLik=lik, model="yule"));
		}
	} else { # at the moment, only BD will get through. but eventually constrained models also
		return(.get.optimal.model.flavour(z=z.tip, sp=sp, model=model, criterion=criterion));
	}
}


# prefitting ensure that clades are intact. this can enable shortcuts for certain configurations/models.
.prefit.medusa <- function (node, z, desc, sp, model, shiftCut, criterion) {
	zz <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = shiftCut, extract = TRUE)$z;
	if (all(zz[,"n.t"] == 1, na.rm=T)) {
		fit <- .get.optimal.simple(z = zz, sp = sp, model = model, criterion = criterion);
	} else {
		fit <- .get.optimal.model.flavour(z = zz, sp = sp, model = model, criterion = criterion);
	}
	return(fit);
}

.get.optimal.simple <- function (z, sp, model, criterion) {
	fit.bd <- NULL;
	fit.yule <- NULL;
	fit <- NULL;

	if (model == "yule" | model == "mixed") {
		fit.yule <- .simple.yule.medusa(z = z);
		fit.yule$model <- "yule";
	}
	if (model == "bd" | model == "mixed") {
		if (is.na(sp[2])) {sp[2] <- 0.5;}
		fit.bd <- .fit.partition.medusa(z = z, sp = sp, model = "bd");
		fit.bd$model <- "bd";
	}
## Figure out which model fits best
	if (is.null(fit.bd)) {
		fit <- fit.yule;
	} else if (is.null(fit.yule)) {
		fit <- fit.bd;
	} else {
## Considering both models
		fit <- .get.best.partial.model(fit1 = fit.yule, fit2 = fit.bd, z = z, criterion = criterion);
	}
	return(fit);
}

# solution known when subtree is completely sampled
.simple.yule.medusa <- function (z) {
	n.int <- sum(is.na(z[,"n.t"])); # the number of speciation events
	if (n.int == 0) {
		return(list(par=c(0, NA), lnLik=0));
	}
	sum.t <- sum(z[,"t.len"]);
	r <- (n.int) / sum.t;
	lik <- n.int * log(r) - r * sum.t;
	par <- c(r, NA);
	return(list(par=c(r, NA), lnLik=lik));
}

## Split the edge matrix 'z' by adding a partition rooted at node 'node'.
##   Note: in original MEDUSA parlance, this is cutAtStem=T.
## The list 'desc' is a list of descendants (see make.cache.medusa, above).
## Returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
## This is where 'shiftCut' matters
## Arg 'extract' determines whether only the new clade should be returned. Used in prefitting.
#medusa.split
.split.z.at.node.medusa <- function (node, z, desc, shiftCut, extract=FALSE) {
	descendants <- NULL;
	if (shiftCut == "stem") {
		descendants <- desc$stem;
	} else {
		descendants <- desc$node;
	}
	part <- z[, "partition"];
	base <- min(part[z[, 1] == node | z[, 2] == node]);
	tag <- max(part) + 1;
	i <- descendants[[node]];
	idx <- i[part[i] == base];
	z[idx, "partition"] <- tag;
	#z[which(z["dec"] == node), "partition"] <- tag; # Possible to have several edges to consider

	if (extract) {z <- z[idx,,drop = FALSE];}

	return(list(z = z, affected = c(unique(part[idx]), tag)));
}

## 'fit' contains parameter values from previous model, used to initialize subsequent model.
## Pass in pre-fitted values for pendant edges and virgin nodes (in 'prefit'); DON'T recalculate.
## Need to consider the possibility of birth-death, yule, or mixed models.
## Need to consider where shft is placed (shiftCut). Placement affects not only new clade, but
## also the size of the split clade. Only relevant if shiftCut = "both".
## fit1 model is already logged; only need to record fit2 model, and only non-prefitted nodes
#medusa.ml.update
.update.fit.medusa <- function (node, z, desc, fit, prefit, num.tips, root.node, model, criterion, shiftCut) {
	## various combinations possible
	fit1.stem <- NULL;
	fit1.node <- NULL;
	fit2.stem <- NULL;
	fit2.node <- NULL;
	cut.at <- NULL;

	sp <- NULL;
	aff <- NULL;
	op <- fit$par;
	cut.at <- NULL;

	fit1 <- NULL;
	fit2 <- NULL;

	if (shiftCut == "stem" || shiftCut == "both") {
	## First, diminshed clade
		obj.stem <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "stem");
		z.stem <- obj.stem$z;
		aff <- obj.stem$affected;

	## Ensure that neither partition is empty; can occur with "node" or "both" cutting. If so, kill it.
		if (sum(z.stem[,"partition"] == aff[1]) == 0 || sum(z.stem[,"partition"] == aff[2]) == 0) {
			fit$lnLik <- -Inf;
			fit$aic <- Inf;
			fit$aicc <- Inf;
			return(fit);
		}

	## Everything is cool; proceed.
		sp <- op[aff[1], ]; # Use previously fit parameter values from clade that is currently being split

		dimClade <- z.stem[z.stem[,"partition"] == aff[1],,drop = FALSE];

	## first, consider diminshed clade. may result in a clade that has been cached previously
	## check if diminished clade conforms to a cached clade
		x <- which(dimClade[,"anc"] == min(dimClade[,"anc"]));
		if (length(x) == 1) { # a cut-at-stem scenario
			y <- as.numeric(dimClade[x, "dec"]);
			if (length(unique(dimClade[(dimClade[,"dec"] < root.node),"dec"])) == num.tips[[y]]) {
				if (y < root.node) {
					fit1.stem <- prefit$tips[[y]];
				} else {
					fit1.stem <- prefit$virgin.nodes$stem[[y - root.node]];
				}
			}
		}

		if (is.null(fit1.stem)) { # not cached
			fit1.stem <- .get.optimal.model.flavour(z = dimClade, sp = sp, model = model, criterion = criterion);
		}

	## Second, new clade
		if (node < root.node) { # tip, already calculated
			fit2.stem <- prefit$tips[[node]];
		} else if (length(unique(z.stem[(z.stem[, "partition"] == aff[2] & z.stem[, "dec"] < root.node), "dec"])) == num.tips[[node]]) { # virgin node, already calculated
			fit2.stem <- prefit$virgin.nodes$stem[[node - root.node]];
		} else { # novel shift
			newClade <- z.stem[z.stem[,"partition"] == aff[2],,drop = FALSE]; # only subset if necessary
			fit2.stem <- .get.optimal.model.flavour(z=newClade, sp=sp, model=model, criterion=criterion);
		}
	}

	if ((shiftCut == "node" || shiftCut == "both")  && (node > root.node)) {
		## First, diminshed clade
		obj.node <- .split.z.at.node.medusa(node = node, z = z, desc = desc, shiftCut = "node");
		z.node <- obj.node$z;
		aff <- obj.node$affected;
		sp <- op[aff[1], ]; # Use previously fit parameter values from clade that is currently being split

		dimClade <- z.node[z.node[,"partition"] == aff[1],,drop = FALSE];

		fit1.node <- .get.optimal.model.flavour(z = dimClade, sp = sp, model = model, criterion = criterion);

		if (length(unique(z.node[(z.node[, "partition"] == aff[2] & z.node[, "dec"] < root.node), "dec"])) == num.tips[[node]]) { # virgin node
			fit2.node <- prefit$virgin.nodes$node[[node - root.node]];
		} else { # novel shift
			newClade <- z.node[z.node[,"partition"] == aff[2],,drop = FALSE];
			fit2.node <- .get.optimal.model.flavour(z=newClade, sp=sp, model=model, criterion=criterion);
		}
	}

	## Now, figure out which shift position is optimal
	if (is.null(fit2.node)) {
		fit1 <- fit1.stem;
		fit2 <- fit2.stem;
		cut.at <- "stem";
	} else if (is.null(fit1.stem)) {
		fit1 <- fit1.node;
		fit2 <- fit2.node;
		cut.at <- "node";
	} else {
		## Considering both places for a shift
		stem.lik <- (fit1.stem$lnLik + fit2.stem$lnLik);
		stem.par <- rbind(fit1.stem$par, fit2.stem$par);
		stem.val <- list(lnLik = stem.lik, par = stem.par);
		stem.fit <- .calculate.modelfit.medusa(fit = stem.val, z = z);

		node.lik <- (fit1.node$lnLik + fit2.node$lnLik);
		node.par <- rbind(fit1.node$par, fit2.node$par);
		node.val <- list(lnLik = node.lik, par = node.par);
		node.fit <- .calculate.modelfit.medusa(fit = node.val, z = z);

		if (stem.fit[[criterion]] < node.fit[[criterion]]) {
			fit1 <- fit1.stem;
			fit2 <- fit2.stem;
			cut.at <- "stem";
		} else {
			fit1 <- fit1.node;
			fit2 <- fit2.node;
			cut.at <- "node";
		}
	}
	op[aff[1], ] <- fit1$par # Replace parameters with new values for diminished clade
	fit$model[aff[1]] <- fit1$model; # update altered model

	fit$par <- rbind(op, fit2$par);
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik); # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node);
	fit$lnLik <- sum(fit$lnLik.part);

	model.fit <- .calculate.modelfit.medusa(fit = fit, z = z);

	fit$aic <- model.fit$aic;
	fit$aicc <- model.fit$aicc;
	fit$num.par <- model.fit$k;

	fit$cut.at <- c(fit$cut.at, cut.at);
	fit$model <- c(fit$model, fit2$model);

	return(fit);
}






### NEED TO PUT IN ALL OTHER MODELS HERE!!! ###

## sp = initializing values for r & epsilon
## Default values should never be used (except for first model), as the values from the previous model are passed in
#medusa.ml.fit.partition
#.fit.partition.medusa <- function (partition, z, sp = c(0.1, 0.05), model) {
.fit.partition.medusa <- function (z, sp = c(0.05, 0.5), model) {
	# Construct likelihood function:
	lik <- .lik.partition.medusa(partition = z, model = model);
	foo <- function (x) {-lik(pars = exp(x));} # work with parameters in log-space to preserve precision

	# starting parameters for bd come from previously-fit models (or default for base model)
	if (model == "bd") {
		fit <- optim(fn = foo, par = log(sp), method = "N"); # last argument connotes maximization
		return(list(par = exp(fit$par), lnLik = -fit$value));
	}

## grab values to get an intelligent upper bound on b or r
	node.richness  <- sum(z[,"n.t"], na.rm=TRUE);
	depth <- max(z[,"t.0"]);

	if (model == "yule") {
		#fit <- optimize(f = foo, interval = c(-25, 1));
		#par <- c(exp(fit$minimum), NA);

		n.int <- sum(is.na(z[,"n.t"]));

		if ((n.int == 0) && all(z[,"n.t"] == 1, na.rm = TRUE)) {
			return(list(par=c(0, NA), lnLik=-Inf));
		}

		maxVal <- (log(node.richness) / depth) * 5;
		if (node.richness <= 1) {maxVal <- 1e-5;}

		suppressWarnings(fit <- optimize(f=foo, interval=c(-25, log(maxVal)))); # don't really need to suppress here
		par <- c(exp(fit$minimum), NA);

		#test <- FALSE;
		#if (par[1]/maxVal > 0.95) {
		#	print(z);
		#	cat("node.richness = ", node.richness, "; depth = ", depth, "\n", sep="");
		#	test <- TRUE;
		#}
		while (par[1]/maxVal > 0.95) { # crash against boundary; rarely used
			maxVal <- par[1] * 3;
		#	cat("Hit boundary. Increasing maxVal from ", par[1], " to ", maxVal, "\n", sep="");
			suppressWarnings(fit <- optimize(f=foo, interval=c(log(par[1]/2), log(maxVal)))); # don't really need to suppress here
			par <- c(exp(fit$minimum), NA);
		}
		#if (test) cat("Final par = ", par[1], "\n", sep="");

		return(list(par = par, lnLik = -fit$objective));
	}
	# and all of the constrained models ...

}

## make.lik.medusa.part: generate a likelihood function for a single partition.
#make.lik.medusa.part
## lots of alternative constrained models are missing - add back in
.lik.partition.medusa <- function (partition, model) {
	# Handle internal and pendant edges separately
	is.int <- is.na(partition[, "n.t"]);
	is.pend <- !is.int;

	n.int <- sum(is.int);
	n.pend <- sum(is.pend);

	if (n.int + n.pend != length(partition[, 1])) {
		stop("You messed up, yo.");
	}

	## Internal and pendant calculations differ; split'em up
	int <- partition[is.int, , drop = FALSE];
	pend <- partition[is.pend, , drop = FALSE];

	sum.int.t.len <- sum(int[, "t.len"]); # Simply sum all internal edges
	int.t.0 <- int[, "t.0"];

	# 'n.0' = Foote's 'a', initial diversity; 'n.t' = Foote's 'n', final diversity
	pend.n.0 <- pend[, "n.0"]; # Foote's 'a': initial diversity
	pend.n.t <- pend[, "n.t"]; # Foote's 'n': final diversity
	pend.t.len <- pend[, "t.len"];

	# User may pass in epsilon; don't change it, just estimate r
	f <- function (pars) {
		if (model == "bd") {
			r <- pars[1];
			epsilon <- pars[2];

			if (r < 0 || epsilon <= 0 || epsilon >= 1) {
				return(-Inf);
			}
		} else {
			r <- pars[1];
			epsilon <- 0;
			if (r < 0) {
				return(-Inf);
			}
		}
		l.int <- numeric();
		l.pend <- numeric();

		if (n.int == 0) {
			l.int <- 0;
		} else {
			## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))));
		}

		if (n.pend == 0) {
			l.pend <- 0;
		} else {
## Calculations are from the following:
## Rabosky et al. 2007. Proc. Roy. Soc. 274: 2915-2923.
## Foote et al. 1999. Science. 283: 1310-1314
## Raup. 1985. Paleobiology 11: 42-52 [Foote et al. correct the equation [A18] where a > 1]
## Bailey. 1964. The Elements Of Stochastic Processes, With Applications To The Natural Sciences
## Kendall. 1948. Ann. Math. Stat. 19: 1-15.
##
## A = probability of extinction of one lineage over time 't'
## B = A * (lambda/mu)
##
## When there is a single lineage at time 0 (a = 1), the calculation is
##   log(1 - A) + log(1 - B) + (n - 1)*log(B)
## but this is conditioned on survival by dividing by (1-A)
## (subtracting log(1 - A) on a log scale) which cancels to give:
##   log(1 - B) + (n - 1)*log(B)
##      - for n.t == 1, reduces further to log(1-B)
##
## A = mu*(exp((lambda - mu)*t) - 1)) / (lambda*exp((lambda - mu)*t) - mu)
##  let r = (lambda - mu); ert = exp((lambda - mu)*t)
## A = mu*(ert - 1)/(lambda*ert - mu)
##
## B = A * (lambda/mu)
##   = [mu*(ert - 1)/(lambda*ert - mu)] * (lambda/mu)
##   = (lambda*(ert - 1))/(lambda*ert - mu)
##   = (lambda*(ert - 1))/(lambda(ert - mu/lambda))
##   = (ert - 1) / (ert - epsilon)

## All pendant nodes begin with richness '1'; calculations simple.
# i.pend.n.t.1 <- which(pend.n.t == 1)   # calculations even simpler: log(1-B)
# i.pend.n.t.n1 <- which(pend.n.t != 1)

			ert <- exp(r * pend.t.len);
			B <- (ert - 1)/(ert - epsilon); # Equivalently: B <- (bert - b) / (bert - d)

			l.pend <- sum(log(1 - B) + (pend.n.t - 1) * log(B));
		}
		return(l.int + l.pend);
	}
}

## 'fit' contains '$par' and '$lnlik'
.calculate.modelfit.medusa <- function (fit, z) {
	## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant)
	# Since each edge defines a node (i.e. an 'observation'), need only add root node as final obervation
	n <- (length(z[, 1]) + 1);

	# Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
	## Models where all parameters are estimated (i.e. BD model):
	# 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model

	# Determine number of piecewise models currently involved
	if (length(fit$par) < 3) { # i.e. base model
		num.models <- 1;
	} else {
		num.models <- length(fit$par[, 1]);
	}

	# Updated for more general models: check how many parameter values != NA
	k <- sum(!is.na(fit$par)) + (num.models - 1); # number of estimated parameters + number of breaks

	ll <- list(k = k, lnL = fit$lnLik);
	return(.aic(ll, n));
	#return(.aic(ll$lnL, n, k));
}


## Functions for model removal ##

# ## Consider removing previously-fit rate shifts
# #.back.step.medusa <- function (currentModel, z, step, model, fixPar, criterion) {
.back.step.medusa <- function (currentModel, z, step, model, criterion) {
## As a first step, only consider removing entire shifts. Later deal with individual parameters.
	z.opt <- z;
	bestModel <- currentModel;
	bestScore <- as.numeric(bestModel[criterion]);
	allDeletedShifts <- NULL;
	bestRemoved <- NULL;
	improve <- T;

	while (improve) { # may be possible to remove > 1 previously fit shift
		allDeletedShifts <- c(allDeletedShifts, bestRemoved);
		currentModel <- bestModel;
		z <- z.opt;
		cuts <- bestModel$cut.at;
		nodes <- bestModel$split.at;
		pars <- bestModel$par;
		numModels <- length(bestModel$par)/2;
		improve <- F;

		if (numModels > 2) {
			for (i in 2:(numModels - 1)) { # don't waste time removing last shift
				fitModel <- currentModel;
				obj <- .dissolve.split.medusa(z, cut = cuts[i], node = nodes[i], aff = i);
				aff <- obj$affected;
				z.temp <- obj$z[obj$z[,"partition"] == aff,,drop=FALSE];

		## set par to weighted mean of 2 affected partitions
			## updated to reflect total path length rather than number of tips
			## really only influences weighted parameter starting values
				weights <- c(sum(z[which(z[,"partition"] == aff),"t.len"]), sum(z[which(z[,"partition"] == i),"t.len"]));
				sp <- c(weighted.mean(pars[c(aff,i),1], weights), weighted.mean(pars[c(aff,i),2], weights, na.rm=T));
				#fit <- getOptimalModelFlavour(z=z.temp, sp=sp, model=model, fixPar=fixPar, criterion=criterion);
				fit <- .get.optimal.model.flavour(z=z.temp, sp=sp, model=model, criterion=criterion);

		## Update fit values
				fitModel$par[aff,] <- fit$par;
				fitModel$par <- fitModel$par[-i,];
				fitModel$lnLik.part[aff] <- fit$lnLik;
				fitModel$lnLik.part <- fitModel$lnLik.part[-i];
				fitModel$lnLik <- sum(fitModel$lnLik.part);
				model.fit <- .calculate.modelfit.medusa(fit=fitModel, z=z);
				fitModel$aic <- model.fit$aic;
				fitModel$aicc <- model.fit$aicc;
				fitModel$num.par <- model.fit$k;

				if (fitModel[criterion] < bestScore) {
					#cat("Found a better fit by removing a model (split at node #", fitModel$split.at[i], ")!\n", sep="");
					fitModel$split.at <- fitModel$split.at[-i];
					fitModel$model[aff] <- fit$model;
					fitModel$model <- fitModel$model[-i];
					fitModel$cut.at <- fitModel$cut.at[-i];
					bestModel <- fitModel;
					bestScore <- as.numeric(fitModel[criterion]);
					z.opt <- .updateZ(z=obj$z, deletedPart=i);
					bestRemoved <- nodes[i];
					improve <- T;
				}
			}
			if (improve) {step <- rbind(step, c("remove", bestRemoved));}
		}
	}
	return(list(fit=bestModel, z=z.opt, step=step, remove=bestRemoved));
}

## Remove previously-fit rate shift
.dissolve.split.medusa <- function (z, cut, node, aff) {
## Grab ancestral branch partition membership
	anc <- z[which(z[,"dec"] == node)];
	root <- min(z[,"anc"]);
	tag <- NULL;

	if (cut == "node") {
		tag <- as.numeric(z[which(z[,"dec"] == node),"partition"]);
	} else if (cut == "stem" && anc > root) {
		tag <- as.numeric(z[which(z[,"dec"] == anc),"partition"]);
	} else if (cut == "stem" && anc == root) { # need to take other side of root
		dec <- z[which(z[,"anc"] == root),"dec"];
		tag <- as.numeric(z[which(z[,"dec"] == dec[which(dec != node)]),"partition"]); # ug. li.
	}

	idx <- which(z[,"partition"] == aff);
	z[idx,"partition"] <- tag;

	return(list(z=z, affected=tag));
}

## Only print if model improves AIC score
.print.removed.shifts <- function (remove) {
	for (i in 1:length(remove)) {
		cat("  Removing shift at node #", remove[i], "\n", sep="");
	}
}

## reduce partition IDs to reflect dissolved split
.updateZ <- function (z, deletedPart) {
	idx <- z[,"partition"] > deletedPart;
	z[idx,"partition"] <- z[idx,"partition"] - 1;
	return(z);
}










##########################################################################################
# Printing / plotting / summary
##########################################################################################

## Prints out a table of likelihoods, and parameters.
.summary.modelfit.medusa <- function (optModel) {
	modelSize <- length(optModel$split.at);

	summ <- data.frame(cbind(seq(1:modelSize), optModel$split.at, optModel$cut.at, optModel$model,
		signif(optModel$lnLik.part, digits=7)), signif(optModel$par, digits=6), stringsAsFactors = FALSE);
	colnames(summ) <- c("Model.ID", "Shift.Node", "Cut.At", "Model", "Ln.Lik.part", "r", "epsilon");

	return(summ);
}

## Print out summary matrix with liks, parameter vals, and confidence intervals (if present)
print.medusa <- function (x, ...) {
	n.taxa <- sum(x$zSummary[,"n.t"], na.rm=TRUE);
	n.tips <- length(x$cache$phy$tip.label);
	if (n.taxa == n.tips) {
		cat("\nOptimal MEDUSA model for tree with ", n.tips, " taxa.\n\n", sep="");
	} else {
		cat("\nOptimal MEDUSA model for tree with ", n.tips, " tips representing ", n.taxa, " taxa.\n\n", sep="");
	}
	print(x$summary);
	cat("\n");
	if (!is.null(x$summary$r.low)) {
		cat("95% confidence intervals on parameter values calculated from profile likelihoods\n");
	}
}

print.multimedusa <- function (x) {
	n.trees <- length(x) - 1;
	n.taxa <- sum(x[[1]]$zSummary[,"n.t"], na.rm=TRUE);
	n.tips <- length(x[[1]]$cache$phy$tip.label);
	if (n.taxa == n.tips) {
		cat("\nMEDUSA results for ", n.trees, " trees, each with ", n.taxa, " taxa.\n\n", sep="");
	} else {
		cat("\nMEDUSA results for ", n.trees, " trees, each with ", n.tips, " tips representing ", n.taxa, " taxa.\n\n", sep="");
	}
}

# should add in the ability to plot terminal richnesses (if any != 1)
#plot.medusa <- function (x, partitions = list(cex = 2, bg = "gray", alpha = 0.75, col = "black", lwd = 1), ...) {
#plot.medusa <- function (x, cex = 0.5, time = TRUE, bg = "gray", alpha = 0.75, col = "black", lwd = 1, ...) {
# TODO: (1) come up with custom colour-scheme, (2) plot richnesses if any > 1.
plot.medusa <- function (x, cex = 0.5, time = TRUE, ...) {
	z <- x$zSummary;
	shift <- list(cex = 1, bg = "gray", alpha = 0.75, col = "black", lwd = 1);
	phy <- x$cache$phy;

	mm <- match(phy$edge[,2], z[,"dec"]);
	edge.color <- z[mm, "partition"];

	plot(phy, edge.color=edge.color, cex=0.5, ...);
	if (time) axisPhylo();
	shifts <- z[(idx <- !is.na(z[, "shift"])), "dec"];
	if (length(shifts)) {
		ss <- z[idx, "shift"];
		ww <- character(nrow(phy$edge));
		ww[match(shifts, phy$edge[, 2])] <- ss;
		xx <- numeric(nrow(phy$edge));
		xx[phy$edge[, 2] %in% shifts] <- 1;
		edgelabels.auteur(NULL, frame = "circle", cex = ifelse(xx == 1, shift$cex, 1e-10), pch = ifelse(xx == 1, 21, NA),
		bg = .transparency(shift$bg, shift$alpha), col = shift$col, lwd = shift$lwd);
		edgelabels.auteur(ww, frame = "none", cex = ifelse(xx == 1, shift$cex/3, 1e-10), pch = NA, col = shift$col);
	}
}

## Get confidence intervals on parameter values using profile likelihoods
## This is only used for single-tree analyses. Could be added for multimedusa, but I don't really see the point.
## passed-in parameters parm are MLEs stored in a matrix
.get.profile.likelihoods <- function (res, crit=1.92) {
	parm <- res$model$par;
	models <- res$model$model;
	#z <- res$model$z;
	z <- res$zSummary;

	prof.par <- matrix(nrow=length(parm[,1]), ncol=4);
	colnames(prof.par) <- c("r.low", "r.high", "eps.low", "eps.high")
	inc <- 0.05;

	# get this out of for-loop format
	for (i in 1:length(parm[,1])) {
		#cat("Model",i,"\n")
		model <- models[i]
		sp <- as.numeric(parm[i,]);
		new.part <- z[z[,"partition"] == i,,drop=FALSE];
		lik <- .lik.partition.medusa(partition=new.part, model=model);

		if (model == "yule") {
			par <- sp[1];
			maxLik <- NULL;

			if (par == 0) {
				maxLik <- 0;
			} else {
				maxLik <- lik(par);
				if (maxLik == -Inf) maxLik <- 0; # correct for -Inf at boundary lambda == 0
			}
			threshold <- function (x) lik(x) - maxLik + crit; # find roots on either side of maxLik

	## need intelligent bounds
			if (par != 0) {
				low.bound <- par - par/2;
				up.bound <- par + par/2;
			} else {
				low.bound <- par;
				up.bound <- par + inc/2;
			}
			if (low.bound != 0) {
				while (threshold(low.bound) > 0) {
					low.bound <- low.bound - inc;
				}
			}
			while (threshold(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}

			if (low.bound <= 0) {
				prof.par[i,1] <- 0;
			} else {
				prof.par[i,1] <- uniroot(threshold, lower=low.bound, upper=par)$root;
			}
			if (par == 0) par <- 1e-10; # avoid -Inf at boundary
			prof.par[i,2] <- uniroot(threshold, lower=par, upper=up.bound)$root;

		} else if (model == "bd") {
			par1 <- sp[1]; par2 <- sp[2];
			maxLik <- lik(sp);

	## first, r
			thresholdR <- function (x) lik(c(x, par2)) - maxLik + crit;

			low.bound <- par1 - par1/2;
			up.bound <- par1 + par1/2;

			while (thresholdR(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thresholdR(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			if (low.bound <= 0) low.bound <- 0;

			prof.par[i,1] <- uniroot(thresholdR, lower=low.bound, upper=par1)$root;
			prof.par[i,2] <- uniroot(thresholdR, lower=par1, upper=up.bound)$root;

	## now, epsilon
			thresholdE <- function (x) lik(c(par1, x)) - maxLik + crit;

			low.bound <- par2 - par2/2;
			up.bound <- par2 + par2/2;

			while (thresholdE(low.bound) > 0) {
				low.bound <- low.bound - inc;
			}
			while (thresholdE(up.bound) > 0) {
				up.bound <- up.bound + inc;
			}
			if (low.bound < 0) low.bound <- 0;
			if (up.bound > 1) up.bound <- 1;

			if (low.bound == 0) {
				prof.par[i,3] <- 0;
			} else {
				prof.par[i,3] <- uniroot(thresholdE, lower=0, upper=par2)$root;
			}
			if (up.bound == 1) {
				prof.par[i,4] <- 1;
			} else {
				prof.par[i,4] <- uniroot(thresholdE, lower=par2, upper=up.bound)$root;
			}
			if (prof.par[i,3] == par2) prof.par[i,3] <- 0; # precision problem?!? check optimization ***
		}
	}
	# remove columns that are not used
	idx <- which(apply(prof.par, MARGIN=2, function (x) all(is.na(x))));
	if (length(idx) > 0) {
		prof.par <- prof.par[,-idx];
	}
	prof.par <- as.data.frame(round(prof.par, digits=7));
	return(prof.par);
}
