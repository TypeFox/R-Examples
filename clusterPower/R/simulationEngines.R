## all power functions have
## OUTPUTS
##   ~ results matrix with cols
##          "dataset", "beta.est", "beta.cil", "beta.cih", "reject.null", "pval.permute"


power.sim.normal <- function(n.sim=10,
			     effect.size,
			     alpha=0.05,
			     n.clusters,
			     n.periods,
			     cluster.size,
			     btw.clust.var,
			     indiv.var=NULL,
			     ICC=NULL,
			     period.effect,
			     period.var,
			     estimation.function,
			     permute=FALSE,
			     verbose=FALSE)
{
	## validation
	incl.period.effect <- TRUE
	if(period.var<0)
		stop("period.var must be positive")
	if(is.null(ICC) & is.null(indiv.var))
		stop("Either ICC or indiv.var must be specified.")
	if(length(period.effect) != n.periods & length(period.effect) != 1)
		stop("Length of period.effect must equal n.periods or 1.")
	if(length(period.effect) > 1 & period.var > 0)
		stop("cannot specify both period variance and exact period effects.")
	if(period.var==0 & length(period.effect)==1){
		## message("**MSG**: Period effect taken as a constant across entire study.")
		incl.period.effect <- FALSE
	}
	if(n.clusters > 10 & grepl("fixed.effect", as.character(substitute(estimation.function)))) {
		warning("Power estimates from fixed effect models with large numbers of clusters may be unstable.")
	}

	## calculate total number of observations in entire study
	if(length(cluster.size)==1){
		## message("**MSG**: Treating cluster.size as the same for all cluster-periods.")
		cluster.size <- rep(cluster.size, n.clusters*n.periods)
	} else {
		## message("**MSG**: cluster.size values assumed fixed for a given cluster across all periods.")
	}
	n.obs <- sum(cluster.size)

	############################
	## set up data structures ##
	############################

	## dataset-level data
	results <- matrix(ncol=6, nrow=n.sim)
	colnames(results) <- c("dataset", "beta.est", "beta.cil", "beta.cih",
			       "reject.null", "pval.permute")
	results[,"dataset"] <- 1:n.sim

	## observation-level data
	tmp <- make.base.data(n.obs=n.obs,
			      n.clusters=n.clusters,
			      cluster.size=cluster.size,
			      n.periods=n.periods)
	sim.dat.base <- tmp$sim.dat.base
	design.mat <- tmp$design.mat

	message(paste("starting simulation ::", Sys.time()))
	for (i in 1:n.sim){
		sim.dat <- sim.dat.base

		## generate dataset
		clust.effects <- rnorm(n.clusters, mean=0, sd=sqrt(btw.clust.var))
		period.effect <- rnorm(n.periods,
				       mean=period.effect,
				       sd=sqrt(period.var))
		full.beta <- c(effect.size, clust.effects, period.effect)
		mean.y <- design.mat %*% full.beta
		sim.dat[,"mean.y"] <- mean.y
		if(is.null(indiv.var)) {
			indiv.var <- btw.clust.var*(1/ICC - 1)
		}
		noise <- rnorm(n.obs, 0, sd=sqrt(indiv.var))
		sim.dat[,"y"] <- mean.y + noise

		sim.dat <- data.frame(sim.dat, at.risk.time=1)
		sim.dat$clust <- factor(sim.dat$clust)
		sim.dat$per <- factor(sim.dat$per)


		## estimate tx effect
		tmp <- estimation.function(dat=sim.dat,
					   incl.period.effect=incl.period.effect,
					   outcome.type="gaussian",
					   alpha=alpha)
		results[i,c("beta.est", "beta.cil", "beta.cih")] <- tmp

		## make permutation inference about tx effect
		if(permute){
			## RUN PERMUTATION TEST
		}

		if(i%%10==0 & verbose) message(paste("iteration", i, "::", Sys.time()))
	}
	results[,"reject.null"] <- (sign(results[,"beta.cil"]) ==
				    sign(results[,"beta.cih"]))

	power <- sum(results[,"reject.null"])/n.sim
	if(verbose) message(paste("power =", round(power,3)))

	permute.power <- ifelse(permute, sum(results[,"pval.permute"]<alpha)/n.sim, NA)

	out <- list(results=results,
		    power=power,
		    permute.power=permute.power,
		    sample.data=sim.dat)
	return(out)
}


power.sim.binomial <- function(n.sim=10,
			     effect.size,
			     alpha=0.05,
			     n.clusters,
			     n.periods,
			     cluster.size,
			     btw.clust.var,
			     period.effect,
			     period.var,
			     estimation.function,
			     permute=FALSE,
			     verbose=FALSE)
{

	## validation
	incl.period.effect <- TRUE
	if(period.var<0)
		stop("period.var must be positive")
	if(length(period.effect) != n.periods & length(period.effect) != 1)
		stop("Length of period.effect must equal n.periods or 1.")
	if(length(period.effect) > 1 & period.var > 0)
		stop("cannot specify both period variance and exact period effects.")
	if(period.var==0 & length(period.effect)==1){
		## message("**MSG**: Period effect taken as a constant across entire study.")
		incl.period.effect <- FALSE
	}
	if(n.clusters > 10 & grepl("fixed.effect", as.character(substitute(estimation.function)))) {
		warning("Power estimates from fixed effect models with large numbers of clusters may be unstable.")
	}


	## calculate total number of observations in entire study
	if(length(cluster.size)==1){
		## message("**MSG**: Treating cluster.size as the same for all cluster-periods.")
		cluster.size <- rep(cluster.size, n.clusters*n.periods)
	} else {
		if(length(cluster.size)!=n.clusters)
			stop("cluster sizes not equal to number of clusters")
		## message("**MSG**: cluster.size values assumed fixed for a given cluster across all periods.")
	}
	n.obs <- sum(cluster.size)

	############################
	## set up data structures ##
	############################

	## dataset-level data
	results <- matrix(ncol=6, nrow=n.sim)
	colnames(results) <- c("dataset", "beta.est", "beta.cil", "beta.cih",
			       "reject.null", "pval.permute")
	results[,"dataset"] <- 1:n.sim

	## observation-level data
	tmp <- make.base.data(n.obs=n.obs,
			      n.clusters=n.clusters,
			      cluster.size=cluster.size,
			      n.periods=n.periods)
	sim.dat.base <- tmp$sim.dat.base
	design.mat <- tmp$design.mat

	message(paste("starting simulation ::", Sys.time()))
	for (i in 1:n.sim){
		sim.dat <- sim.dat.base

		## generate dataset
		clust.effects <- rnorm(n.clusters, mean=0, sd=sqrt(btw.clust.var))
		period.effect <- rnorm(n.periods,
				       mean=period.effect,
				       sd=sqrt(period.var))
		full.beta <- c(effect.size, clust.effects, period.effect)
		mean.y <- design.mat %*% full.beta
		sim.dat[,"mean.y"] <- mean.y
		sim.dat[,"y"] <- rbinom(nrow(sim.dat), size=1, prob=expit(mean.y))

		## using "at risk time" column to count cluster sizes
		sim.dat <- data.frame(sim.dat, at.risk.time=1)
		sim.dat$clust <- factor(sim.dat$clust)
		sim.dat$per <- factor(sim.dat$per)


		## estimate tx effect
		tmp <- estimation.function(dat=sim.dat,
					   incl.period.effect=incl.period.effect,
					   outcome.type="binomial",
					   alpha=alpha)
		results[i,c("beta.est", "beta.cil", "beta.cih")] <- tmp

		## make permutation inference about tx effect
		if(permute){
			## RUN PERMUTATION TEST
		}

		if(i%%10==0 & verbose) message(paste("iteration", i, "::", Sys.time()))
	}
	results[,"reject.null"] <- (sign(results[,"beta.cil"]) ==
				    sign(results[,"beta.cih"]))

	power <- sum(results[,"reject.null"])/n.sim
	if(verbose) message(paste("power =", round(power,3)))

	permute.power <- ifelse(permute, sum(results[,"pval.permute"]<alpha)/n.sim, NA)

	out <- list(results=results,
		    power=power,
		    permute.power=permute.power,
		    sample.data=sim.dat)
	return(out)
}


power.sim.poisson <- function(n.sim=10,
			      effect.size,
			      alpha=0.05,
			      n.clusters,
			      n.periods,
			      cluster.size,
			      btw.clust.var,
			      period.effect,
			      period.var,
			      estimation.function,
			      at.risk.params,
			      permute=FALSE,
			      verbose=FALSE)
{

	## validation
	incl.period.effect <- TRUE
	if(period.var<0)
		stop("period.var must be positive")
	if(length(period.effect) != n.periods & length(period.effect) != 1)
		stop("Length of period.effect must equal n.periods or 1.")
	if(length(period.effect) > 1 & period.var > 0)
		stop("cannot specify both period variance and exact period effects.")
	if(period.var==0 & length(period.effect)==1){
		## message("**MSG**: Period effect taken as a constant across entire study.")
		incl.period.effect <- FALSE
	}
	if(length(at.risk.params)==1)
		## message("**MSG**: at risk time assumed to be constant for all cluster-periods.")
	if(length(at.risk.params)!=1 & length(at.risk.params)!=2)
		stop("at.risk.params must have length 1 or 2")
	if(n.clusters > 10 & grepl("fixed.effect", as.character(substitute(estimation.function)))) {
		warning("Power estimates from fixed effect models with large numbers of clusters may be unstable.")
	}



	## calculate total number of observations in entire study
	if(length(cluster.size)==1){
		## message("**MSG**: Treating cluster.size as the same for all cluster-periods.")
		cluster.size <- rep(cluster.size, n.clusters*n.periods)
	} else {
		## message("**MSG**: cluster.size values assumed fixed for a given cluster across all periods.")
		cluster.size <- rep(cluster.size, each=n.periods)
	}
	n.obs <- sum(cluster.size)

	############################
	## set up data structures ##
	############################

	## dataset-level data
	results <- matrix(ncol=6, nrow=n.sim)
	colnames(results) <- c("dataset", "beta.est", "beta.cil", "beta.cih",
			       "reject.null", "pval.permute")
	results[,"dataset"] <- 1:n.sim

	## observation-level data
	tmp <- make.base.data(n.obs=n.obs,
			      n.clusters=n.clusters,
			      cluster.size=cluster.size,
			      n.periods=n.periods)
	sim.dat.base <- tmp$sim.dat.base
	design.mat <- tmp$design.mat

	message(paste("starting simulation ::", Sys.time()))
	for (i in 1:n.sim){
		sim.dat <- sim.dat.base

		## generate dataset
		clust.effects <- rnorm(n.clusters, mean=0, sd=sqrt(btw.clust.var))
		period.effect <- rnorm(n.periods,
				       mean=period.effect,
				       sd=sqrt(period.var))
		if(length(at.risk.params)==1){
			at.risk.time <- rep(at.risk.params, n.obs)
		} else {
			at.risk.time <- 1 + rnbinom(n.obs,
						    size=at.risk.params[2],
						    mu=at.risk.params[1])
		}


		full.beta <- c(effect.size, clust.effects, period.effect)
		mean.y <- design.mat %*% full.beta + log(at.risk.time)
		sim.dat[,"mean.y"] <- mean.y
		sim.dat[,"y"] <- rpois(nrow(sim.dat), exp(mean.y))

		sim.dat <- data.frame(sim.dat, at.risk.time=at.risk.time)
		sim.dat$clust <- factor(sim.dat$clust)
		sim.dat$per <- factor(sim.dat$per)


		## estimate tx effect
		tmp <- estimation.function(dat=sim.dat,
					   incl.period.effect=incl.period.effect,
					   outcome.type="poisson",
					   alpha=alpha)
		results[i,c("beta.est", "beta.cil", "beta.cih")] <- tmp

		## make permutation inference about tx effect
		if(permute){
			## RUN PERMUTATION TEST
		}

		if(i%%10==0 & verbose) message(paste("iteration", i, "::", Sys.time()))
	}
	results[,"reject.null"] <- (sign(results[,"beta.cil"]) ==
				    sign(results[,"beta.cih"]))

	power <- sum(results[,"reject.null"])/n.sim
	if(verbose) message(paste("power =", round(power,3)))

	permute.power <- ifelse(permute, sum(results[,"pval.permute"]<alpha)/n.sim, NA)

	out <- list(results=results,
		    power=power,
		    permute.power=permute.power,
		    sample.data=sim.dat)
	return(out)
}



make.base.data <- function(n.obs, n.clusters,
			   cluster.size, n.periods){
	## #####################
	## make base data     ##
	## #####################
	sim.dat.base <- matrix(ncol=6, nrow=n.obs)
	colnames(sim.dat.base) <- c("id", "clust", "per",
				    "trt", "mean.y", "y")
	sim.dat.base[,"id"] <- 1:n.obs

	## NOTE: length(cluster.size) should = n.periods*n.clusters,
	##       where cluster.size[1] = size of cluster 1 during period 1
	##       where cluster.size[2] = size of cluster 1 during period 2
	##       where cluster.size[3] = size of cluster 2 during period 1 (assuming 2 periods)
	##       ...
	sim.dat.base[,"clust"] <- rep(1:n.clusters, times=cluster.size, each=n.periods)

	sim.dat.base[,"per"] <- rep(rep(1:n.periods, times=n.clusters), times=cluster.size)

	## assign treatment times
	clust.list <- 1:n.clusters
	clust.seq.A <- sample(clust.list, round(n.clusters/2), replace=FALSE)
	clust.seq.B <- clust.list[which(!( clust.list %in% clust.seq.A))]
	trt.seq.A <- seq(1, n.periods, by=2)
	trt.seq.B <- ifelse(n.periods==1, 0, seq(2, n.periods, by=2))

	## generate index of cluster-periods with treatment
	trt.idx.A <- which(sim.dat.base[,"clust"] %in% clust.seq.A &
			   sim.dat.base[,"per"] %in% trt.seq.A)
	trt.idx.B <- which(sim.dat.base[,"clust"] %in% clust.seq.B &
			   sim.dat.base[,"per"] %in% trt.seq.B)

	## assign 1s indicating treatment received
	sim.dat.base[,"trt"] <- 0
	sim.dat.base[trt.idx.A,"trt"] <- 1
	sim.dat.base[trt.idx.B,"trt"] <- 1

	## #####################
	## make model matrix  ##
	## #####################
	X <- matrix(0, ncol=n.clusters+n.periods+2, nrow=nrow(sim.dat.base))
	colnames(X) <- c("id",
			 paste("clust.", 1:n.clusters, sep=""),
			 paste("per.", 1:n.periods, sep=""),
			 "trt")
	X[,"id"] <- sim.dat.base[,"id"]
	X[,"trt"] <- sim.dat.base[,"trt"]

	## loop through periods and clusters to create dummy variables
	for(i in 1:n.clusters){
		tmp.col <- paste("clust.", i, sep="")
		X[which(sim.dat.base[,"clust"]==i),tmp.col] <- 1
	}
	for(i in 1:n.periods){
		tmp.col <- paste("per.", i, sep="")
		X[which(sim.dat.base[,"per"]==i),tmp.col] <- 1
	}

	trt.col <- which(substr(colnames(X), 1, 3)=="trt")
	clust.cols <- which(substr(colnames(X), 1, 5)=="clust")
	per.cols <- which(substr(colnames(X), 1, 3)=="per")
	design.mat <- X[,c(trt.col, clust.cols, per.cols)]

	out <- list(design.mat=design.mat,
		    sim.dat.base=sim.dat.base)
	return(out)

}


logit <- function(p) log(p/(1-p))
expit <- function(x) exp(x)/(1+exp(x))

## converts mixed effect model parameters to hayes k parameters, for log-linear mixed effect models only
mixed.eff.params <- function(pi0, btw.clust.var, Tk) {
        e <- Tk * exp(pi0) * exp(btw.clust.var/2)
        v <- e + e^2*(exp(btw.clust.var)-1)
        k <- sqrt(v)/e
        return(c(expectation=e, variance=v, hayes.k=k))
}

## implements power calculation based on Hayes (1999) formulas and the coef of variation, k
hayes.power.poisson <- function(n.clusters, period.effect, btw.clust.var, at.risk.params, cluster.size, effect.size, alpha=.05) {
	Tk <- at.risk.params*cluster.size
        z.a <- qnorm(alpha/2, lower.tail=FALSE)
        l0 <- exp(period.effect)
        l1 <- exp(period.effect)*exp(effect.size)
        ## calculate Hayes metrics
        me.pars <- mixed.eff.params(pi0=period.effect, btw.clust.var=btw.clust.var, Tk=Tk)
        obs.k <- me.pars["hayes.k"]
        ## from formula 2 in Hayes et al.
        z.b <- sqrt((n.clusters-1) * (l0-l1)^2 / ( (l0+l1)/Tk + obs.k^2 * (l0^2+l1^2)) )-z.a
	beta <- unname(pnorm(z.b, lower.tail=FALSE))
        return(1-beta)
}
