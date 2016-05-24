######################################
## code for sample size simulations ##
## Nicholas Reich
## August 2010
######################################

## this is the central user function
## N = overall sample size
## med = median of log normal distribution
## disp = dispersion of log normal distribution
## percentile = what percentile(s) should be simulated
## nsim = how many simulations to run
## exact data = T/F should simulation be done with exact data
## pct.type.A = percent type A data, will be rounded up to nearest integer
## exp.win.dat = a vector of exposure window lengths to sample from
## verb = whether to print 10 iteration counts during the course of the sim


##' @name precision.simulation
##' @aliases precision.simulation.exact
##' @aliases precision.simulation.coarse
##' @aliases generate.coarse.data
##'   
##' @title Simulate incubation period analyses with coarse data
##'   
##' @description These functions simulate coarse incubation period data sets and
##'   analyze them.  The goal is for these simulations to provide evidence for
##'   how much information a given dataset contains about a characteristic of
##'   the incubation period distribution.
##'   
##' @param N Overall sample size for the datasets to be simulated.
##' @param med Median for the assumed log normal distribution of the incubation
##'   periods.
##' @param disp Dispersion for the assumed log normal distribution of the
##'   incubation periods.
##' @param percentile Percentile of the incubation period distribution which we
##'   want to estimate.
##' @param nsim Number of datasets to analyze in the simulation.
##' @param exact.data Either TRUE/FALSE.  Incidates whether the data generated
##'   should be coarsened at all.  If TRUE, pct.type.A and exp.win.dat are
##'   ignored.
##' @param pct.type.A Percent of the N observations that are assumed to be type
##'   A data.  If N*pct.type.A is not an integer, it will be rounded to the
##'   nearest integer.
##' @param exp.win.dat A vector of exposure window lengths.  Defaults to the
##'   observed window lengths from Lessler et al. (see below).
##' @param verb If TRUE, a message with the system time and iteration number
##'   will be printed ten times during the simulation run.
##'   
##'   
##' @rdname precision.simulation
##' @return The \code{precision.simulation} functions return a matrix with four
##'   columns and nsim rows.  The "ests" column gives the estimated percentiles
##'   for the incubation period distribution.  The "SE" column gives the
##'   standard error for the estimate.  The "conv" column is 1 if the doubly
##'   interval-censored likelihood maximization converged.  Otherwise, it is 0.
##'   The "bias" column gives the estimated percentile - true percentile. The
##'   \code{generate.coarse.data} function returns a matrix with data suitable
##'   for analysis by the \code{dic.fit} function.
##' @export
precision.simulation <- function(N,
			       med=2,
			       disp=1.3,
			       percentile=.5,
			       nsim=100,
			       exact.data=FALSE,
			       pct.type.A=.5,
 			       exp.win.dat=NULL,
			       verb=FALSE) {
	## logic check
	if(percentile <= 0 | percentile >=1)
		stop("percentile must be between 0 and 1.")
	if(pct.type.A < 0 | pct.type.A >1)
		stop("% of data that is type A must be between 0 and 1.")
	if(is.null(exp.win.dat)){
		if(verb) message("NYC exposure window data used")
		exp.win.dat <- get(data(exp.win.lengths, envir = environment()))
	}

	## TODO: add default exp.win.dat to package and load if NULL

	if(exact.data) {
		out <- precision.simulation.exact(N=N,
						med=med,
						disp=disp,
						percentile=percentile,
						nsim=nsim,
						verb=verb)
	} else {
		out <- precision.simulation.coarse(N=N,
						 med=med,
						 disp=disp,
						 percentile=percentile,
						 nsim=nsim,
						 pct.type.A=pct.type.A,
						 exp.win.dat=exp.win.dat,
						 verb=verb)
	}

	target <- qlnorm(percentile, log(med), log(disp))
	bias <- out[,"ests"]-target
	out <- cbind(out, bias)

	return(out)
}

##' @rdname precision.simulation
##' @export
precision.simulation.exact <- function(N,
                                       med,
                                       disp,
                                       percentile,
                                       nsim,
                                       verb) {
    storage <- matrix(NA, ncol=3, nrow=nsim)
    colnames(storage) <- c("ests", "SE", "conv")

    data <- matrix(0, ncol=6, nrow=nsim*N)
    colnames(data) <- c("dataset.id", "EL", "ER", "SL", "SR", "type")
    data[,"dataset.id"] <- rep(1:nsim, each=N)
    data[,"SL"] <- rlnorm(N*nsim, meanlog=log(med), sdlog=log(disp))
    data[,"SR"] <- data[,"SL"]
    data[,"type"] <- 2

    for(i in 1:nsim){
        tmp.dat <- data[which(data[,"dataset.id"]==i),]
        tmp.fit <- dic.fit(tmp.dat, ptiles=percentile)
        if(tmp.fit$conv==1){
            row.name <- paste("p", round(percentile*100), sep="")
            which.row <-
                which(rownames(tmp.fit$ests)==row.name)[1]
            storage[i,c("ests", "SE")] <-
                tmp.fit$ests[which.row, c("est", "StdErr")]
        } else {
            storage[i,c("ests", "SE")] <- NA
        }
        storage[i,"conv"] <- tmp.fit$conv
        if(verb & i%%(round(nsim/10))==0)
			print(paste("iteration",i,"complete ::", Sys.time()))

    }
    return(storage)
}

##' @rdname precision.simulation
##' @export
precision.simulation.coarse <- function(N,
				      med,
				      disp,
				      percentile,
				      nsim,
				      pct.type.A,
				      exp.win.dat,
				      verb) {
	## create storage
	storage <- matrix(NA, ncol=3, nrow=nsim)
	colnames(storage) <- c("ests", "SE", "conv")

	for(i in 1:nsim){
		tmp.dat <- generate.coarse.data(N=N,
						med=med,
						disp=disp,
						pct.type.A=pct.type.A,
						exp.win.dat=exp.win.dat)
		tmp.fit <- dic.fit(tmp.dat, ptiles=percentile)
		if(tmp.fit$conv==1){
			row.name <- paste("p", round(percentile*100), sep="")
			which.row <-
                            which(rownames(tmp.fit$ests)==row.name)[1]
			storage[i,c("ests", "SE")] <-
				tmp.fit$ests[which.row, c("est", "StdErr")]
		} else {
			storage[i,c("ests", "SE")] <- NA
		}
		storage[i,"conv"] <- tmp.fit$conv
		if(verb & i%%(round(nsim/10))==0)
			print(paste("iteration",i,"complete ::", Sys.time()))
	}
	return(storage)
}


##' @rdname precision.simulation
##' @export
generate.coarse.data <- function(N, med, disp, pct.type.A, exp.win.dat) {

 	n.type.A <- round(N*pct.type.A)
	n.type.B <- N-n.type.A


	E <- runif(N, 10, 11)
	T <- rlnorm(N, log(med), log(disp))
	S <- T + E
	SR <- ceiling(S)
	SL <- floor(S)

	## generate window types
	##   0 = short with bounded ER = type A
	##   1 = "long" with no ER = type B
	win.type <- rep(0:1, times=c(n.type.A, n.type.B))

	## generate window lengths,
	potential.lengths <- sample(exp.win.dat, size=N, replace=TRUE)
	win.length <- 1*(win.type==0) + potential.lengths*(win.type>0)

	## fix data
	ER <- ceiling(E)*(win.type==0) + SR*(win.type==1)
	EL <- pmin(ER-win.length, floor(E))

	## return the data
	cbind(EL, E, ER, SL, S, SR, win.length, win.type, type=0)
}
