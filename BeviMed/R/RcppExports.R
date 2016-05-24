#' R interface to BeviMed c++ MCMC procedure 
#'
#' Allows other functions in the package to call the c++ function passing arguments more succinctly and by name.
#'
#' @template samples_per_chain
#' @param y Logical vector of subject affectedness status.
#' @param block_starts Integer vector of k 0-indexed start positions (with respect to \code{cases} and \code{counts}) for contiguous blocks relating to the k variants.
#' @param block_ends Integer vector of (exclusive) k 0-indexed end positions.
#' @param cases 0 based vector of case indices with respect to y.
#' @param counts Vector of variant counts.
#' @param min_ac Minimum allele count required for pathogenic configuration.
#' @param q_shape Beta distribution parameterisation of benign variant configuration rate of affection, q.
#' @param p_shape Beta distribution parameterisation of pathogenic variant configuration rate of affection, p.
#' @param omega_shape Beta distribution of global rate of pathogenicty of variants in gene given pathogenicity of gene, omega.
#' @template temperatures
#' @param Z0_matrix Matrix of logicals, where the rows are used as an initial Zs for the chains.
#' @param estimate_omega Logical value determining whether to estimate the parameter omega.
#' @param logit_omegas Numeric vector of logit omega values, one value per chain.
#' @param logit_omega_proposal_sds Numeric vector of proposal standard deviations for Metropolis-Hastings sampling of logit omega parameter, one value per chain.
#' @template variant_weights
#' @template estimate_phi
#' @param log_phis Numeric vector of log phi values, one value per chain.
#' @template log_phi_mean
#' @template log_phi_sd
#' @param log_phi_proposal_sds Numeric vector of proposal standard deviations for Metropolis-Hastings sampling of log phi parameter, one value per chain.
#' @param chain_swaps_per_cycle Number of chain swaps to propose per update cycle.
#' @param annealing Logical value determining whether to anneal the chains, e.g. for optimisation.
#' @template tandem_variant_updates
#' @param case_variant_block_starts 0-indexed start positions for contiguous blocks of variants in \code{case_variants}.
#' @param case_variant_block_ends As \code{case_variant_block_starts} for (exclusive) stop positions.
#' @param case_variants Integer vector giving variant numbers (0-based, i.e. between 0 and k-1). Used to pick pairs of variants for tandem updates from.
#' @template store_Z
#' @template burn
#' @param check Logical value indicating whether to perform validation on the arguments before calling the c++ function. 
#' @return Object of class \code{BeviMed}, containing the output of the MCMC sampling.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib BeviMed
call_cpp <- function(
	samples_per_chain,
	y,
	block_starts,
	block_ends,
	cases,
	counts,
	min_ac,
	q_shape,
	p_shape,
	omega_shape,
	temperatures,
	Z0_matrix,
	estimate_omega,
	logit_omegas,
	logit_omega_proposal_sds,
	variant_weights,
	estimate_phi,
	log_phis,
	log_phi_mean,
	log_phi_sd,
	log_phi_proposal_sds,
	chain_swaps_per_cycle,
	annealing,
	tandem_variant_updates,
	case_variant_block_starts,
	case_variant_block_ends,
	case_variants,
	store_Z,
	burn=0,
	check=TRUE
) {
	if (check) {
		stopifnot(length(omega_shape) == 2 & min(omega_shape) > 0)
		stopifnot(length(q_shape) == 2 & min(q_shape) > 0)
		stopifnot(length(p_shape) == 2 & min(p_shape) > 0)
		stopifnot(is.matrix(Z0_matrix))
		stopifnot(length(temperatures) == nrow(Z0_matrix))
		stopifnot(length(temperatures) == nrow(logit_omegas))
		stopifnot(length(temperatures) == nrow(log_phis))
		stopifnot(length(temperatures) == length(logit_omega_proposal_sds))
		stopifnot(length(temperatures) == length(log_phi_proposal_sds))
		stopifnot(estimate_phi | all(log_phis == 0))
		stopifnot(length(block_ends) == length(block_starts))
		stopifnot(block_ends[length(block_ends)] == length(cases))
		stopifnot(length(cases) == length(counts))
		if (tandem_variant_updates > 0 & length(unique(case_variants)) < 2) stop("Must have more than 1 variant to select from if making tandem updates")
	}

	raw <- .Call(
		"R_parallel_tempered_markov_chain",
		samples_per_chain,
		y,
		block_starts,
		block_ends,
		cases,
		counts,
		min_ac,
		q_shape[1],
		q_shape[2],
		p_shape[1],
		p_shape[2],
		omega_shape[1],
		omega_shape[2],
		Z0_matrix,
		estimate_omega,
		logit_omegas,
		logit_omega_proposal_sds,
		variant_weights,
		estimate_phi,
		log_phis,
		log_phi_mean,
		log_phi_sd,
		log_phi_proposal_sds,
		temperatures,
		chain_swaps_per_cycle,
		annealing,
		tandem_variant_updates,
		case_variant_block_starts,
		case_variant_block_ends,
		case_variants,
		store_Z,
		PACKAGE="BeviMed"
	)

	if (burn > 0) {
		raw <- c(
			lapply(
				raw[c(
					"y_log_lik",
					"y_log_lik_t_equals_1",
					"Z",
					"logit_omega",
					"log_phi"
				)],
				function(x) x[-seq(length.out=burn),]
			),
			lapply(
				raw[c("swap_accept", "swap_at_temperature")],
				function(x) x[-seq(length.out=chain_swaps_per_cycle * burn)]
			),
			raw[c("terminal_Z", "terminal_log_phi", "terminal_logit_omega")]
		)
	}

	structure(
		c(
			raw,
			list(
				omega_shape=omega_shape,
				p_shape=p_shape,
				q_shape=q_shape,
				variant_weights=variant_weights,
				temperatures=temperatures,
				estimate_phi=estimate_phi,
				estimate_omega=estimate_omega,
				y=y,
				variant_table=data.frame(
					variant=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 1:length(block_starts), block_ends-block_starts)),
					case=cases+1,
					count=counts
				),
				n=length(y),
				k=ncol(Z0_matrix)
			)
		),
		class="BeviMed"
	)
}

#' Calculate the Marginal Likelihood under v by summation over power posterior likelihood exptectances
#'
#' @template y_log_lik_t_equals_1_traces
#' @param temperatures Numeric vector of temperatures used to produce \code{y_log_lik_t_equals_1_traces}.
#' @return Numeric value of estimated log marginal likelihood.
#' @export
power_posteriors_ML_sum <- function(y_log_lik_t_equals_1_traces, temperatures) {
	sum(mapply(FUN=function(y_lik, t_diff) { log(mean(exp(t_diff*(y_lik)))) }, split(t(y_log_lik_t_equals_1_traces)[-length(temperatures),], seq(length(temperatures)-1)), diff(temperatures)))
}

#' Concatenate objects of class \code{BeviMed}
#'
#' This function could be used to stitch together consecutive chains to create one larger sampled set of states from the MCMC procedure.
#' @param ... BeviMed objects
#' @return \code{BeviMed} object.
#' @importFrom stats setNames
stack_BeviMeds <- function(...) {
	bevis <- list(...)
	stopifnot(all(sapply(bevis, class) == "BeviMed"))

	structure(
		c(
			lapply(
				setNames(nm=c(
					"y_log_lik",
					"y_log_lik_t_equals_1",
					"Z",
					"logit_omega",
					"log_phi"
				)),
				function(trace_name) do.call(what=rbind, lapply(bevis, "[[", trace_name))
			),
			bevis[[length(bevis)]][-which(names(bevis[[length(bevis)]]) %in% c("y_log_lik", "y_log_lik_t_equals_1", "Z", "logit_omega", "log_phi"))]
		),
		class="BeviMed"
	)
}

#' Estimate confidence interval for estimated marginal likelihood by simulation
#'
#' @template temperatures
#' @template y_log_lik_t_equals_1_traces
#' @template confidence
#' @template simulations
#' @export
#' @importFrom stats rt quantile
estimate_confidence_interval <- function(
	temperatures,
	y_log_lik_t_equals_1_traces,
	confidence=0.95,
	simulations=1000
) {
	#we want all blocks to have the same length
	a <- b <- as.integer(sqrt(nrow(y_log_lik_t_equals_1_traces)))
	if (a < 2) stop("Longer sample-block lengths are required to make confidence interval estimation")

	#with columns for batches and rows for temperatures...
	batch_means <- simplify2array(tapply(X=(nrow(y_log_lik_t_equals_1_traces)-a*b+1):nrow(y_log_lik_t_equals_1_traces), INDEX=gl(n=a, k=b), FUN=function(rows) apply(MARGIN=1, FUN=mean, X=exp(t(y_log_lik_t_equals_1_traces[rows,-ncol(y_log_lik_t_equals_1_traces),drop=FALSE]) * diff(temperatures)))))

	overall_mean <- apply(MARGIN=1, FUN=mean, X=exp(t(y_log_lik_t_equals_1_traces[,-ncol(y_log_lik_t_equals_1_traces)]) * diff(temperatures)))

	estimate_var_by_temp <- mapply(FUN=function(overall, batch) b / (a-1) * sum((batch-overall)^2), overall_mean, split(batch_means, seq(nrow(batch_means))))

	samples <- mapply(SIMPLIFY=TRUE, FUN=function(est_mean, est_var) est_mean + (rt(df=a-1, n=simulations) * sqrt(est_var) / sqrt(a)), overall_mean, estimate_var_by_temp)

	simulated_MLs <- apply(samples, 1, function(expectance_at_temp_i) sum(log(ifelse(expectance_at_temp_i < 0, 0, ifelse(expectance_at_temp_i > 1, 1, expectance_at_temp_i)))))

	quantile(probs=c((1-confidence)/2,1-(1-confidence)/2), simulated_MLs)
}

#' Apply the MCMC algorithm in blocks until conditions are met
#'
#' Sample blocks of a given size until either the estimated log marginal likelihood falls within a given confidence interval, there is sufficient confidence that the log Bayes factor over the null model is at most a certain quantity, or a certain number of blocks have been sampled.
#' @template y
#' @param blocks_remaining Maximum number of blocks left before termination.
#' @param start_Zs Initial (logical) Z-matrix.
#' @param start_logit_omegas Initial values of logit_omega (numeric vector - one value per chain).
#' @param start_log_phis Initial values of log_phi (numeric vector - one value per chain).
#' @template temperatures
#' @param tolerance Maximum width for confidence_interval of log marginal likelihood to allow before stopping the chain.
#' @template confidence
#' @template simulations
#' @template marginal_likelihood_only
#' @param quit_if_highest_possible_BF_less_than Numeric value used to determine whether to stop the sampling after consecutive blocks. If we are confident (to the extent given by the parameter \code{confidence}) that log Bayes factor of v against n is under this value, we stop the sampling as soon as possible.
#' @template y_log_lik_t_equals_1_traces
#' @param full_block_traces List of outputs of calls to MCMC routine. 
#' @template verbose
#' @param ... Other arguments passed to \code{\link{call_cpp}}
#' @return If \code{marginal_likelihood_only == TRUE}, then a numeric value giving the log marginal likelihood, otherwise an object of class \code{BeviMed}.
stop_chain <- function(
	y,
	blocks_remaining,
	start_Zs,
	start_logit_omegas,
	start_log_phis,
	temperatures,
	tolerance=1,
	confidence=0.95,
	simulations=1000,
	marginal_likelihood_only=TRUE,
	quit_if_highest_possible_BF_less_than=-Inf,
	y_log_lik_t_equals_1_traces=matrix(ncol=length(temperatures),nrow=0),
	full_block_traces=list(),
	verbose=FALSE,
	...
) {
	if (verbose) cat("Sampling up to ", blocks_remaining, " more blocks to get marginal likelihood within tolerance of ", tolerance, "\n", sep="\n")

	mc <- call_cpp(
		Z0_matrix=start_Zs,
		logit_omegas=start_logit_omegas,
		log_phis=start_log_phis,
		temperatures=temperatures,
		y=y,
		...
	)

	y_ll <- rbind(
		y_log_lik_t_equals_1_traces,
		mc[["y_log_lik_t_equals_1"]]
	)
	
	confidence_interval <- estimate_confidence_interval(temperatures, y_ll, confidence=confidence, simulations=simulations)

	if (verbose) {
		cat("Confidence interval:\n")
		print(round(digits=2, confidence_interval))
	}

	if (diff(confidence_interval) < tolerance | (confidence_interval[2]-`log P_n`(y)) < quit_if_highest_possible_BF_less_than | blocks_remaining <= 1) {
		if (marginal_likelihood_only) { 
			power_posteriors_ML_sum(y_ll, temperatures) 
		} else {
			do.call(what=stack_BeviMeds, c(full_block_traces, list(mc)))
		}
	} else {
		stop_chain(
			y=y,
			blocks_remaining=blocks_remaining - 1,
			tolerance=tolerance,
			confidence=confidence,
			start_Zs=mc[["terminal_Z"]],
			start_logit_omegas=mc[["terminal_logit_omega"]],
			start_log_phis=mc[["terminal_log_phi"]],
			marginal_likelihood_only=marginal_likelihood_only,
			temperatures=temperatures,
			quit_if_highest_possible_BF_less_than=quit_if_highest_possible_BF_less_than,
			y_log_lik_t_equals_1_traces=y_ll,
			full_block_traces=c(full_block_traces, list(mc)),
			...
		)
	}
}

#' Tune the proposal standard deviations for the Metropolis-Hastings updates of either phi or omega
#'
#' @param tune_for Character vector, where only the first element is used, determining which variable to tune the proposal SDs for.
#' @param initial_proposal_sds Numeric vector with the initial values of the proposal SDs.
#' @param target_acceptance_range Numeric vector of length 2 where the first element is the lower bound for the acceptance interval and the second is the upper bound.
#' @param other_param_proposal_sd The proposal SD to use for \code{log_phi} when tuning \code{logit_omega} or vice versa.
#' @param max_tuning_cycles Maximum number of tuning cycles to perform before returning the proposal SDs as they are.
#' @param initial_rate Initial rate at which to mutate the proposal SDs.
#' @param rate_decay Geometric rate of decay for size of proposal SD mutation with each successive tuning cycle.
#' @template verbose
#' @param ... Other arguments to be passed to \code{\link{call_cpp}}.
#' @return Numeric vector of proposal SDs for the different temperature chains.
#' @export
tune_proposal_sds <- function(tune_for=c("logit_omega", "log_phi"), initial_proposal_sds, target_acceptance_range=c(0.3,0.7), other_param_proposal_sd=0.7, max_tuning_cycles=10, initial_rate=1, rate_decay=1.2, verbose=FALSE, ...) {
	stopifnot(all(tune_for %in% c("logit_omega", "log_phi")))

	if (verbose) {
		cat("Tuning proposal standard deviations for ", tune_for[1], " targeting acceptance rate range (", target_acceptance_range[1], ",", target_acceptance_range[2], ")\n", sep="")
	}

	current_proposal_sds <- initial_proposal_sds
	acceptances <- rep(-Inf, length(current_proposal_sds))
	cycle <- 0
	while (
		cycle <= max_tuning_cycles 
		& ( any(acceptances < target_acceptance_range[1])
		   |any(acceptances > target_acceptance_range[2]))
	) {
		if (verbose) cat("Tuning cycle ", cycle, "\n\tTest proposal SDs:\n\t\t", paste0(collapse=" : ", round(digits=2, current_proposal_sds)), "\n", sep="")
		out <- if (tune_for[1] == "logit_omega") call_cpp(logit_omega_proposal_sds=current_proposal_sds, log_phi_proposal_sds=rep(other_param_proposal_sd, length(current_proposal_sds)), ...) else call_cpp(log_phi_proposal_sds=current_proposal_sds, logit_omega_proposal_sds=rep(other_param_proposal_sd, length(current_proposal_sds)), ...)

		acceptances <- apply(out[[tune_for[1]]], 2, function(var_vals) mean(var_vals[-length(var_vals)] != var_vals[-1]))
		current_proposal_sds <- mapply(SIMPLIFY=TRUE, FUN=function(prop_sd, acc_rate) if (acc_rate < target_acceptance_range[1] | acc_rate > target_acceptance_range[2]) { match.fun(if (acc_rate < target_acceptance_range[1]) "/" else "*")(prop_sd, (1 + initial_rate * rate_decay ^ (-cycle+1))) } else { prop_sd }, current_proposal_sds, acceptances)

		if (verbose) cat("\tAcceptance rates:\n\t\t", paste0(collapse=" : ", round(digits=2, acceptances)), "\n", sep="")
	}
	if (verbose) cat("Terminating\n")
	current_proposal_sds
}

#' Tune temperatures using interval bisection to minimimise Kullback-Liebler divergence between adjacent power posteriors
#'
#' @param number_of_temperatures Integer value giving number of tuned temperatures (including 0 and 1) to obtain.
#' @param return_temperatures Logical value determining whether to return just the numeric vector of tuned temperatures or to return the \code{BeviMed}-classed object containing the output of the MCMC sampling.
#' @param ... Other arguments to pass to \code{call_cpp}.
#' @return If \code{return_temperatures == TRUE}, a numeric vector of tuned temperatures, otherwise an object of class \code{BeviMed}.
#' @export
#' @importFrom stats var
tune_temperatures <- function(
	number_of_temperatures,
	return_temperatures=FALSE,
	...
) {
#could add option for parallel tempering the tuning runs... i.e. discard all but last - obviously less efficient as in this case we won't be using the other temperature chains...
	temperatures <- 0:1
	chains <- lapply(temperatures, function(t) call_cpp(
		temperatures=t,
		...
	))

	E <- sapply(chains, function(ch) mean(ch[["y_log_lik_t_equals_1"]]))
	V <- sapply(chains, function(ch) var(ch[["y_log_lik_t_equals_1"]]))

	while (length(temperatures) <= number_of_temperatures) {
		areas <- diff(temperatures) * diff(E)
		largest <- which.max(areas)
		t_pair <- c(largest, largest+1)
		t_intersect <- mean(temperatures[t_pair])
		temperatures <- c(temperatures[1:largest], t_intersect, temperatures[(largest+1):length(temperatures)])
		chains <- c(chains[1:largest], list(call_cpp(temperatures=t_intersect, ...)), chains[(largest+1):length(chains)])

		E <- sapply(chains, function(ch) mean(ch[["y_log_lik_t_equals_1"]]))
		V <- sapply(chains, function(ch) var(ch[["y_log_lik_t_equals_1"]]))
	}
	
	if (return_temperatures) temperatures else list(chains=chains, E=E, V=V, temperatures=temperatures)
}

#' Perform inference under model variant-level model v for y using MCMC sampling
#'
#' This function is the main user-interface to the underlying c++ MCMC sampling routine.
#' @template y
#' @template G_matrix
#' @template min_ac
#' @template q_shape
#' @template p_shape
#' @template omega_shape 
#' @template samples_per_chain
#' @param stop_early Logical value determining whether to attempt to stop the sampling as soon as certain conditions are met (i.e. either the estimated marginal log likelihood lies within a certain confidence interval, or we are sufficiently confidence that the log Bayes factor against the null model is sufficiently low).
#' @param max_samples Maximum number of samples to take before terminating the sampling.
#' @template burn
#' @template temperatures
#' @param tune_temperatures Integer value - if greater than 0, the \code{temperatures} argument is ignored, and instead \code{tune_temperatures} tuned temperatures are used instead.
#' @template marginal_likelihood_only
#' @template store_Z
#' @param swaps Number of swaps between adjacent tempered chains to perform per update cycle.
#' @param optimise_Z0 Logical value determining whether to use a simulated annealing optimisation run to tune the initial values of \code{Z}.
#' @param tune_omega_and_phi_proposal_sd Logical value determining whether the proposal SDs of the Metropolis-Hastings estimated parameters should be tuned for a target acceptance range (see \code{\link{tune_proposal_sds}},
#' @param tune_block_size Integer value giving number of samples to draw when estimatating the acceptance rate of the omega/phi proposals in call to \code{\link{tune_proposal_sds}}.
#' @template variant_weights
#' @template estimate_phi
#' @template log_phi_mean
#' @template log_phi_sd
#' @template tandem_variant_updates
#' @param ... Other arguments to be passed to \code{\link{stop_chain}} and/or \code{\link{tune_proposal_sds}}.
#' @return If \code{marginal_likelihood_only == TRUE}, then a numeric value giving the log marginal likelihood, otherwise an object of class \code{BeviMed}.
#' @export
#' @importFrom stats rnorm runif rbeta
bevimed <- function(
	y,
	G,
	min_ac=1L,
	q_shape=c(2,100),
	p_shape=c(10, 2),
	omega_shape=if (min_ac == 1L) c(2, 9) else c(2, 2),
	samples_per_chain=2000,
	stop_early=FALSE,
	max_samples=5000,
	burn=as.integer(samples_per_chain/10),
	temperatures=(0:10/10)^2,
	tune_temperatures=0,
	marginal_likelihood_only=FALSE,
	store_Z=!marginal_likelihood_only,
	swaps=as.integer(length(temperatures)/2),
	optimise_Z0=FALSE,
	tune_omega_and_phi_proposal_sd=FALSE,
	tune_block_size=100,
	variant_weights=NULL,
	estimate_phi=!is.null(variant_weights), 
	log_phi_mean=-1.5,
	log_phi_sd=1,
	tandem_variant_updates=if (min_ac == 1) 0 else sum(y),
	...
) {
	stopifnot(is.matrix(G))
	stopifnot(ncol(G)==length(y))
	stopifnot(is.numeric(G))
	stopifnot(is.logical(y))

	normalised_weights <- if (is.null(variant_weights)) integer(nrow(G)) else variant_weights - mean(variant_weights)

	counts <- as.integer(t(G))
	variants <- rep(1:nrow(G), each=ncol(G))
	cases <- rep(1:ncol(G), times=nrow(G))

	block_ends <- cumsum(lapply(split(counts, variants), function(cnts) sum(cnts > 0)))
	block_starts <- c(0, block_ends[-length(block_ends)])

	y1_cases_with_more_than_1_variant <- apply(G > 0, 2, sum) > 1
	y1_variants <- lapply(split(t(G[,y1_cases_with_more_than_1_variant] > 0), seq(sum(y1_cases_with_more_than_1_variant))), which)
	y1_block_ends <- cumsum(lapply(y1_variants, length))
	y1_block_starts <- c(0, y1_block_ends[-length(y1_block_ends)])
	adjusted_tvu <- if (sum(y1_cases_with_more_than_1_variant) > 0) tandem_variant_updates else 0

	initial_log_phi_proposal_sd <- 0.5
	initial_logit_omega_proposal_sd <- 1
	estimate_omega <- !is.null(variant_weights)

	reused_arguments <- list(
		y=y,
		block_starts=block_starts,
		block_ends=block_ends,
		cases=cases[counts > 0]-1,
		counts=counts[counts > 0],
		min_ac=min_ac,
		q_shape=q_shape,
		p_shape=p_shape,
		omega_shape=omega_shape,
		estimate_omega=estimate_omega,
		variant_weights=normalised_weights,
		estimate_phi=estimate_phi,
		log_phi_mean=log_phi_mean,
		log_phi_sd=log_phi_sd,
		chain_swaps_per_cycle=swaps,
		tandem_variant_updates=adjusted_tvu,
		case_variant_block_starts=y1_block_starts,
		case_variant_block_ends=y1_block_ends,
		case_variants=unlist(y1_variants)-1
	)

	if (tune_temperatures > 0) {
		temperatures <- do.call(
			what=BeviMed::tune_temperatures,
			c(
				reused_arguments,
				list(
					samples_per_chain=samples_per_chain,
					number_of_temperatures=tune_temperatures,
					return_temperatures=TRUE,
					Z0_matrix=matrix(runif(nrow(G)) < omega_shape[1]/sum(omega_shape), nrow=1, ncol=nrow(G)),
					logit_omegas=local({ w <- rbeta(n=1, shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					logit_omega_proposal_sds=rep(initial_logit_omega_proposal_sd, 1),
					log_phis=if (estimate_phi) rnorm(n=1, mean=log_phi_mean, sd=log_phi_sd) else rep(0, 1),
					log_phi_proposal_sds=rep(initial_log_phi_proposal_sd, 1), 
					annealing=FALSE,
					store_Z=FALSE
				)
			)
		)
	}

	initial_Z <- if (optimise_Z0) {
		do.call(what=call_cpp, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				Z0_matrix=matrix(runif(nrow(G)*length(temperatures)) < omega_shape[1]/sum(omega_shape), nrow=length(temperatures), ncol=nrow(G)),
				logit_omegas=local({ w <- rbeta(n=temperatures, shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
				logit_omega_proposal_sds=rep(initial_logit_omega_proposal_sd, length(temperatures)),
				log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
				log_phi_proposal_sds=rep(initial_log_phi_proposal_sd, length(temperatures)), 
				temperatures=rep(1, length(temperatures)),
				chain_swaps_per_cycle=swaps,
				annealing=TRUE,
				store_Z=FALSE
			)
		))[["terminal_Z"]]
	} else {
		matrix(runif(nrow(G) * length(temperatures)) < omega_shape[1]/sum(omega_shape), nrow=length(temperatures), ncol=nrow(G))
	}

	proposal_sds <- lapply(
		setNames(nm=c("logit_omega", "log_phi")), 
		FUN=if (tune_omega_and_phi_proposal_sd & (estimate_phi | estimate_omega)) { 
			function(tune_for) do.call(what=tune_proposal_sds, c(
				reused_arguments, 
				list(
					initial_proposal_sds=rep(if (tune_for == "log_phi") initial_log_phi_proposal_sd else initial_logit_omega_proposal_sd, length(temperatures)),
					samples_per_chain=tune_block_size,
					burn=0,
					Z0_matrix=initial_Z,
					logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
					temperatures=temperatures,
					annealing=FALSE,
					store_Z=FALSE
				),
				list(...)[intersect(names(list(...)), names(formals(tune_proposal_sds)))]
			))
		} else { 
			function(tune_for) rep(if (tune_for == "log_phi") initial_log_phi_proposal_sd else initial_logit_omega_proposal_sd, length(temperatures))
		})

	if (stop_early) {
		burn <- do.call(what=call_cpp, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				burn=0,
				Z0_matrix=initial_Z,
				logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
				logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
				log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
				log_phi_proposal_sds=proposal_sds[["log_phi"]], 
				temperatures=temperatures,
				annealing=FALSE,
				store_Z=FALSE
			)
		))

		do.call(what=stop_chain, c(
			reused_arguments,
			list(
				samples_per_chain=samples_per_chain,
				y_log_lik_t_equals_1=matrix(ncol=length(temperatures),nrow=0),
				burn=0,
				marginal_likelihood_only=marginal_likelihood_only,
				start_Zs=burn[["terminal_Z"]],
				start_logit_omegas=burn[["terminal_logit_omega"]],
				start_log_phis=burn[["terminal_log_phi"]],
				temperatures=temperatures,
				blocks_remaining=max(1, as.integer(max_samples/samples_per_chain)),
				logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
				log_phi_proposal_sds=proposal_sds[["log_phi"]], 
				annealing=FALSE,
				store_Z=store_Z
			),
			list(...)[intersect(names(list(...)), names(formals(stop_chain)))]
		))
	} else {
		chains <- do.call(
			what=call_cpp,
			c(
				reused_arguments,
				list(
					samples_per_chain=samples_per_chain,
					burn=burn,
					Z0_matrix=initial_Z,
					logit_omegas=local({ w <- rbeta(n=length(temperatures), shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
					logit_omega_proposal_sds=proposal_sds[["logit_omega"]],
					log_phis=if (estimate_phi) rnorm(n=length(temperatures), mean=log_phi_mean, sd=log_phi_sd) else rep(0, length(temperatures)),
					log_phi_proposal_sds=proposal_sds[["log_phi"]], 
					temperatures=temperatures,
					annealing=FALSE,
					store_Z=store_Z
				)
			)
		)

		if (!marginal_likelihood_only) chains else power_posteriors_ML_sum(chains[["y_log_lik_t_equals_1"]], temperatures)
	}
}

#' Generate samples from the posterior under the variant-level pathogenicity model
#' 
#' Samples the posterior of the variant-level pathogenicity model using a Gibbs sampling based MCMC routine.
#' 
#' @template y
#' @template G_matrix
#' @template min_ac
#' @template q_shape
#' @template p_shape
#' @template omega_shape
#' @template variant_weights
#' @template estimate_phi
#' @param samples Number of samples to draw.
#' @template burn
#' @param just_Z Logical value determining whether to only return the matrix of Z-samples.
#' @export
#' @importFrom stats runif rbeta rnorm
sample_posterior_v <- function(
	y,
	G,
	min_ac=1L,
	q_shape=c(2,100),
	p_shape=c(10, 2),
	omega_shape=c(2, 9),
	variant_weights=NULL,
	estimate_phi=!is.null(variant_weights),
	samples=10000,
	burn=as.integer(samples/10),
	just_Z=FALSE
) {
	counts <- as.integer(t(G))
	variants <- rep(1:nrow(G), each=ncol(G))
	cases <- rep(1:ncol(G), times=nrow(G))

	block_ends <- cumsum(lapply(split(counts, variants), function(cnts) sum(cnts > 0)))
	block_starts <- c(0, block_ends[-length(block_ends)])

	#note - we *must* estimate omega if we are to use weights as we can't integrate it out
	estimate_omega <- !is.null(variant_weights)

	out <- call_cpp(
		samples_per_chain=samples,
		y=y,
		block_starts=block_starts,
		block_ends=block_ends,
		cases=cases[counts > 0]-1,
		counts=counts[counts > 0],
		min_ac=min_ac,
		q_shape=q_shape,
		p_shape=p_shape,
		omega_shape=omega_shape,
		Z0_matrix=matrix(nrow=1, runif(n=nrow(G)) < omega_shape[1]/sum(omega_shape)),
		temperatures=1.0,
		estimate_omega=estimate_omega,
		logit_omegas=local({ w <- rbeta(n=1, shape1=omega_shape[1], shape2=omega_shape[2]); log(w)-log(1-w) }),
		logit_omega_proposal_sds=1.0,
		variant_weights=if (is.null(variant_weights)) integer(nrow(G)) else variant_weights - mean(variant_weights),
		estimate_phi=estimate_phi,
		log_phis=if (estimate_phi) rnorm(n=1, mean=-1.5, sd=1) else integer(1),
		log_phi_mean=-1.5,
		log_phi_sd=1,
		log_phi_proposal_sds=1,
		chain_swaps_per_cycle=0,
		annealing=FALSE,
		tandem_variant_updates=0,
		case_variant_block_starts=integer(0),
		case_variant_block_ends=integer(0),
		case_variants=integer(0),
		store_Z=TRUE
	)

	if (just_Z)
		"[["(out, "Z")
	else
		out[c("Z", "logit_omega", "log_phi")[c(TRUE, rep(estimate_omega, 2))]]
}

#' Calculate marginal probability of observed genotypes under 'null' model
#' 
#' Marginal probability calculated exactly by integration.
#' @template y
#' @template q_shape
#' @export
`log P_n` <- function(
	y,
	q_shape=c(2, 100)
) {
	lbeta(sum(y)+q_shape[1], length(y)-sum(y)+q_shape[2]) -
	lbeta(q_shape[1], q_shape[2])
}

#' Calculate log Bayes factor comparing the variant-level model for case-control status
#' 
#' @template y
#' @param ... Other arguments to pass to \code{\link{bevimed}}.
#' @return Log Bayes factor.
#' @export
`log BF` <- function(
	y,
	...
) {
	`log P_v`(
		y=y,
		...
	) - 
	`log P_n`(y)
}

#' Calculate region or variant-level probability of pathogencity given prior probability on pathogenicity of region
#'
#' @param prior Numeric value giving prior for probability of pathogenicity of region.
#' @param variant_level Logical value determining whether to return the probabilities of pathogenicity for the individual variants or for the region as a whole.
#' @template y
#' @template G_matrix
#' @template min_ac
#' @param ... Other arguments to pass to \code{\link{bevimed}}. 
#' @return Probabilities of pathogenicity.
#' @export
probability_pathogenic <- function(
	prior=0.05,
	variant_level=FALSE,
	y,
	G,
	min_ac=1,
	...
) {
	stopifnot(min_ac > 0)
	stopifnot(is.matrix(G) & is.numeric(G))
	stopifnot(ncol(G) == length(y))

	sampling_summary <- summary(bevimed(
		marginal_likelihood_only=FALSE,
		y=y,
		G=G,
		min_ac=min_ac,
		...
	))

	log_bf <- sampling_summary[["ML_v"]]-sampling_summary[["ML_n"]]

	posterior <- exp(log_bf) * prior / (exp(log_bf) * prior + 1 - prior)

	if (variant_level) {
		sampling_summary[["Z"]]	* posterior
	} else {
		posterior
	}
}

#' Calculate marginal probability of observed genotypes under variant-level pathogenicity model
#' 
#' @template y
#' @template G_matrix
#' @template min_ac
#' @param ... Other arguments to pass to \code{\link{bevimed}}.
#' @export
`log P_v` <- function(
	y,
	G,
	min_ac=1,
	...
) {
	stopifnot(min_ac > 0)
	stopifnot(is.matrix(G) & is.numeric(G))
	stopifnot(ncol(G) == length(y))

	bevimed(
		marginal_likelihood_only=TRUE,
		y=y,
		G=G,
		min_ac=min_ac,
		...
	)
}

#' Calculate marginal probability of observed genotypes under 'pathogenic region' model
#' 
#' Marginal probability calculated exactly by integration.
#' @template y
#' @param G Integer/logical vector of genotypes by individual corresponding to case-control label \code{y} giving the 'rare variant counts'/'presence of rare variant indicators'.
#' @template min_ac
#' @template q_shape
#' @template p_shape
#' @export
`log P_r` <- function(
	y,
	G,
	min_ac=1L,
	q_shape=c(2, 4),
	p_shape=c(1, 1)
) {
	stopifnot(is.vector(G))
	stopifnot(length(y) == length(G))
	stopifnot(is.logical(y))

	`G'` <- if (is.logical(G)) G else G >= min_ac

	ll_pathogenic <- 
		+ lbeta(sum(y[`G'`])+p_shape[1], sum(`G'`)-sum(y[`G'`])+p_shape[2])
		- lbeta(p_shape[1], p_shape[2])

	ll_benign <- 
		+ lbeta(sum(y[!`G'`])+q_shape[1], sum(!`G'`)-sum(y[!`G'`])+q_shape[2])
		- lbeta(q_shape[1], q_shape[2])

	ll_pathogenic+ll_benign
}

#' Calculate marginal probability of observed genotypes under variant-level pathogenicity model by summing over marginals condintional on pathogenic variant combinations
#' 
#' @template y
#' @template G_matrix
#' @template min_ac
#' @template q_shape
#' @template p_shape
#' @template omega_shape
#' @param sum_over_variants Subset of variants for whose power set to calculate the direct sum over.
#' @export
lower_bound_P_v_by_direct_sum <- function(
	y,
	G,
	min_ac=1,
	q_shape=c(2,100),
	p_shape=c(10, 2),
	omega_shape=c(2, 9),
	sum_over_variants=1:nrow(G)
) {
	sum(exp(apply(
		as.matrix(do.call(what=expand.grid, rep(list(c(F,T)), length(sum_over_variants)))),
		1,
		function(Z) {
			x_i <- apply(G[sum_over_variants[Z],,drop=FALSE], 2, function(variant_config) sum(variant_config) >= min_ac)
			sig11 <- sum(x_i & y)
			sig10 <- sum(x_i & !y)
			sig01 <- sum(!x_i & y)
			sig00 <- sum(!x_i & !y)

			(
				+ lbeta(sig11+p_shape[1], sig10+p_shape[2])
				+ lbeta(sig01+q_shape[1], sig00+q_shape[2])
				+ lbeta(sum(Z)+omega_shape[1], nrow(G)-sum(Z)+omega_shape[2])
			) -
			(
			 	+ lbeta(p_shape[1], p_shape[2])
			 	+ lbeta(q_shape[1], q_shape[2])
			 	+ lbeta(omega_shape[1], omega_shape[2])
			)
		}
	)))
}
