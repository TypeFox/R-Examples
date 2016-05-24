tune_proposal_sd_conf <- list(
	initial_rate=1,
	rate_decay=1.2,
	max_tuning_cycles=10,
	target_acceptance=c(0.3, 0.5),
	pre_burn=1000,
	cycle_length=500
)

sim_reg_parameters <- c("gamma", "alpha_star", "alpha", "log_beta", "phi", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g")
sim_reg_parameters_type <- c("logical", "numeric", "numeric", "numeric", "list", "numeric", "numeric", "numeric", "numeric")
sim_reg_parameters_gamma <- c(NA, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
sim_reg_likelihood_traces <- paste0(c(paste0("total", c(paste0("_gamma",0:1), "")), "y", sim_reg_parameters), "_likelihood")
sim_reg_parameter_traces <- c("phi_vector","s_x","s_phi","s",sim_reg_parameters)
sim_reg_parameter_acceptance_traces <- paste0(sim_reg_parameters, "_accept")
sim_reg_all_traces <- unique(c(sim_reg_parameter_traces, sim_reg_parameter_acceptance_traces, sim_reg_likelihood_traces, "iteration"))

#' Get final sample from \code{sim_reg_samples} object
#'
#' @template samples
final_sample <- function(samples) {
	c(lapply(X=samples[sim_reg_parameters[-which(sim_reg_parameters == "phi")]], function(x) x[length(x)]), list(phi=samples$phi_vector[nrow(samples$phi_vector),]))
}

#' Get boolean vectors indicating whether proposed values in \code{sim_reg} Markov chain were accepted
#'
#' @template samples
#' @importFrom stats setNames
acceptance_traces <- function(samples) c(
	list(phi_accept=sapply(1:(nrow(samples$phi_vector)-1), function(x) !identical(samples$phi_vector[x,], samples$phi_vector[x+1,]))),
	lapply(
		setNames(sim_reg_parameters[sim_reg_parameters_type=="numeric"], nm=paste0(sim_reg_parameters[sim_reg_parameters_type=="numeric"], "_accept")),
		function(x) samples[[x]][1:(length(samples[[x]])-1)] != samples[[x]][2:length(samples[[x]])]
	)
)

#' @importFrom ontologyIndex get_term_info_content get_term_descendancy_matrix get_ancestors
#' @importFrom ontologySimilarity get_term_sim_mat
#' @importFrom stats sd runif rnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
tune_proposal_sds <- function(
	tune_parameters=Filter(f=Negate(is.na), x=sim_reg_parameters[sim_reg_parameters_type == "numeric" & sim_reg_parameters_gamma]),
	verbose=FALSE,
	ontology=NULL,
	y,
	x=NULL, 
	g=rep(0, length(y)),

	information_content=get_term_info_content(ontology, term_sets=x),

	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	term_sim_mat=prune_sim_mat(ontology, get_term_sim_mat(ontology, information_content, term_descendancy_matrix=term_descendancy_matrix)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1,

	gamma=(runif(1) < gamma_prior_prob)[1],
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(term_descendancy_matrix), size=3, replace=TRUE)-1,
	logit_mean_f=rnorm(n=1, mean=logit_mean_f_mean, sd=logit_mean_f_sd),
	log_alpha_plus_beta_f=rnorm(n=1, mean=log_alpha_plus_beta_f_mean, sd=log_alpha_plus_beta_f_sd),
	logit_mean_g=rnorm(n=1, mean=logit_mean_g_mean, sd=logit_mean_g_sd),
	log_alpha_plus_beta_g=rnorm(n=1, mean=log_alpha_plus_beta_g_mean, sd=log_alpha_plus_beta_g_sd),

	gamma_prior_prob=1,
	alpha_star_mean=0,
	alpha_mean=0,
	alpha_star_sd=5,
	alpha_sd=5,
	log_beta_mean=2,
	log_beta_sd=1,
	logit_mean_f_mean=1,
	logit_mean_f_sd=1,
	log_alpha_plus_beta_f_mean=2,
	log_alpha_plus_beta_f_sd=1,
	logit_mean_g_mean=0,
	logit_mean_g_sd=1.5,
	log_alpha_plus_beta_g_mean=2,
	log_alpha_plus_beta_g_sd=1,
	pseudo_alpha_star_mean=0,
	pseudo_alpha_mean=0,
	pseudo_alpha_star_sd=5,
	pseudo_alpha_sd=5,
	pseudo_log_beta_mean=2,
	pseudo_log_beta_sd=2,
	pseudo_logit_mean_f_mean=2,
	pseudo_logit_mean_f_sd=2,
	pseudo_log_alpha_plus_beta_f_mean=2,
	pseudo_log_alpha_plus_beta_f_sd=2,
	pseudo_logit_mean_g_mean=2,
	pseudo_logit_mean_g_sd=2,
	pseudo_log_alpha_plus_beta_g_mean=2,
	pseudo_log_alpha_plus_beta_g_sd=2,
	alpha_star_proposal_sd=2,
	alpha_proposal_sd=2,
	log_beta_proposal_sd=2,
	logit_mean_f_proposal_sd=2,
	log_alpha_plus_beta_f_proposal_sd=2,
	logit_mean_g_proposal_sd=2,
	log_alpha_plus_beta_g_proposal_sd=2,
	phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	phi_num_leaves_geometric_rate=1,
	lit_sims=setNames(rep(1, ncol(term_sim_mat)), colnames(term_sim_mat))
) {
	dashes <- paste0(collapse="",rep("-",getOption("width")))

	if (verbose) cat(dashes, "\n", "Tuning proposal variances for parameters:", "\n", paste("\t", tune_parameters, "\n", collapse="", sep=""), sep="")
	
	if (verbose) cat(dashes, "\n", "Discarding initial update cycles", "\n", sep="")

	burn_trace <- .Call(
		"R_sim_reg",
		PACKAGE="SimReg",
		tune_proposal_sd_conf$pre_burn,
		1L,
		FALSE,
		FALSE,
		term_sim_mat,
		term_ids,
		case_ids,
		y,
		g,

		gamma,
		alpha_star,
		alpha,
		log_beta,
		phi,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g,

		gamma_prior_prob,
		alpha_star_mean,
		alpha_mean,
		alpha_star_sd,
		alpha_sd,
		log_beta_mean,
		log_beta_sd,
		logit_mean_f_mean,
		logit_mean_f_sd,
		log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd,
		logit_mean_g_mean,
		logit_mean_g_sd,
		log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean,
		pseudo_alpha_mean,
		pseudo_alpha_star_sd,
		pseudo_alpha_sd,
		pseudo_log_beta_mean,
		pseudo_log_beta_sd,
		pseudo_logit_mean_f_mean,
		pseudo_logit_mean_f_sd,
		pseudo_log_alpha_plus_beta_f_mean,
		pseudo_log_alpha_plus_beta_f_sd,
		pseudo_logit_mean_g_mean,
		pseudo_logit_mean_g_sd,
		pseudo_log_alpha_plus_beta_g_mean,
		pseudo_log_alpha_plus_beta_g_sd,
		pseudo_phi_marginal_prior,
		alpha_star_proposal_sd,
		alpha_proposal_sd,
		log_beta_proposal_sd,
		logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd,
		phi_jumps,

		lit_sims,
		term_descendancy_matrix,
		phi_num_leaves_geometric_rate,
		FALSE,
		FALSE,
		0,
		0,
		0,
		1
	)

	#if (verbose) print(burn_time)

	start_params <- final_sample(burn_trace)
	
	cycle <- 1
	proposal_sds <- lapply(setNames(nm=tune_parameters), function(x) get(paste0(x, "_proposal_sd")))
	acceptances <- sapply(acceptance_traces(burn_trace)[paste0(tune_parameters, "_accept")], mean)
	while (cycle <= tune_proposal_sd_conf$max_tuning_cycles & (any(acceptances < min(tune_proposal_sd_conf$target_acceptance))|any(acceptances > max(tune_proposal_sd_conf$target_acceptance)))) {
		proposal_sds <- mapply(SIMPLIFY=FALSE, FUN=function(prop_sd, acc_rate) if (acc_rate < tune_proposal_sd_conf$target_acceptance[1] | acc_rate > tune_proposal_sd_conf$target_acceptance[2]) { match.fun(if (acc_rate < tune_proposal_sd_conf$target_acceptance[1]) "/" else "*")(prop_sd, (1 + tune_proposal_sd_conf$initial_rate * tune_proposal_sd_conf$rate_decay ^ (-cycle+1))) } else { prop_sd }, proposal_sds, acceptances)

		cycle_trace <- .Call(
			"R_sim_reg",
			PACKAGE="SimReg",
			tune_proposal_sd_conf$cycle_length,
			1L,
			FALSE,
			FALSE,
			term_sim_mat,
			term_ids,
			case_ids,
			y,
			g,

			start_params$gamma,
			start_params$alpha_star,
			start_params$alpha,
			start_params$log_beta,
			start_params$phi,
			start_params$logit_mean_f,
			start_params$log_alpha_plus_beta_f,
			start_params$logit_mean_g,
			start_params$log_alpha_plus_beta_g,

			gamma_prior_prob,
			alpha_star_mean,
			alpha_mean,
			alpha_star_sd,
			alpha_sd,
			log_beta_mean,
			log_beta_sd,
			logit_mean_f_mean,
			logit_mean_f_sd,
			log_alpha_plus_beta_f_mean,
			log_alpha_plus_beta_f_sd,
			logit_mean_g_mean,
			logit_mean_g_sd,
			log_alpha_plus_beta_g_mean,
			log_alpha_plus_beta_g_sd,
			pseudo_alpha_star_mean,
			pseudo_alpha_mean,
			pseudo_alpha_star_sd,
			pseudo_alpha_sd,
			pseudo_log_beta_mean,
			pseudo_log_beta_sd,
			pseudo_logit_mean_f_mean,
			pseudo_logit_mean_f_sd,
			pseudo_log_alpha_plus_beta_f_mean,
			pseudo_log_alpha_plus_beta_f_sd,
			pseudo_logit_mean_g_mean,
			pseudo_logit_mean_g_sd,
			pseudo_log_alpha_plus_beta_g_mean,
			pseudo_log_alpha_plus_beta_g_sd,
			pseudo_phi_marginal_prior,
			if ("alpha_star" %in% tune_parameters) proposal_sds[["alpha_star"]] else alpha_star_proposal_sd,
			if ("alpha" %in% tune_parameters) proposal_sds[["alpha"]] else alpha_proposal_sd,
			if ("log_beta" %in% tune_parameters) proposal_sds[["log_beta"]] else log_beta_proposal_sd,
			if ("logit_mean_f" %in% tune_parameters) proposal_sds[["logit_mean_f"]] else logit_mean_f_proposal_sd,
			if ("log_alpha_plus_beta_f" %in% tune_parameters) proposal_sds[["log_alpha_plus_beta_f"]] else log_alpha_plus_beta_f_proposal_sd,
			if ("logit_mean_g" %in% tune_parameters) proposal_sds[["logit_mean_g"]] else logit_mean_g_proposal_sd,
			if ("log_alpha_plus_beta_g" %in% tune_parameters) proposal_sds[["log_alpha_plus_beta_g"]] else log_alpha_plus_beta_g_proposal_sd,
			phi_jumps,

			lit_sims,
			term_descendancy_matrix,
			phi_num_leaves_geometric_rate,
			FALSE,
			FALSE,
			0,
			0,
			0,
			1
		)

		start_params <- final_sample(cycle_trace)
		acceptances <- sapply(acceptance_traces(cycle_trace)[paste0(tune_parameters, "_accept")], mean)

		if (verbose) {
			cat("Cycle ", cycle, "\n", sep="")
			print(data.frame(parameter=tune_parameters, proposal_sd=simplify2array(proposal_sds), acceptance_rate=acceptances),row.names=FALSE)
			cat("\n", sep="")
		}

		cycle <- cycle + 1
	}

	if (verbose) cat(dashes, "\n", sep="")

	proposal_sds
}

#' @importFrom ontologyIndex get_term_info_content get_term_descendancy_matrix get_ancestors
#' @importFrom ontologySimilarity get_term_sim_mat
#' @importFrom stats sd runif rnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
sim_reg_no_pseudopriors <- function(
	ontology=NULL,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=10000,
	thin=1,
	burn=2000,
	record_sims=FALSE,
	record_model_likelihoods=FALSE,

	information_content=get_term_info_content(ontology, term_sets=x),

	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	term_sim_mat=prune_sim_mat(ontology, get_term_sim_mat(ontology, information_content, term_descendancy_matrix=term_descendancy_matrix)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1,

	gamma=(runif(1) < gamma_prior_prob)[1],
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(term_descendancy_matrix), size=3, replace=TRUE)-1,
	logit_mean_f=rnorm(n=1, mean=logit_mean_f_mean, sd=logit_mean_f_sd),
	log_alpha_plus_beta_f=rnorm(n=1, mean=log_alpha_plus_beta_f_mean, sd=log_alpha_plus_beta_f_sd),
	logit_mean_g=rnorm(n=1, mean=logit_mean_g_mean, sd=logit_mean_g_sd),
	log_alpha_plus_beta_g=rnorm(n=1, mean=log_alpha_plus_beta_g_mean, sd=log_alpha_plus_beta_g_sd),

	gamma_prior_prob=0.05,
	alpha_star_mean=0,
	alpha_mean=0,
	alpha_star_sd=5,
	alpha_sd=5,
	log_beta_mean=2,
	log_beta_sd=1,
	logit_mean_f_mean=1,
	logit_mean_f_sd=1,
	log_alpha_plus_beta_f_mean=2,
	log_alpha_plus_beta_f_sd=1,
	logit_mean_g_mean=0,
	logit_mean_g_sd=1.5,
	log_alpha_plus_beta_g_mean=2,
	log_alpha_plus_beta_g_sd=1,
	pseudo_alpha_star_mean=0,
	pseudo_alpha_mean=0,
	pseudo_alpha_star_sd=5,
	pseudo_alpha_sd=5,
	pseudo_log_beta_mean=2,
	pseudo_log_beta_sd=2,
	pseudo_logit_mean_f_mean=2,
	pseudo_logit_mean_f_sd=2,
	pseudo_log_alpha_plus_beta_f_mean=2,
	pseudo_log_alpha_plus_beta_f_sd=2,
	pseudo_logit_mean_g_mean=2,
	pseudo_logit_mean_g_sd=2,
	pseudo_log_alpha_plus_beta_g_mean=2,
	pseudo_log_alpha_plus_beta_g_sd=2,
	alpha_star_proposal_sd=2,
	alpha_proposal_sd=2,
	log_beta_proposal_sd=2,
	logit_mean_f_proposal_sd=2,
	log_alpha_plus_beta_f_proposal_sd=2,
	logit_mean_g_proposal_sd=2,
	log_alpha_plus_beta_g_proposal_sd=2,
	phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	phi_num_leaves_geometric_rate=1,
	lit_sims=setNames(rep(1, ncol(term_sim_mat)), colnames(term_sim_mat)),
	check_args=TRUE
) {
	if (check_args) {
		stopifnot(burn < (its/thin))
		stopifnot(length(y) > 0)
		stopifnot(class(y) == "logical")
		stopifnot(length(y) == length(y))
		stopifnot(!any(is.na(y)))
		stopifnot(!any(is.na(g)))
		if (is.null(x)) {
			stopifnot(length(term_ids) == length(case_ids))
			stopifnot(length(term_ids) > 0)
			if (length(term_ids > 0)) {
				stopifnot(!any(is.na(term_ids)))
				stopifnot(max(term_ids) < nrow(term_sim_mat))
				stopifnot(min(term_ids) >= 0L)
				stopifnot(!any(is.na(case_ids)))
				stopifnot(min(case_ids) >= 0L)
				stopifnot(max(case_ids) < length(y))
			}
		} else {
			stopifnot(class(x) == "list")
			if (sum(sapply(x, length) == 0) > 0) stop(paste("All elements of x must contain at least 1 term"))
			stopifnot(!any(sapply(x, is.null)))
			stopifnot(!all(unlist(x) %in% rownames(term_sim_mat)))
			stopifnot(length(x) == length(y))
		}
		stopifnot(identical(dim(term_descendancy_matrix), dim(term_sim_mat)))
		stopifnot(!any(is.na(as.numeric(lit_sims))))
		stopifnot(all(names(lit_sims) %in% rownames(term_descendancy_matrix)))
		stopifnot(identical(rownames(term_descendancy_matrix), colnames(term_sim_mat)))
		stopifnot(identical(rownames(term_sim_mat), colnames(term_sim_mat)))
		stopifnot(identical(rownames(term_descendancy_matrix), rownames(term_sim_mat)))
	}

	result <- .Call(
		"R_sim_reg",
		PACKAGE="SimReg",
		its,
		thin,
		record_sims,
		record_model_likelihoods,
		term_sim_mat,
		term_ids,
		case_ids,
		y,
		g,

		gamma,
		alpha_star,
		alpha,
		log_beta,
		phi,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g,

		gamma_prior_prob,
		alpha_star_mean,
		alpha_mean,
		alpha_star_sd,
		alpha_sd,
		log_beta_mean,
		log_beta_sd,
		logit_mean_f_mean,
		logit_mean_f_sd,
		log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd,
		logit_mean_g_mean,
		logit_mean_g_sd,
		log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean,
		pseudo_alpha_mean,
		pseudo_alpha_star_sd,
		pseudo_alpha_sd,
		pseudo_log_beta_mean,
		pseudo_log_beta_sd,
		pseudo_logit_mean_f_mean,
		pseudo_logit_mean_f_sd,
		pseudo_log_alpha_plus_beta_f_mean,
		pseudo_log_alpha_plus_beta_f_sd,
		pseudo_logit_mean_g_mean,
		pseudo_logit_mean_g_sd,
		pseudo_log_alpha_plus_beta_g_mean,
		pseudo_log_alpha_plus_beta_g_sd,
		pseudo_phi_marginal_prior,
		alpha_star_proposal_sd,
		alpha_proposal_sd,
		log_beta_proposal_sd,
		logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd,
		phi_jumps,

		lit_sims,
		term_descendancy_matrix,
		phi_num_leaves_geometric_rate,
		FALSE,
		FALSE,
		0,
		0,
		0,
		1
	)

	#flatten likelihood traces into main results list
	result <- c(result[-which(names(result) == "likelihoods")], result[["likelihoods"]])

	#burn initial iterations
	result <- lapply(result, function(y) if (class(y) == "matrix") { if (nrow(y) == its) y[(burn+1):(nrow(y)),,drop=FALSE] else y } else { y[(burn+1):(length(y))] })

	#calculate acceptances here, based on the burn chains and before modifying phi
	acceptances <- acceptance_traces(result)

	#make phi_vector slot using term IDs instead of integer references
	result$phi_vector <- apply(result$phi_vector, 2, function(x) colnames(term_descendancy_matrix)[x+1])

	#make phi slot a list of minimal sets
	result$phi <- mapply(SIMPLIFY=FALSE, FUN="[", split(result$phi_vector, seq(nrow(result$phi_vector))), split(as_row_leaves(term_descendancy_matrix, result$phi_vector), seq(nrow(result$phi_vector))))

	c(	
		list(
			priors=list(
				gamma_prior_prob=gamma_prior_prob,
				alpha_star_mean=alpha_star_mean,
				alpha_mean=alpha_mean,
				alpha_star_sd=alpha_star_sd,
				alpha_sd=alpha_sd,
				log_beta_mean=log_beta_mean,
				log_beta_sd=log_beta_sd,
				logit_mean_f_mean=logit_mean_f_mean,
				logit_mean_f_sd=logit_mean_f_sd,
				log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
				log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
				logit_mean_g_mean=logit_mean_g_mean,
				logit_mean_g_sd=logit_mean_g_sd,
				log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
				log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
				pseudo_alpha_star_mean=pseudo_alpha_star_mean,
				pseudo_alpha_star_sd=pseudo_alpha_star_sd,
				pseudo_alpha_mean=pseudo_alpha_mean,
				pseudo_alpha_sd=pseudo_alpha_sd,
				pseudo_log_beta_mean=pseudo_log_beta_mean,
				pseudo_log_beta_sd=pseudo_log_beta_sd,
				pseudo_logit_mean_f_mean=pseudo_logit_mean_f_mean,
				pseudo_logit_mean_f_sd=pseudo_logit_mean_f_sd,
				pseudo_log_alpha_plus_beta_f_mean=pseudo_log_alpha_plus_beta_f_mean,
				pseudo_log_alpha_plus_beta_f_sd=pseudo_log_alpha_plus_beta_f_sd,
				pseudo_logit_mean_g_mean=pseudo_logit_mean_g_mean,
				pseudo_logit_mean_g_sd=pseudo_logit_mean_g_sd,
				pseudo_log_alpha_plus_beta_g_mean=pseudo_log_alpha_plus_beta_g_mean,
				pseudo_log_alpha_plus_beta_g_sd=pseudo_log_alpha_plus_beta_g_sd,
				pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
				alpha_star_proposal_sd=alpha_star_proposal_sd,
				alpha_proposal_sd=alpha_proposal_sd,
				log_beta_proposal_sd=log_beta_proposal_sd,
				logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
				log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
				logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
				log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd
			)
		), 
		result, 
		acceptances
	)
}

#' Similarity regression
#'
#' Performns Bayesian `similarity regression' on given binary genotype \code{y} (logical vector) against ontological term sets \code{x} (list of character vectors of term IDs). This could, for example, be a \code{list} of character vectors of HPO term IDs representing case phenotypes. It returns an object of class `sim_reg_samples` which is a list of traces for the sampled parameters. The results can be summarised with `summary`. Of particular interest are the estimated mean posteriors of \code{gamma} (the model selection indicator, thus giving an estimate of the probability of an association under the model assumptions - stored in the `mean_posterior_gamma' slot in the result, i.e. \code{result$mean_posterior_gamma} (which can also be calculated \code{mean(result$gamma)}), and the characteristic ontological profile phi (which can be visualised by the functions \code{\link{phi_plot}}, \code{\link{term_pair_marginals_plot}}, and \code{\link{term_marginals}}).
#'
#' @template ontology
#' @param y Logical vector of genotypes (typically TRUE for rare genotype, FALSE for common genotype).
#' @template x
#' @param g Genotype log odds offset per individual.
#' @param its Number of update cycles to perform .
#' @param thin Factor by which to thin resultant chains of parameter samples.
#' @param record_sims Logical indicating whether to record trace of similarities.
#' @param record_model_likelihoods Record likelihood of parameters under both models.
#' @param tune_proposals Logical value determining whether to adaptively tune proposal variances for \code{sim_reg} numeric parameters.
#' @param verbose Logical value determining whether to print progress of execution.
#' @param information_content Numeric vector, named by HPO IDs, containing the information content of corresponding terms.
#' @template term_descendancy_matrix
#' @template term_sim_mat
#' @template case_ids
#' @template term_ids
#' @param return_tuning_runs Logical indicating whether to return the MCMC output of the tuning phase of the inference procedure.
#' @param tuning_its Number of update cycles to perform in the tuning phase of the inference procedure.
#' @param tuning_burn Number of update cycles to discard in tuning phase.
#' @param burn Number of update cycles to discard .
#' @param tune_phi_pseudoprior Logical value determining whether tuned pseudoprior for phi is used in main Markov chain.
#' @param gamma Initial value of model selection indicator gamma..
#' @param alpha_star  Initial value of alpha_star, the rate of observing the rare genotype y = 1 under gamma = 0, i.e. the no association model .
#' @param alpha Initial value of alpha, the background rate of observing the rare genotype under gamma = 1.
#' @param log_beta Initial value of log_beta, the log of the effect size of onotological similarity.
#' @param phi Character vector of HPO term IDs giving the initial value of phi, the characteristic phenotype.
#' @param logit_mean_f Initial value of logit_mean_f.
#' @param log_alpha_plus_beta_f Initial value of log_alpha_plus_beta_f.
#' @param logit_mean_g Initial value of logit_mean_g.
#' @param log_alpha_plus_beta_g Initial value of log_alpha_plus_beta_g.
#' @param gamma_prior_prob Prior probability of gamma = 1.
#' @param alpha_star_mean Prior mean of alpha_star given gamma = 0.
#' @param alpha_mean Prior mean of alpha given gamma = 1.
#' @param alpha_star_sd Prior sd of alpha_star given gamma = 0.
#' @param alpha_sd Prior sd of alpha given gamma = 1.
#' @param log_beta_mean Prior mean of log_beta given gamma = 1.
#' @param log_beta_sd Prior sd of log_beta given gamma = 1.
#' @param logit_mean_f_mean Prior mean of logit_mean_f given gamma = 1.
#' @param logit_mean_f_sd Prior sd of logit_mean_f given gamma = 1.
#' @param log_alpha_plus_beta_f_mean Prior mean of log_alpha_plus_beta_f given gamma = 1.
#' @param log_alpha_plus_beta_f_sd Prior sd of log_alpha_plus_beta_f given gamma = 1.
#' @param logit_mean_g_mean Prior mean of logit_mean_g given gamma = 1.
#' @param logit_mean_g_sd Prior sd of logit_mean_g given gamma = 1.
#' @param log_alpha_plus_beta_g_mean Prior mean of log_alpha_plus_beta_g given gamma = 1.
#' @param log_alpha_plus_beta_g_sd Prior sd of log_alpha_plus_beta_g given gamma = 1.
#' @param alpha_star_proposal_sd Proposal sd of local jumps in MH updates of alpha_star used during inference.
#' @param alpha_proposal_sd Proposal sd of local jumps in MH updates of alpha used during inference.
#' @param log_beta_proposal_sd Proposal sd of local jumps in MH updates of log_beta used during inference.
#' @param logit_mean_f_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_f used during inference.
#' @param log_alpha_plus_beta_f_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_f used during inference.
#' @param logit_mean_g_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_g used during inference.
#' @param log_alpha_plus_beta_g_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_g used during inference.
#' @param phi_jumps Vector of HPO term IDs to be used as jumping distribution for proposal replacements of terms in phi during inference given gamma = 1.
#' @param pseudo_phi_marginal_prior Vector of HPO term IDs to be used as prior distribution on marginal probability of single term in phi given gamma = 0.
#' @param phi_num_leaves_geometric_rate Geometric parameter for truncated geometric distribution on number of leaf terms in phi.
#' @param lit_sims Numeric vector of similarities (greater than 0) of literature phenotype to individual terms (named by term ID).
#' @param favour_gamma1_factor Value by which to multiply odds of \code{P(gamma=1)/P(gamma=0)} in order to encourage better mixing and higher accuracy for a given number of iterations. Defaults to 1.
#' @param check_args Logical value determining whether arguments are checked for consistency.
#' @return List (by parameter) of vectors of consecutive parameter samples from MCMC inference.
#' @examples
#' \dontrun{
#' set.seed(0)
#' data(hpo)
#' disease_terms <- c("HP:0005537", "HP:0000729", "HP:0001873")
#' all_terms <- get_ancestors(hpo, 
#'	c(disease_terms, sample(hpo$id, size=50)))
#' y <- c(rep(FALSE, 96), rep(TRUE, 3))
#' x <- lapply(y, function(.y) minimal_set(
#'	hpo, if (!.y) sample(all_terms, size=3) else 
#'		c(sample(all_terms, size=1), disease_terms[runif(n=3) < 0.8])))
#' sim_reg_out <- sim_reg(ontology=hpo, x=x, y=y)
#' mean(sim_reg_out$gamma)
#' phi_plot(hpo, 
#'	sim_reg_out$phi[sim_reg_out$gamma])
#' }
#' @export
#' @importFrom stats rnorm runif
#' @importFrom ontologyIndex get_term_descendancy_matrix
sim_reg <- function(
	ontology=NULL,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=20000,
	thin=1,
	record_sims=FALSE,
	record_model_likelihoods=FALSE,

	tune_proposals=TRUE, 
	verbose=FALSE,

	information_content=get_term_info_content(ontology, term_sets=x),
	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	term_sim_mat=prune_sim_mat(ontology, get_term_sim_mat(ontology, information_content, term_descendancy_matrix=term_descendancy_matrix)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1,

	return_tuning_runs=FALSE,
	tuning_its=its,
	tuning_burn=as.integer(tuning_its/5),
	burn=as.integer(its/5),
	tune_phi_pseudoprior=TRUE,

	gamma=(runif(1) < gamma_prior_prob),
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(term_descendancy_matrix), size=3, replace=TRUE)-1,
	logit_mean_f=rnorm(n=1, mean=logit_mean_f_mean, sd=logit_mean_f_sd),
	log_alpha_plus_beta_f=rnorm(n=1, mean=log_alpha_plus_beta_f_mean, sd=log_alpha_plus_beta_f_sd),
	logit_mean_g=rnorm(n=1, mean=logit_mean_g_mean, sd=logit_mean_g_sd),
	log_alpha_plus_beta_g=rnorm(n=1, mean=log_alpha_plus_beta_g_mean, sd=log_alpha_plus_beta_g_sd),

	gamma_prior_prob=0.05,
	alpha_star_mean=0,
	alpha_mean=0,
	alpha_star_sd=5,
	alpha_sd=5,
	log_beta_mean=2,
	log_beta_sd=1,
	logit_mean_f_mean=1,
	logit_mean_f_sd=1,
	log_alpha_plus_beta_f_mean=2,
	log_alpha_plus_beta_f_sd=1,
	logit_mean_g_mean=0,
	logit_mean_g_sd=1.5,
	log_alpha_plus_beta_g_mean=2,
	log_alpha_plus_beta_g_sd=1,
	alpha_star_proposal_sd=2,
	alpha_proposal_sd=2,
	log_beta_proposal_sd=2,
	logit_mean_f_proposal_sd=2,
	log_alpha_plus_beta_f_proposal_sd=2,
	logit_mean_g_proposal_sd=2,
	log_alpha_plus_beta_g_proposal_sd=2,

	phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	pseudo_phi_marginal_prior=phi_jumps,

	phi_num_leaves_geometric_rate=1,
	lit_sims=setNames(rep(1, ncol(term_sim_mat)), colnames(term_sim_mat)),
	favour_gamma1_factor=1,
	check_args=TRUE
) {
	if (is.null(lit_sims)) { lit_sims <- rep(1, ncol(term_sim_mat)) }
	else if (is.null(names(lit_sims))) { if (length(lit_sims) != ncol(term_sim_mat)) stop("'lit_sims' parameter must have length = ncol(term_sim_mat)") }
	else { if (!identical(names(lit_sims), colnames(term_sim_mat))) lit_sims <- ifelse(colnames(term_sim_mat) %in% names(lit_sims), lit_sims[colnames(term_sim_mat)], 1) }

	if (is.character(term_ids)) 
		term_ids <- match(term_ids, colnames(term_sim_mat))-1

	if (is.character(phi_jumps)) 
		phi_jumps <- match(phi_jumps, colnames(term_sim_mat))-1

	if (is.character(pseudo_phi_marginal_prior)) 
		pseudo_phi_marginal_prior <- match(pseudo_phi_marginal_prior, colnames(term_sim_mat))-1

	proposal_sds <- if (!tune_proposals) { lapply(setNames(nm=sim_reg_parameters[sim_reg_parameters_type == "numeric"]), function(x) get(paste0(x, "_proposal_sd"))) } else {
		g0_proposal_sds <- tune_proposal_sds(
			tune_parameters="alpha_star",
			verbose=verbose,
			ontology=ontology,
			y=y,
			x=x, 
			g=g,

			information_content=information_content,

			term_descendancy_matrix=term_descendancy_matrix,
			term_sim_mat=term_sim_mat,
			case_ids=case_ids,
			term_ids=term_ids,

			gamma=FALSE,
			alpha_star=alpha_star,
			alpha=alpha,
			log_beta=log_beta,
			phi=phi,
			logit_mean_f=logit_mean_f,
			log_alpha_plus_beta_f=log_alpha_plus_beta_f,
			logit_mean_g=logit_mean_g,
			log_alpha_plus_beta_g=log_alpha_plus_beta_g,

			gamma_prior_prob=0,
			alpha_star_mean=alpha_star_mean,
			alpha_mean=alpha_mean,
			alpha_star_sd=alpha_star_sd,
			alpha_sd=alpha_sd,
			log_beta_mean=log_beta_mean,
			log_beta_sd=log_beta_sd,
			logit_mean_f_mean=logit_mean_f_mean,
			logit_mean_f_sd=logit_mean_f_sd,
			log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
			log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
			logit_mean_g_mean=logit_mean_g_mean,
			logit_mean_g_sd=logit_mean_g_sd,
			log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
			log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
			alpha_star_proposal_sd=alpha_star_proposal_sd,
			alpha_proposal_sd=alpha_proposal_sd,
			log_beta_proposal_sd=log_beta_proposal_sd,
			logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
			log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
			logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
			log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd,
			phi_jumps=phi_jumps,
			pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
			phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
			lit_sims=lit_sims
		)
		
		g1_proposal_sds <- tune_proposal_sds(
			tune_parameters=Filter(f=Negate(is.na), x=sim_reg_parameters[sim_reg_parameters_type == "numeric" & sim_reg_parameters_gamma]),
			verbose=verbose,
			ontology=ontology,
			y=y,
			x=x, 
			g=g,

			information_content=information_content,

			term_descendancy_matrix=term_descendancy_matrix,
			term_sim_mat=term_sim_mat,
			case_ids=case_ids,
			term_ids=term_ids,

			gamma=TRUE,
			alpha_star=alpha_star,
			alpha=alpha,
			log_beta=log_beta,
			phi=phi,
			logit_mean_f=logit_mean_f,
			log_alpha_plus_beta_f=log_alpha_plus_beta_f,
			logit_mean_g=logit_mean_g,
			log_alpha_plus_beta_g=log_alpha_plus_beta_g,

			gamma_prior_prob=1,
			alpha_star_mean=alpha_star_mean,
			alpha_mean=alpha_mean,
			alpha_star_sd=alpha_star_sd,
			alpha_sd=alpha_sd,
			log_beta_mean=log_beta_mean,
			log_beta_sd=log_beta_sd,
			logit_mean_f_mean=logit_mean_f_mean,
			logit_mean_f_sd=logit_mean_f_sd,
			log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
			log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
			logit_mean_g_mean=logit_mean_g_mean,
			logit_mean_g_sd=logit_mean_g_sd,
			log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
			log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
			alpha_star_proposal_sd=alpha_star_proposal_sd,
			alpha_proposal_sd=alpha_proposal_sd,
			log_beta_proposal_sd=log_beta_proposal_sd,
			logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
			log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
			logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
			log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd,
			phi_jumps=phi_jumps,
			pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
			phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
			lit_sims=lit_sims
		)

		c(g1_proposal_sds, g0_proposal_sds)
	}

	dashes <- paste0(collapse="",rep("-",getOption("width")))
	if (verbose) cat(dashes, "\n", "Tuning pseudoprior distributions for gamma=0 model", "\n", sep="")
	
	null.time <- system.time(null.out <- sim_reg_no_pseudopriors(
		ontology=ontology,
		term_descendancy_matrix=term_descendancy_matrix,
		term_sim_mat=term_sim_mat,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=tuning_its,
		thin=thin,
		record_sims=record_sims,
		record_model_likelihoods=record_model_likelihoods,
		burn=tuning_burn,

		gamma=FALSE,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=0,
		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		alpha_star_proposal_sd=proposal_sds[["alpha_star"]],
		alpha_proposal_sd=proposal_sds[["alpha"]],
		log_beta_proposal_sd=proposal_sds[["log_beta"]],
		logit_mean_f_proposal_sd=proposal_sds[["logit_mean_f"]],
		log_alpha_plus_beta_f_proposal_sd=proposal_sds[["log_alpha_plus_beta_f"]],
		logit_mean_g_proposal_sd=proposal_sds[["logit_mean_g"]],
		log_alpha_plus_beta_g_proposal_sd=proposal_sds[["log_alpha_plus_beta_g"]],
		phi_jumps=phi_jumps,
		pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
		
		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims,
		check_args=check_args
	))

	if (verbose) print(null.time)

	if (verbose) cat(dashes, "\n", "Tuning pseudoprior distributions for gamma=1 model", "\n", sep="")
	pheno.time <- system.time(pheno.out <- sim_reg_no_pseudopriors(
		ontology=ontology,
		term_descendancy_matrix=term_descendancy_matrix,
		term_sim_mat=term_sim_mat,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=tuning_its,
		thin=thin,
		record_sims=record_sims,
		record_model_likelihoods=record_model_likelihoods,
		burn=tuning_burn,

		gamma=TRUE,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=1,
		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		alpha_star_proposal_sd=proposal_sds[["alpha_star"]],
		alpha_proposal_sd=proposal_sds[["alpha"]],
		log_beta_proposal_sd=proposal_sds[["log_beta"]],
		logit_mean_f_proposal_sd=proposal_sds[["logit_mean_f"]],
		log_alpha_plus_beta_f_proposal_sd=proposal_sds[["log_alpha_plus_beta_f"]],
		logit_mean_g_proposal_sd=proposal_sds[["logit_mean_g"]],
		log_alpha_plus_beta_g_proposal_sd=proposal_sds[["log_alpha_plus_beta_g"]],
		phi_jumps=phi_jumps,
		pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims,
		check_args=check_args
	))
	
	if (verbose) print(pheno.time)

	g_prior <- gamma_prior_prob * favour_gamma1_factor / (1-gamma_prior_prob+favour_gamma1_factor * gamma_prior_prob)

	if (verbose) cat(dashes, "\n", "Sampling parameters from Markov chain", "\n", sep="")
	main.time <- system.time(result <- sim_reg_no_pseudopriors(
		ontology=ontology,
		term_descendancy_matrix=term_descendancy_matrix,
		term_sim_mat=term_sim_mat,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=its,
		thin=thin,
		record_sims=record_sims,
		record_model_likelihoods=record_model_likelihoods,
		burn=burn,

		gamma=gamma,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=g_prior,

		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean=mean(null.out$alpha_star[1:length(null.out$alpha_star)]),
		pseudo_alpha_star_sd=sd(null.out$alpha_star[1:length(null.out$alpha_star)]),
		pseudo_alpha_mean=mean(pheno.out$alpha[1:length(pheno.out$alpha)]),
		pseudo_alpha_sd=sd(pheno.out$alpha[1:length(pheno.out$alpha)]),
		pseudo_log_beta_mean=mean(pheno.out$log_beta[1:length(pheno.out$log_beta)]),
		pseudo_log_beta_sd=sd(pheno.out$log_beta[1:length(pheno.out$log_beta)]),
		pseudo_logit_mean_f_mean=mean(pheno.out$logit_mean_f[1:length(pheno.out$logit_mean_f)]),
		pseudo_logit_mean_f_sd=sd(pheno.out$logit_mean_f[1:length(pheno.out$logit_mean_f)]),
		pseudo_log_alpha_plus_beta_f_mean=mean(pheno.out$log_alpha_plus_beta_f[1:length(pheno.out$log_alpha_plus_beta_f)]),
		pseudo_log_alpha_plus_beta_f_sd=sd(pheno.out$log_alpha_plus_beta_f[1:length(pheno.out$log_alpha_plus_beta_f)]),
		pseudo_logit_mean_g_mean=mean(pheno.out$logit_mean_g[1:length(pheno.out$logit_mean_g)]),
		pseudo_logit_mean_g_sd=sd(pheno.out$logit_mean_g[1:length(pheno.out$logit_mean_g)]),
		pseudo_log_alpha_plus_beta_g_mean=mean(pheno.out$log_alpha_plus_beta_g[1:length(pheno.out$log_alpha_plus_beta_g)]),
		pseudo_log_alpha_plus_beta_g_sd=sd(pheno.out$log_alpha_plus_beta_g[1:length(pheno.out$log_alpha_plus_beta_g)]),
		alpha_star_proposal_sd=proposal_sds[["alpha_star"]],
		alpha_proposal_sd=proposal_sds[["alpha"]],
		log_beta_proposal_sd=proposal_sds[["log_beta"]],
		logit_mean_f_proposal_sd=proposal_sds[["logit_mean_f"]],
		log_alpha_plus_beta_f_proposal_sd=proposal_sds[["log_alpha_plus_beta_f"]],
		logit_mean_g_proposal_sd=proposal_sds[["logit_mean_g"]],
		log_alpha_plus_beta_g_proposal_sd=proposal_sds[["log_alpha_plus_beta_g"]],

		#pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
		#phi_jumps=phi_jumps,
		pseudo_phi_marginal_prior=if (tune_phi_pseudoprior) c(0:(ncol(term_sim_mat)-1), match(as.character(pheno.out[["phi_vector"]]), colnames(term_sim_mat))-1) else pseudo_phi_marginal_prior,
		phi_jumps=phi_jumps,#c(0:(ncol(term_sim_mat)-1), match(as.character(pheno.out[["phi_vector"]]), colnames(term_sim_mat))-1),

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims,
		check_args=check_args
	))

	if (verbose) print(main.time)

	if (return_tuning_runs) {
		result$tune_gamma0 <- null.out
		result$tune_gamma1 <- pheno.out
	}

	bf <- (1-g_prior) * mean(result$gamma) / (1-mean(result$gamma)) / g_prior
	result$mean_posterior_gamma <- if (mean(result$gamma) == 1) { 1 } else { if (mean(result$gamma) == 0) 0 else bf * gamma_prior_prob / (bf * gamma_prior_prob + 1 - gamma_prior_prob) }
	result$proposal_sds <- proposal_sds

	if (verbose) cat(dashes, "\n", "Done. Mean posterior gamma = ", result$mean_posterior_gamma, "\n", dashes, "\n", sep="")
	class(result) <- "sim_reg_samples"
	result
}


