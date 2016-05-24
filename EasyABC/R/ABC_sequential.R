## FUNCTION ABC_sequential: Sequential ABC methods (Beaumont et al. 2009, Drovandi
## & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
ABC_sequential <- function(method, model, prior, nb_simul, summary_stat_target, prior_test = NULL, 
    n_cluster = 1, use_seed = FALSE, verbose = FALSE, dist_weights=NULL, ...) {
    ## checking errors in the inputs
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior)) 
        stop("'prior' is missing")
    data = .wrap_constants_in_model(prior, model, use_seed)
    prior = data$new_prior
    model = data$new_model
    prior = .process_prior(prior)
    if (!is.null(prior_test)) 
        .check_prior_test(length(prior), prior_test)
    if (missing(nb_simul)) 
        stop("'nb_simul' is missing")
    if (missing(summary_stat_target)) 
        stop("'summary_stat_target' is missing")
    if (!any(method == c("Beaumont", "Drovandi", "Delmoral", "Lenormand", "Emulation"))) {
        stop("Method must be Beaumont, Drovandi, Delmoral, Lenormand or Emulation")
    }
    if (!is.vector(nb_simul)) 
        stop("'nb_simul' has to be a number.")
    if (length(nb_simul) > 1) 
        stop("'nb_simul' has to be a number.")
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1.")
    nb_simul = floor(nb_simul)
    if (!is.vector(summary_stat_target)) 
        stop("'summary_stat_target' has to be a vector.")
    if (!is.vector(n_cluster)) 
        stop("'n_cluster' has to be a number.")
    if (length(n_cluster) > 1) 
        stop("'n_cluster' has to be a number.")
    if (n_cluster < 1) 
        stop("'n_cluster' has to be a positive number.")
    n_cluster = floor(n_cluster)
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    if (!is.null(dist_weights) && length(dist_weights)!=length(summary_stat_target)) {
        stop("'dist_weights' has to be the same length than 'summary_stat_target'")
    }
    sequential = NULL
    if (n_cluster == 1) {
        sequential = .ABC_sequential(method, model, prior, prior_test, nb_simul, 
            summary_stat_target, use_seed, verbose, dist_weights=dist_weights, ...)
    } else {
        if (method=="Emulation") {
            stop("'Emulation' method isn't yet available in 'cluster' mode'")
        }
        if (use_seed == FALSE) {
            stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
        }
        sequential = .ABC_sequential_cluster(method, model, prior, prior_test, nb_simul, 
            summary_stat_target, n_cluster, use_seed, verbose, dist_weights=dist_weights, ...)
    }
    sequential
} 
