## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
ABC_emulation <- function(model, prior, nb_design_pts, nb_simul, prior_test = NULL, summary_stat_target = NULL, 
    emulator_span = 50, tol = NULL, use_seed = FALSE, seed_count = 0, n_cluster = 1, verbose = FALSE, 
    progress_bar = FALSE) {
    ## checking errors in the inputs
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior)) 
        stop("'prior' is missing")
    if (!is.null(prior_test)) 
        .check_prior_test(length(prior), prior_test)
    data = .wrap_constants_in_model(prior, model, use_seed)
    prior = data$new_prior
    model = data$new_model
    prior = .process_prior(prior)
    if (missing(nb_simul)) 
        stop("'nb_simul' is missing")
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1")
    if (missing(nb_design_pts)) 
        stop("'nb_design_pts' is missing")
    if (nb_design_pts > nb_simul) 
        stop("'nb_design_pts' must be a number less than 'nb_simul'")
    if ((!is.null(summary_stat_target)) && (!is.vector(summary_stat_target))) {
        stop("'summary_stat_target' has to be a number")
    }
    if ((!is.null(tol)) && (!is.vector(tol))) {
        stop("'tol' has to be a number")
    }
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.vector(seed_count)) 
        stop("'seed_count' has to be a number")
    if (length(seed_count) > 1) 
        stop("'seed_count' has to be a number")
    if (seed_count < 0) 
        stop("'seed_count' has to be a positive number")
    if (!is.vector(n_cluster)) 
        stop("'n_cluster' has to be a number.")
    if (length(n_cluster) > 1) 
        stop("'n_cluster' has to be a number.")
    if (n_cluster < 1) 
        stop("'n_cluster' has to be a positive number.")
    n_cluster = floor(n_cluster)
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    if (!is.logical(progress_bar)) 
        stop("'progress_bar' has to be boolean")
    nb_simul = floor(nb_simul)
    seed_count = floor(seed_count)
    if ((!is.null(summary_stat_target)) && (is.null(tol))) {
        stop("'tol' is missing")
    }
    span = min(1, emulator_span/nb_design_pts)

    if (n_cluster != 1) {
        stop("'Emulation' method isn't yet available in 'cluster' mode'")
    } else {
        rejection1 = .ABC_rejection(model, prior, prior_test, nb_design_pts, use_seed, 
            seed_count, verbose, progress_bar)
        
        if (use_seed) {
            model_emulator = function(parameters) {
                .emulator_locreg_deg2(parameters[2:length(parameters)], rejection1$param, rejection1$stats, span)
            }
        } else {
            model_emulator = function(parameters) {
                .emulator_locreg_deg2(parameters, rejection1$param, rejection1$stats, span)
            }
        }
        rejection = .ABC_rejection(model_emulator, prior, prior_test, nb_simul, use_seed, 
            seed_count + nb_design_pts, verbose, progress_bar)
        res = NULL
        if (is.null(summary_stat_target)) {
            res = list(param = rejection$param, stats = rejection$stats, weights = rejection$weights, 
                stats_normalization = rejection$stats_normalization, nsim = rejection$nsim, 
                computime = rejection$computime)
        } else {
            options(warn = -1)
            rej = abc(summary_stat_target, rejection$param, rejection$stats, tol, method = "rejection")
            options(warn = 0)
            nr = dim(rej$unadj.values)[1]
            res = list(param = rej$unadj.values, stats = rej$ss, weights = array(1/nr, 
                nr), stats_normalization = rejection$stats_normalization, nsim = rejection$nsim, 
                nrec = nr, computime = rejection$computime)
        }
        res
    }
} 
