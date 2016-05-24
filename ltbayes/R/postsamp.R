postsamp <-
function(fmodel, y, method = "metrop", control = list(), ...) 
{
		if (!is.vector(y) & !is.matrix(y))
				stop("responses must be a vector or matrix")	
		if (!(method %in% c("metrop","adapt")))
			stop("only 'metrop' and 'adapt' methods are available")
		if (method == "metrop") {
			control.mcmc <- list(initial = 0, nbatch = 1000, blen = 1, nspac = 1, scale = 1)
			control.mcmc[names(control)] <- control
			f <- function(z, ...) fmodel(z, ...)$post
			out <- metrop(f, initial = control.mcmc$initial,
				nbatch = control.mcmc$nbatch, blen = control.mcmc$blen, 
				nspac = control.mcmc$nspac, scale = control.mcmc$scale, y = y, ...)
		}
		if (method == "adapt") {
			control.mcmc <- list(pars = 0, prop_sigma = NULL, par_names = NULL, 
				iterations = 1000, burn_in = 0, adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE)
			control.mcmc[names(control)] <- control
			f <- function(z, ...) fmodel(z, ...)$post
			out <- Metro_Hastings(f, pars = control.mcmc$pars, prop_sigma = control.mcmc$prop_sigma,
				par_names = control.mcmc$par_names, iterations = control.mcmc$iterations, 
				burn_in = control.mcmc$burn_in, adapt_par = control.mcmc$adapt_par, 
				quiet = control.mcmc$quiet, y = y, ...)
		}
    return(out)
}
