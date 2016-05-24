## -- simulate() *method* {since Sept.2005}  (instead of function):
## -- following simulate.lm() in .../R/src/library/stats/R/lm.R
simulate.vlmc <-
    function(object, nsim = 1, seed = NULL,
             n, n.start = 64 * object$size[["context"]],
             integer.return = FALSE, keep.RSeed = TRUE, ...)
{
    ## Author: Martin Maechler, Date: 10 Apr 2000
    if(!is.vlmc(object))
	stop("first argument must be a \"vlmc\" object; see ?vlmc")
    if(sys.nframe() < 2 || !identical(sys.call(-1)[[1]], quote(simulate)))
	warning("Calling simulate.vlmc(<vlmc_obj>, ..) explicitly is deprecated;\n",
	      "  Use     simulate     (<vlmc_obj>, ..) instead!")
    if(length(list(...))) warning("ignoring extraneous arguments")
    ivlmc <- object $ vlmc

    ## back compatibility: till Oct.2014, only 'n = *' worked correctly
    if(!missing(n)) {
	if(missing(nsim)) {
	    warning("The argument 'n' is deprecated; use 'nsim' instead.")
	    nsim <- n
	} else
	    stop("Both 'n' and 'nsim' have been specified; use 'nsim' only, now")
    }
    if(0 > (nsim <- as.integer(nsim)))
        stop("required output-length nsim must be >= 0")
    n <- nsim + (n.start <- as.integer(n.start))
    if(keep.RSeed) { # behave as generic '?simulate' says
        if(!exists(".Random.seed", envir = .GlobalEnv))
            runif(1)                 # initialize the RNG if necessary
        if(is.null(seed))
            RNGstate <- .Random.seed
        else {
            R.seed <- .Random.seed
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
    } ## else: back-compatible to older 'simulate.vlmc'

    ## FIXME: also pass n.start -- and drop the first n.start in C code
    iy <- .Call(vlmc_sim, ivlmc, n)
    if(n.start) iy <- iy[-seq_len(n.start)]
    if(length(name <- deparse(substitute(object))) > 1)
	name <- paste(name[1], "...")
    structure(if(integer.return) iy else strsplit(object$alpha, NULL)[[1]][1 + iy],
              name=name, seed = if(keep.RSeed) RNGstate,
              class = "simulate.vlmc")
}

print.simulate.vlmc <- function (x, ...) {
    cat("simulate(", attr(x,"name"), ", ", length(x), ", ..):\n", sep="")
    print(as.vector(x), ...)
    if(!is.null(seed <- attr(x, "seed"))) { cat(" seed:"); str(seed, vec.len=2) }
    invisible(x)
}
