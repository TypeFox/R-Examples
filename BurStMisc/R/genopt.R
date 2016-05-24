"genopt" <-
function (fun, population, lower = -Inf, upper = Inf, scale = dcontrol["eps"], 
    add.args = NULL, control = genopt.control(...), ...) 
{
    if(!exists(".Random.seed")) runif(1)
    random.seed <- .Random.seed
    if (is.character(fun)) 
        fun <- get(fun, mode = "function")
    fun.args <- c(list(NULL), add.args)
    go.rectify <- function(pars, lower, upper) {
        pars[pars < lower] <- lower[pars < lower]
        pars[pars > upper] <- upper[pars > upper]
        pars
    }
    if (is.list(population)) {
        objective <- population$objective
        funevals <- population$funevals
        population <- population$population
        popsize <- ncol(population)
        if (is.null(popsize) || length(objective) != popsize) 
            stop("bad input population")
        if (!is.numeric(funevals) || is.na(funevals)) {
            funevals <- 0
            warning("funevals starting at 0")
        }
    }
    else {
        if (!is.matrix(population)) 
            stop("bad input population")
        popsize <- ncol(population)
        objective <- numeric(popsize)
        npar <- nrow(population)
        lower <- rep(lower, length = npar)
        upper <- rep(upper, length = npar)
        if (any(upper < lower)) 
            stop("upper element smaller than lower")
        for (i in 1:popsize) {
            population[, i] <- fun.args[[1]] <- go.rectify(population[, 
                i], lower, upper)
            objective[i] <- do.call("fun", fun.args)
        }
        funevals <- popsize
    }
    icontrol <- control$icontrol
    dcontrol <- control$dcontrol
    trace <- icontrol["trace"]
    maxeval <- icontrol["maxeval"]
    minobj <- min(objective)
    npar <- nrow(population)
    if (trace) {
        cat("objectives go from", format(minobj), "to", format(max(objective)),
            "\n")
    }
    if (icontrol["random.n"]) {
        par.range <- apply(population, 1, range)
        par.range[2, par.range[2, ] == par.range[1, ]] <- par.range[2, 
            par.range[2, ] == par.range[1, ]] + dcontrol["scale.min"]
        maxobj <- max(objective)
	randn <- min(icontrol["random.n"], maxeval - funevals)
        for (i in 1:randn) {
            fun.args[[1]] <- runif(npar, par.range[1, ], par.range[2, 
                ])
            this.obj <- do.call("fun", fun.args)
            if (this.obj < maxobj) {
                maxind <- order(objective)[popsize]
                population[, maxind] <- fun.args[[1]]
                objective[maxind] <- this.obj
                maxobj <- max(objective)
            }
        }
	funevals <- funevals + randn
        if (trace) {
            cat("objectives go from", format(minobj), "to", format(maxobj), 
                "\n")
        }
    }
    njit <- icontrol["jitters.n"]
    lower <- rep(lower, length = npar)
    upper <- rep(upper, length = npar)
    if (any(upper < lower)) 
        stop("upper element smaller than lower")
    scale[scale < dcontrol["scale.min"]] <- dcontrol["scale.min"]
    scale <- rep(scale, length = npar)
    prob <- dcontrol["prob"]
    prob <- c(prob, 1 - prob)
    for (i in 1:icontrol["births"]) {
        if (funevals >= maxeval) 
            break
        parents <- sample(popsize, 2)
        child <- population[, parents[1]]
        cloc <- sample(c(TRUE, FALSE), npar, replace = TRUE, prob = prob)
        if (all(cloc)) 
            cloc[sample(npar, 1)] <- FALSE
        else if (all(!cloc)) 
            cloc[sample(npar, 1)] <- TRUE
        child[cloc] <- population[cloc, parents[2]]
        fun.args[[1]] <- child
        child.obj <- do.call("fun", fun.args)
        funevals <- funevals + 1
        parent.obj <- objective[parents]
        survive <- child.obj < max(parent.obj)
        if (trace) {
            cat(i, "parents:", parent.obj, "child:", format(child.obj), 
                if (survive) 
                  "(improve)", "\n")
        }
        if (survive || (child.obj == parent.obj[1] && child.obj == 
            parent.obj[2])) {
            if (parent.obj[1] > parent.obj[2]) 
                out <- parents[1]
            else out <- parents[2]
            population[, out] <- child
            objective[out] <- child.obj
            if (trace && child.obj < minobj) {
                minobj <- child.obj
                cat("new minimum\n")
            }
            for (i in seq(length = njit)) {
                fun.args[[1]] <- jchild <- go.rectify(rnorm(npar, 
                  child, scale), lower, upper)
                jchild.obj <- do.call("fun", fun.args)
                funevals <- funevals + 1
                if (jchild.obj < child.obj) {
                  child <- population[, out] <- jchild
                  child.obj <- objective[out] <- jchild.obj
                  if (trace) {
                    cat("jitter successsful:", format(jchild.obj), 
                      "\n")
                    if (jchild.obj < minobj) {
                      cat("new minimum\n")
                      minobj <- jchild.obj
                    }
                  }
                }
            }
        }
    }
    ord <- order(objective)
    list(population = population[, ord], objective = objective[ord], 
        funevals = funevals, random.seed = random.seed, call = match.call())
}
