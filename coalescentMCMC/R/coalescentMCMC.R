## coalescentMCMC.R (2013-12-13)

##   Run MCMC for Coalescent Trees

## Copyright 2012-2013 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

coalescentMCMC <-
    function(x, ntrees = 3000, burnin = 1000, frequency = 1,
             tree0 = NULL, model = NULL, printevery = 100)
{
    if (packageVersion("phangorn") >= "1.99.5") {
        edQt <- phangorn::edQt
        lli <- phangorn::lli
        pml.fit <- phangorn::pml.fit
        pml.free <- phangorn::pml.free
        pml.init <- phangorn::pml.init
    }

    on.exit({
        if (packageVersion("phangorn") >= "1.99.5") pml.free()
        if (k < nOut) {
            if (!k) stop("burn-in period not yet finished")
            TREES <- TREES[seq_len(k)]
            LL <- LL[seq_len(i)]
            params <- params[seq_len(i), ]
            warning(paste("MCMC interrupted after", i, "generations"))
        }
        ## compress the list of trees:
        attr(TREES, "TipLabel") <- TREES[[1L]]$tip.label
        lapply(TREES, function(x) x$tip.label <- NULL)
        class(TREES) <- "multiPhylo"

        suffix <- 1
        list.trees <- .get.list.trees()
        if (l <- length(list.trees))
            suffix <- 1 + as.numeric(sub("TREES_", "", list.trees[l]))
        suffix <- sprintf("%03d", suffix)
        assign(paste("TREES", suffix, sep = "_"), TREES,
               envir = .coalescentMCMCenv)

        i <- i - 1L

        MCMCstats <- get("MCMCstats", envir = .coalescentMCMCenv)
        MCMCstats[[suffix]] <- c(k, burnin, frequency, i, j)
        assign("MCMCstats", MCMCstats, envir = .coalescentMCMCenv)

        LL <- cbind(LL, params)
        colnames(LL) <- c("logLik", para.nms)
        LL <- mcmc(LL, start = 1, end = i)
        return(LL)
    })

    verbose <- as.logical(printevery)

    if (is.null(tree0)) {
        d <- dist.dna(x, "JC69")
        tree0 <- as.phylo(hclust(d, "average"))
    }

    X <- phyDat(x)
    n <- length(tree0$tip.label)
    nodeMax <- 2*n - 1
    nOut <- ntrees
    nOut2 <- ntrees * frequency + burnin

    getlogLik <- function(phy, X) pml(phy, X)$logLik

    if (packageVersion("phangorn") >= "1.99.5") {
        INV <- Matrix(lli(X, tree0), sparse = TRUE)
        ll.0 <- numeric(attr(X, "nr"))
        ## by Klaus (ensures that tip labels of tree and data have same order):
        X <- subset(X, tree0$tip.label)
        ##
        bf <- rep(0.25, 4)
        eig <- edQt()
        pml.init(X)
        getlogLik <- function(phy, X) {
            phy <- reorder(phy, "postorder")
            pml.fit(phy, X, bf = bf, eig = eig, INV = INV, ll.0 = ll.0)
        }
    }

    TREES <- vector("list", nOut)
    LL <- numeric(nOut2)
    TREES[[1L]] <- tree0
    lnL0 <- getlogLik(tree0, X)
    LL[1L] <- lnL0

    if (is.null(model)) {
        np <- 1L
        para.nms <- "theta"
        ## quantities to calculate THETA:
        two2n <- 2:n
        K4theta <- length(two2n)
        tmp <- two2n * (two2n - 1) # == 2 * choose(two2n, 2)
        getparams <- function(phy, bt) {
            x4theta <- rev(diff(c(0, sort(bt))))
            sum(x4theta * tmp)/K4theta
        }
        f.theta <- function(t, p) p
    } else {
        switch(model, time = {
            np <- 2L
            para.nms <- c("theta0", "rho")
            getparams <- function(phy, bt) { # 'bt' is not used but is needed to have the same arguments than above
                halfdev <- function(p) {
                    if (any(p <= 0) || any(is.nan(p))) return(1e100)
                    -dcoal.time(phy, p[1], p[2], log = TRUE)
                }
                out <- nlminb(c(0.02, 0), halfdev)
                out$par
            }
            f.theta <- function(t, p) p[1] * exp(p[2] * t)
        }, step = {
            np <- 3L
            para.nms <- c("theta0", "theta1", "tau")
            getparams <- function(phy, bt) {
                halfdev <- function(p) {
                    if (any(p <= 0) || any(is.nan(p))) return(1e100)
                    -dcoal.step(phy, p[1], p[2], p[3], log = TRUE)
                }
                out <- nlminb(c(0.02, 0.02, bt[1]/2), halfdev)
                out$par
            }
            f.theta <- function(t, p) ifelse(t <= p[3], p[1], p[2])
        }, linear = {
            np <- 3L
            para.nms <- c("theta0", "thetaT", "TMRCA")
            getparams <- function(phy, bt) {
                halfdev <- function(p) {
                    if (any(p <= 0) || any(is.nan(p))) return(1e100)
                    -dcoal.linear(phy, p[1], p[2], p[3], log = TRUE)
                }
                out <- nlminb(c(0.02, 0.02, bt[1]), halfdev)
                out$par
            }
            f.theta <- function(t, p) p[1] + t * (p[2] - p[1])/p[3]
        })
    }
    params <- matrix(0, nOut2, np)

    i <- 2L
    j <- 0L # number of accepted moves
    k <- 0L # number of sampled trees

    if (verbose) {
        cat("Running the Markov chain:\n")
        cat("  Number of trees to output:", ntrees, "\n")
        cat("  Burn-in period:", burnin, "\n")
        cat("  Sampling frequency:", frequency, "\n")
        cat("  Number of generations to run:", ntrees * frequency + burnin, "\n")
        cat("Generation    Nb of accepted moves\n")
    }

    bt0 <- branching.times(tree0)
    params[1L, ] <- para0 <- getparams(tree0, bt0)

    nodesToSample <- (n + 2):nodeMax

    while (k < nOut) {
        if (verbose) if (! i %% printevery)
            cat("\r  ", i, "                ", j, "           ")

        ## select one internal node excluding the root:
        target <- sample(nodesToSample, 1L) # target node for rearrangement
        THETA <- f.theta(bt0[target - n], para0) # the value of THETA at this node

        tr.b <- NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)
        ## do TipInterchange() every 10 steps:
        ## tr.b <-
        ##     if (! i %% 10) TipInterchange(tree0, n)
        ##     else NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)

        if (!(i %% frequency) && i > burnin) {
            k <- k + 1L
            TREES[[k]] <- tr.b
        }
        lnL.b <- getlogLik(tr.b, X)
        LL[i] <- lnL.b
        ## calculate theta for the proposed tree:
        bt <- branching.times(tr.b)
        params[i, ] <- para <- getparams(tr.b, bt)
        i <- i + 1L
        ACCEPT <- if (is.na(lnL.b)) FALSE else {
            if (lnL.b >= lnL0) TRUE
            else rbinom(1, 1, exp(lnL.b - lnL0))
        }
        if (ACCEPT) {
            j <- j + 1L
            lnL0 <- lnL.b
            tree0 <- tr.b
            para0 <- para
            bt0 <- bt
        }
    }
    if (verbose) cat("\nDone.\n")
}

.get.list.trees <- function()
    ls(envir = .coalescentMCMCenv, pattern = "^TREES_")

getMCMCtrees <- function(chain = NULL)
{
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (is.null(chain)) {
        if (!l) return(NULL)
        if (l == 1)
            return(get(list.trees, envir = .coalescentMCMCenv))
        ## l > 1:
        cat("Several lists of MCMC trees are stored:\n\n")
        for (i in 1:l) cat(i, ":", list.trees[i], "\n")
        cat("\nReturn which number? ")
        chain <- as.numeric(readLines(n = 1))
    } else {
        if (!l) {
            warning("no list of MCMC trees stored")
            return(NULL)
        }
        if (l < chain) {
            warning("no enough lists of MCMC trees stored")
            return(NULL)
        }
    }
    get(paste("TREES", sprintf("%03d", chain), sep = "_"),
        envir = .coalescentMCMCenv)
}

saveMCMCtrees <- function(destdir = ".", format = "RDS", ...)
{
    format <- match.arg(toupper(format), c("RDS", "NEWICK", "NEXUS"))
    switch(format, RDS = {
        FUN <- saveRDS
        suffix <- ".rds"
    }, NEWICK = {
        FUN <- write.tree
        suffix <- ".tre"
    }, NEXUS = {
        FUN <- write.nexus
        suffix <- ".nex"
    })
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (!l) warning("no list of trees to save") else {
        for (i in 1:l) {
            f <- list.trees[i]
            outfile <- paste(destdir, "/", f, suffix, sep = "")
            FUN(get(f, envir = .coalescentMCMCenv), outfile, ...)
        }
    }
}

cleanMCMCtrees <- function()
    rm(list = .get.list.trees(), envir = .coalescentMCMCenv)

getLastTree <- function(X) X[[length(X)]]

getMCMCstats <- function()
{
    cat("MCMC chain summaries (chains as columns):\n\n")
    get("MCMCstats", envir = .coalescentMCMCenv)
}
