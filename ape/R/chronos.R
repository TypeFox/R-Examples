## chronos.R (2014-05-15)

##   Molecular Dating With Penalized and Maximum Likelihood

## Copyright 2013-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.chronos.ctrl <-
    list(tol = 1e-8, iter.max = 1e4, eval.max = 1e4, nb.rate.cat = 10,
         dual.iter.max = 20)

makeChronosCalib <-
    function(phy, node = "root", age.min = 1, age.max = age.min,
             interactive = FALSE, soft.bounds = FALSE)
{
    n <- Ntip(phy)
    if (interactive) {
        plot(phy)
        cat("Click close to a node and enter the ages (right-click to exit)\n\n")
        node <- integer()
        age.min <- age.max <- numeric()
        repeat {
            ans <- identify(phy, quiet = TRUE)
            if (is.null(ans)) break
            NODE <- ans$nodes
            nodelabels(node = NODE, col = "white", bg = "blue")
            cat("constraints for node ", NODE, sep = "")
            cat("\n  youngest age: ")
            AGE.MIN <- as.numeric(readLines(n = 1))
            cat("  oldest age (ENTER if not applicable): ")
            AGE.MAX <- as.numeric(readLines(n = 1))
            node <- c(node, NODE)
            age.min <- c(age.min, AGE.MIN)
            age.max <- c(age.max, AGE.MAX)
        }
        s <- is.na(age.max)
        if (any(s)) age.max[s] <- age.min[s]
    } else {
        if (identical(node, "root")) node <- n + 1L
    }

    if (any(node <= n))
        stop("node numbers should be greater than the number of tips")

    diff.age <- which(age.max < age.min)
    if (length(diff.age)) {
        msg <- "'old age' less than 'young age' for node"
        if (length(diff.age) > 1) msg <- paste(msg, "s", sep = "")
        stop(paste(msg, paste(node[diff.age], collapse = ", ")))
    }

    data.frame(node, age.min, age.max, soft.bounds = soft.bounds)
}

chronos.control <- function(...)
{
    dots <- list(...)
    x <- .chronos.ctrl
    if (length(dots)) {
        chk.nms <- names(dots) %in% names(x)
        if (any(!chk.nms)) {
            warning("some control parameter names do not match: they were ignored")
            dots <- dots[chk.nms]
        }
        x[names(dots)] <- dots
    }
    x
}

chronos <-
    function(phy, lambda = 1, model = "correlated", quiet = FALSE,
             calibration = makeChronosCalib(phy),
             control = chronos.control())
{
    model <- match.arg(tolower(model), c("correlated", "relaxed", "discrete"))
    n <- Ntip(phy)
    ROOT <- n + 1L
    m <- phy$Nnode
    el <- phy$edge.length
    if (any(el < 0)) stop("some branch lengths are negative")
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    N <- length(e1)
    TIPS <- 1:n
    EDGES <- 1:N

    tol <- control$tol

    node <- calibration$node
    age.min <- calibration$age.min
    age.max <- calibration$age.max

    if (model == "correlated") {
### `basal' contains the indices of the basal edges
### (ie, linked to the root):
        basal <- which(e1 == ROOT)
        Nbasal <- length(basal)

### 'ind1' contains the index of all nonbasal edges, and 'ind2' the
### index of the edges where these edges come from (ie, they contain
### pairs of contiguous edges), eg:

###         ___b___    ind1  ind2
###        |           |   ||   |
### ___a___|           | b || a |
###        |           | c || a |
###        |___c___    |   ||   |

        ind1 <- EDGES[-basal]
        ind2 <- match(e1[EDGES[-basal]], e2)
    }

    age <- numeric(n + m)

### This bit sets 'ini.time' and should result in no negative branch lengths

    if (!quiet) cat("\nSetting initial dates...\n")
    seq.nod <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)

    ii <- 1L
    repeat {
        ini.time <- age
        ini.time[ROOT:(n + m)] <- NA

        ini.time[node] <-
            if (is.null(age.max)) age.min
            else runif(length(node), age.min, age.max) # (age.min + age.max) / 2

        ## if no age given for the root, find one approximately:
        if (is.na(ini.time[ROOT]))
            ini.time[ROOT] <- if (is.null(age.max)) 3 * max(age.min) else 3 * max(age.max)

        ISnotNA.ALL <- unlist(lapply(seq.nod, function(x) sum(!is.na(ini.time[x]))))
        o <- order(ISnotNA.ALL, decreasing = TRUE)

        for (y in seq.nod[o]) {
            ISNA <- is.na(ini.time[y])
            if (any(ISNA)) {
                i <- 2L # we know the 1st value is not NA, so we start at the 2nd one
                while (i <= length(y)) {
                    if (ISNA[i]) { # we stop at the next NA
                        j <- i + 1L
                        while (ISNA[j]) j <- j + 1L # look for the next non-NA
                        nb.val <- j - i
                        by <- (ini.time[y[i - 1L]] - ini.time[y[j]]) / (nb.val + 1)
                        ini.time[y[i:(j - 1L)]] <- ini.time[y[i - 1L]] - by * seq_len(nb.val)
                        i <- j + 1L
                    } else i <- i + 1L
                }
            }
        }
        if (all(ini.time[e1] - ini.time[e2] >= 0)) break
        ii <- ii + 1L
        if (ii > 1000)
            stop("cannot find reasonable starting dates after 1000 tries:
maybe you need to adjust the calibration dates")
    }
### 'ini.time' set

    #ini.time[ROOT:(n+m)] <- branching.times(chr.dis)
    ## ini.time[ROOT:(n+m)] <- ini.time[ROOT:(n+m)] + rnorm(m, 0, 5)
    #print(ini.time)


### Setting 'ini.rate'
    ini.rate <- el/(ini.time[e1] - ini.time[e2])

    if (model == "discrete") {
        Nb.rates <- control$nb.rate.cat
        minmax <- range(ini.rate)
        if (Nb.rates == 1) {
            ini.rate <- sum(minmax)/2
        } else {
            inc <- diff(minmax)/Nb.rates
            ini.rate <- seq(minmax[1] + inc/2, minmax[2] - inc/2, inc)
            ini.freq <- rep(1/Nb.rates, Nb.rates - 1)
            lower.freq <- rep(0, Nb.rates - 1)
            upper.freq <- rep(1, Nb.rates - 1)
        }
    } else Nb.rates <- N
## 'ini.rate' set

### Setting bounds for the node ages

    ## `unknown.ages' will contain the index of the nodes of unknown age:
    unknown.ages <- 1:m + n

    ## initialize vectors for all nodes:
    lower.age <- rep(tol, m)
    upper.age <- rep(1/tol, m)

    lower.age[node - n] <- age.min
    upper.age[node - n] <- age.max

    ## find nodes known within an interval:
    ii <- which(age.min != age.max)
    ## drop them from 'node' since they will be estimated:
    if (length(ii)) {
        node <- node[-ii]
        if (length(node))
            age[node] <- age.min[-ii] # update 'age'
    } else age[node] <- age.min

    ## finally adjust the 3 vectors:
    if (length(node)) {
        unknown.ages <- unknown.ages[n - node] # 'n - node' is simplification for '-(node - n)'
        lower.age <- lower.age[n - node]
        upper.age <- upper.age[n - node]
    }
### Bounds for the node ages set

    ## 'known.ages' contains the index of all nodes
    ## (internal and terminal) of known age:
    known.ages <- c(TIPS, node)

    ## the bounds for the rates:
    lower.rate <- rep(tol, Nb.rates)
    upper.rate <- rep(100 - tol, Nb.rates) # needs to be adjusted to higher values?

### Gradient
    degree_node <- tabulate(phy$edge)
    eta_i <- degree_node[e1]
    eta_i[e2 <= n] <- 1L
    ## eta_i[i] is the number of contiguous branches for branch 'i'

    ## use of a list of indices is slightly faster than an incidence matrix
    ## and takes much less memory (60 Kb vs. 8 Mb for n = 500)
    X <- vector("list", N)
    for (i in EDGES) {
        j <- integer()
        if (e1[i] != ROOT) j <- c(j, which(e2 == e1[i]))
        if (e2[i] >= n) j <- c(j, which(e1 == e2[i]))
        X[[i]] <- j
    }
    ## X is a list whose i-th element gives the indices of the branches
    ## that are contiguous to branch 'i'

    ## D_ki and A_ki are defined in the SI of the paper
    D_ki <- match(unknown.ages, e2)
    A_ki <- lapply(unknown.ages, function(x) which(x == e1))

    gradient.poisson <- function(rate, node.time) {
        age[unknown.ages] <- node.time
        real.edge.length <- age[e1] - age[e2]

        ## gradient for the rates:
        gr <- el/rate - real.edge.length

        ## gradient for the dates:
        tmp <- el/real.edge.length - rate
        tmp2 <- tmp[D_ki]
        tmp2[is.na(tmp2)] <- 0
        gr.dates <- sapply(A_ki, function(x) sum(tmp[x])) - tmp2

        c(gr, gr.dates)
    }

    ## gradient of the penalized lik (must be multiplied by -1 before calling nlminb)
    gradient <-
        switch(model,
               "correlated" =
               function(rate, node.time) {
                   gr <- gradient.poisson(rate, node.time)
                   #if (all(gr == 0)) return(gr)

                   ## contribution of the penalty for the rates:
                   gr[RATE] <- gr[RATE] - lambda * 2 * (eta_i * rate - sapply(X, function(x) sum(rate[x])))
                   ## the contribution of the root variance term:
                   if (Nbasal == 2) { # the simpler formulae if there's a basal dichotomy
                       i <- basal[1]
                       j <- basal[2]
                       gr[i] <- gr[i] - lambda * (rate[i] - rate[j])
                       gr[j] <- gr[j] - lambda * (rate[j] - rate[i])
                   } else { # the general case
                       for (i in 1:Nbasal)
                           j <- basal[i]
                           gr[j] <- gr[j] -
                               lambda*2*(rate[j]*(1 - 1/Nbasal) - sum(rate[basal[-i]])/Nbasal)/(Nbasal - 1)
                   }
                   gr
               },
               "relaxed" =
               function(rate, node.time) {
                   gr <- gradient.poisson(rate, node.time)
                   #if (all(gr == 0)) return(gr)

                   ## contribution of the penalty for the rates:
                   mean.rate <- mean(rate)
                   ## rank(rate)/Nb.rates is the same than ecdf(rate)(rate) but faster
                   gr[RATE] <- gr[RATE] + lambda*2*dgamma(rate, mean.rate)*(rank(rate)/Nb.rates - pgamma(rate, mean.rate))
                   gr
               },
               "discrete" = NULL)

    log.lik.poisson <- function(rate, node.time) {
        age[unknown.ages] <- node.time
        real.edge.length <- age[e1] - age[e2]
        if (isTRUE(any(real.edge.length < 0))) return(-1e100)
        B <- rate * real.edge.length
        sum(el * log(B) - B - lfactorial(el))
    }

### penalized log-likelihood
    penal.loglik <-
        switch(model,
               "correlated" =
               function(rate, node.time) {
                   loglik <- log.lik.poisson(rate, node.time)
                   if (!is.finite(loglik)) return(-1e100)
                   loglik - lambda * (sum((rate[ind1] - rate[ind2])^2)
                                      + var(rate[basal]))
               },
               "relaxed" =
               function(rate, node.time) {
                   loglik <- log.lik.poisson(rate, node.time)
                   if (!is.finite(loglik)) return(-1e100)
                   mu <- mean(rate)
                   ## loglik - lambda * sum((1:N/N - pbeta(sort(rate), mu/(1 + mu), 1))^2) # avec loi beta
                   ## loglik - lambda * sum((1:N/N - pcauchy(sort(rate)))^2) # avec loi Cauchy
                   loglik - lambda * sum((1:N/N - pgamma(sort(rate), mean(rate)))^2) # avec loi Gamma
               },
               "discrete" =
               if (Nb.rates == 1)
                   function(rate, node.time) log.lik.poisson(rate, node.time)
               else function(rate, node.time, freq) {
                   if (isTRUE(sum(freq) > 1)) return(-1e100)
                   rate.freq <- sum(c(freq, 1 - sum(freq)) * rate)
                   log.lik.poisson(rate.freq, node.time)
               })

    opt.ctrl <- list(eval.max = control$eval.max, iter.max = control$iter.max)

    ## the following capitalized vectors give the indices of
    ## the parameters once they are concatenated in 'p'
    RATE <- 1:Nb.rates
    AGE <- Nb.rates + 1:length(unknown.ages)

    if (model == "discrete") {
        if (Nb.rates == 1) {
            start.para <- c(ini.rate, ini.time[unknown.ages])
            f <- function(p) -penal.loglik(p[RATE], p[AGE])
            g <- NULL
            LOW <- c(lower.rate, lower.age)
            UP <- c(upper.rate, upper.age)
        } else {
            FREQ <- length(RATE) + length(AGE) + 1:(Nb.rates - 1)
            start.para <- c(ini.rate, ini.time[unknown.ages], ini.freq)
            f <- function(p) -penal.loglik(p[RATE], p[AGE], p[FREQ])
            g <- NULL
            LOW <- c(lower.rate, lower.age, lower.freq)
            UP <- c(upper.rate, upper.age, upper.freq)
        }
    } else {
        start.para <- c(ini.rate, ini.time[unknown.ages])
        f <- function(p) -penal.loglik(p[RATE], p[AGE])
        g <- function(p) -gradient(p[RATE], p[AGE])
        LOW <- c(lower.rate, lower.age)
        UP <- c(upper.rate, upper.age)
    }

    k <- length(LOW) # number of free parameters

    if (!quiet) cat("Fitting in progress... get a first set of estimates\n")

    out <- nlminb(start.para, f, g,
                  control = opt.ctrl, lower = LOW, upper = UP)

    if (model == "discrete") {
        if (Nb.rates == 1) {
            f.rates <- function(p) -penal.loglik(p, current.ages)
            f.ages <- function(p) -penal.loglik(current.rates, p)
        } else {
            f.rates <- function(p) -penal.loglik(p, current.ages, current.freqs)
            f.ages <- function(p) -penal.loglik(current.rates, p, current.freqs)
            f.freqs <- function(p) -penal.loglik(current.rates, current.ages, p)
            g.freqs <- NULL
        }
        g.rates <- NULL
        g.ages <- NULL
    } else {
        f.rates <- function(p) -penal.loglik(p, current.ages)
        g.rates <- function(p) -gradient(p, current.ages)[RATE]
        f.ages <- function(p) -penal.loglik(current.rates, p)
        g.ages <- function(p) -gradient(current.rates, p)[AGE]
    }

    current.ploglik <- -out$objective
    current.rates <- out$par[RATE]
    current.ages <- out$par[AGE]
    if (model == "discrete" && Nb.rates > 1) current.freqs <- out$par[FREQ]

    dual.iter.max <- control$dual.iter.max
    i <- 0L

    if (!quiet) cat("         Penalised log-lik =", current.ploglik, "\n")

    repeat {
        if (dual.iter.max < 1) break
        if (!quiet) cat("Optimising rates...")
        out.rates <- nlminb(current.rates, f.rates, g.rates,# h.rates,
                            control = list(eval.max = 1000, iter.max = 1000,
                                           step.min = 1e-8, step.max = .1),
                            lower = lower.rate, upper = upper.rate)
        new.rates <- out.rates$par
        if (-out.rates$objective > current.ploglik)
            current.rates <- new.rates

        if (model == "discrete" && Nb.rates > 1) {
            if (!quiet) cat(" frequencies...")
            out.freqs <- nlminb(current.freqs, f.freqs,
                                control = list(eval.max = 1000, iter.max = 1000,
                                               step.min = .001, step.max = .5),
                                lower = lower.freq, upper = upper.freq)
            new.freqs <- out.freqs$par
        }

        if (!quiet) cat(" dates...")
        out.ages <- nlminb(current.ages, f.ages, g.ages,# h.ages,
                           control = list(eval.max = 1000, iter.max = 1000,
                                          step.min = .001, step.max = 100),
                           lower = lower.age, upper = upper.age)
        new.ploglik <- -out.ages$objective

        if (!quiet) cat("", current.ploglik, "\n")

        if (new.ploglik - current.ploglik > 1e-6 && i <= dual.iter.max) {
            current.ploglik <- new.ploglik
            current.rates <- new.rates
            current.ages <- out.ages$par
            if (model == "discrete" && Nb.rates > 1) current.freqs <- new.freqs
            out <- out.ages
            i <- i + 1L
        } else break
    }

    if (!quiet) cat("\nDone.\n")

#    browser()

    if (model == "discrete") {
        rate.freq <-
            if (Nb.rates == 1) current.rates
            else mean(c(current.freqs, 1 - sum(current.freqs)) * current.rates)
        logLik <- log.lik.poisson(rate.freq, current.ages)
        PHIIC <- list(logLik = logLik, k = k, PHIIC = - 2 * logLik + 2 * k)
    } else {
        logLik <- log.lik.poisson(current.rates, current.ages)
        PHI <- switch(model,
                      "correlated" = (current.rates[ind1] - current.rates[ind2])^2 + var(current.rates[basal]),
                      "relaxed" = (1:N/N - pgamma(sort(current.rates), mean(current.rates)))^2) # avec loi Gamma
        PHIIC <- list(logLik = logLik, k = k, lambda = lambda,
                      PHIIC = - 2 * logLik + 2 * k + lambda * svd(PHI)$d)
    }

    attr(phy, "call") <- match.call()
    attr(phy, "ploglik") <- -out$objective
    attr(phy, "rates") <- current.rates #out$par[EDGES]
    if (model == "discrete" && Nb.rates > 1)
        attr(phy, "frequencies") <- current.freqs
    attr(phy, "message") <- out$message
    attr(phy, "PHIIC") <- PHIIC
    age[unknown.ages] <- current.ages #out$par[-EDGES]
    phy$edge.length <- age[e1] - age[e2]
    class(phy) <- c("chronos", class(phy))
    phy
}

print.chronos <- function(x, ...)
{
    cat("\n    Chronogram\n\n")
    cat("Call: ")
    print(attr(x, "call"))
    cat("\n")
    NextMethod("print")
}
