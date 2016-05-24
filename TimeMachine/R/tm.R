# This file is part of the TimeMachine.
# Copyright (C) 2013 Gianluca Campanella <gianluca@campanella.org>
#
# The TimeMachine is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# The TimeMachine is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# theTimeMachine. If not, see <http://www.gnu.org/licenses/>.

setGeneric("pophistory", def=function(x, probs=seq(0, 1, 0.25)) UseMethod("pophistory"))
setGeneric("ttc", def=function(x) UseMethod("ttc"))
setGeneric("ttm", def=function(x) UseMethod("ttm"))

full.transitions <- function(unitary.transitions, loci) {
    unitary.transitions <- as.matrix(unitary.transitions)
    loci <- as.integer(loci)

    .Call("compute_full_transitions", unitary.transitions, loci)
}

stationary.dist <- function(transitions) {
    transitions <- as.matrix(transitions)

    .Call("compute_stationary_distribution", transitions)
}

sample.pop <- function(transitions, pi=NULL, pop.size, mu) {
    transitions <- as.matrix(transitions)
    if (!is.null(pi)) {
        pi <- as.numeric(pi)
    } else {
        pi <- stationary.dist(transitions)
    }
    pop.size <- as.integer(pop.size)
    mu <- as.numeric(mu)

    types <- nrow(transitions)
    population <- integer(types)
    population[sample(types, 1, prob=pi)] <- 2

    while (sum(population) < pop.size) {
        # Sample ancestor type
        ancestor.type <- sample(types, 1, prob=population)

        # Sample type of event
        p.mutation <- mu / (sum(population) - 1 + mu)
        u <- runif(1)
        if (u < p.mutation) {
            # Mutation
            offspring.type <- sample(types, 1, prob=transitions[,ancestor.type])
            population[ancestor.type] <- population[ancestor.type] - 1
            population[offspring.type] <- population[offspring.type] + 1
        } else {
            # Split
            population[ancestor.type] <- population[ancestor.type] + 1
        }
    }

    population
}

tm <- function(transitions, pi=NULL, population, n, mu, samples, threads=NULL) {
    transitions <- as.matrix(transitions)
    if (!is.null(pi)) {
        pi <- as.numeric(pi)
    } else {
        pi <- stationary.dist(transitions)
    }
    population <- as.integer(population)
    n <- as.integer(n)
    mu <- as.numeric(mu)
    samples <- as.integer(samples)
    if (!is.null(threads)) {
        threads <- as.integer(threads)
    }

    structure(.Call("estimate_loglik", transitions, pi, population, n, mu,
                    samples, threads), class="tm")
}

mu.mle <- function(transitions, pi=NULL, population, n, mu.int, samples,
                   threads=NULL, ...) {
    transitions <- as.matrix(transitions)
    if (!is.null(pi)) {
        pi <- as.numeric(pi)
    } else {
        pi <- stationary.dist(transitions)
    }
    population <- as.integer(population)
    n <- as.integer(n)
    mu.int <- as.numeric(mu.int)
    samples <- as.integer(samples)
    threads <- as.integer(threads)

    opt.res <- optimize(function(mu) {
        mean(tm(transitions, pi, population, n, mu, samples, threads))
    }, mu.int, maximum=TRUE, ...)

    list(mu.hat=opt.res$maximum, loglik=opt.res$objective)
}

density.tm <- function(x, ...) {
    density(x$logliks, ...)
}

mean.tm <- function(x, ...) {
    mean(x$logliks)
}

pophistory.tm <- function(x, probs=seq(0, 1, 0.25)) {
    h <- t(apply(x$coalescent.events, 1, quantile, probs=probs))
    rownames(h) <- seq(sum(x$population) - 1, max(2, x$n))
    h
}

print.tm <- function(x, ...) {
    summary(x)
}

quantile.tm <- function(x, ...) {
    quantile(x$logliks, ...)
}

summary.tm <- function(object, ..., digits=max(3, getOption("digits")-3)) {
    n.percentage <- object$n / sum(object$population)# 100
    mean.sim.time <- mean(object$simulation.times)
    mean.loglik <- mean(object)
    loglik.quantiles <- quantile(object, probs=c(0.05, 0.95))

    cat("Parameters\n")
    cat(sprintf("  Mutation rate: %.*f\n", digits, object$mu))
    cat(sprintf("  Target population size: %d (%.*f%%)\n", object$n, digits,
                n.percentage))
    cat(sprintf("  Samples: %d\n", length(object$logliks)))

    cat("Performance\n")
    cat(sprintf("  Mean simulation time: %.*f\n", digits, mean.sim.time))
    cat(sprintf("  Total computation time: %.*f\n", digits, object$total.time))

    cat("Log-likelihood\n")
    cat(sprintf("  Mean: %.*f\n", digits, mean(object)))
    cat(sprintf("  5%%-95%% quantiles: %.*f %.*f\n", digits, loglik.quantiles[1],
                digits, loglik.quantiles[2]))
}

ttc.tm <- function(x) {
    ttc.table <- tabulate(apply(x$coalescent.events, 2, diff))
    structure(ttc.table / sum(ttc.table), names=1:length(ttc.table))
}

ttm.tm <- function(x) {
    ttm.table <- tabulate(do.call(c, apply(x$coalescent.events, 2, function(y) {
        diff(setdiff(0:y[length(y)], y))
    })))
    structure(ttm.table / sum(ttm.table), names=1:length(ttm.table))
}

