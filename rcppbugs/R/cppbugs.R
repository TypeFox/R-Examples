###########################################################################
## Copyright (C) 2012  Whit Armstrong                                    ##
##                                                                       ##
## This program is free software: you can redistribute it and#or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http:##www.gnu.org#licenses#>. ##
###########################################################################

print.mcmc.object <- function(x,digits = NULL, quote = TRUE, na.print = NULL, print.gap = NULL,
    right = FALSE, max = NULL, useSource = TRUE, ...) {

    x <- unclass(x)
    xd <- dim(x)
    attributes(x) <- NULL
    if(!is.null(xd)) {
        dim(x) <- xd
    }
    print.default(x,digits, quote, na.print, print.gap, right, max, useSource, ...)
}

logp <- function(x) {
    .Call("logp",x,parent.frame(),PACKAGE="rcppbugs")
}

## print.mcmc <- function(x) {
##     invisible(.Call("printMCMC",x,PACKAGE="rcppbugs"))
## }

## print.arma <- function(x) {
##     invisible(.Call("printArma",x,PACKAGE="rcppbugs"))
## }

create.model <- function(...) {
    m <- match.call()
    attr(m,"env") <- parent.frame()
    class(m) <- "mcmc.model"
    m
}

run.model <- function(m, iterations, burn, adapt, thin) {
    if(adapt != 0 && adapt < 200) {
        stop("if adapt != 0, it must be at least 200.  Turn adapt off by setting it adapt = 0.")
    }
    .Call("runModel", m, iterations, burn, adapt, thin, PACKAGE="rcppbugs")
}

get.ar <- function(x) {
    if(is.null(attr(x,"acceptance.ratio"))) {
        stop("x is not a 'cppbugs.trace' object.")
    }
    attr(x,"acceptance.ratio")
}

deterministic <- function(f,...) {
    mc <- match.call()
    stopifnot(typeof(eval(mc[[2]]))=="closure")
    ## capture shape/type of result
    x <- do.call(f,list(...))
    attr(x,"distributed") <- "deterministic"
    attr(x,"update.method") <- compiler::compile(f)
    attr(x,"call") <- mc
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

linear <- function(X,b) {
    if(missing(X)) stop("required argument 'X' missing.")
    if(missing(b)) stop("required argument 'b' missing.")
    stopifnot(is.null(dim(b)))
    stopifnot(!is.null(dim(X)))
    stopifnot(length(dim(X))==2L)
    stopifnot(length(b)==ncol(X))
    x <- X %*% b
    attr(x,"distributed") <- "linear.deterministic"
    attr(x,"X") <- substitute(X)
    attr(x,"b") <- substitute(b)
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

linear.grouped <- function(X,b,group) {
    if(missing(X)) stop("required argument 'X' missing.")
    if(missing(b)) stop("required argument 'b' missing.")
    if(missing(group)) stop("required argument 'group' missing.")
    if(is.null(dim(b))) { stop ("'b' must be a matrix.") }
    if(is.null(dim(X))) { stop("'X' must be a matrix.") }
    stopifnot(length(dim(X))==2L)
    if(ncol(b)!=ncol(X)) { stop("the number of columns of X must match the number of columns of b") }
    if(length(group)!=nrow(X)) { stop("the length of 'group' must match the number of rows of 'X'.") }
    x <- as.matrix(apply(X * b[group,],1,sum))
    attr(x,"distributed") <- "linear.grouped.deterministic"
    attr(x,"X") <- substitute(X)
    attr(x,"b") <- substitute(b)
    attr(x,"group") <- substitute(group)
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}


logistic <- function(X,b) {
    if(missing(X)) stop("required argument 'X' missing.")
    if(missing(b)) stop("required argument 'b' missing.")
    stopifnot(is.null(dim(b)))
    stopifnot(!is.null(dim(X)))
    stopifnot(length(dim(X))==2L)
    stopifnot(length(b)==ncol(X))
    x <- 1/(1 + exp(-X %*% b))
    attr(x,"distributed") <- "logistic.deterministic"
    attr(x,"X") <- substitute(X)
    attr(x,"b") <- substitute(b)
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

check.dim.eq <- function(x,hyper) {
    ## check length equality if hyper is not a scalar and x and hyper are both vectors
    if(length(hyper) != 1 && is.null(dim(hyper)) && is.null(dim(x)) && length(x) != length(hyper)) {
        stop("hyperparameter length does not match it's variable (arma does not recycle).")
    }


    ## fail if hyper is not a scalar and x has a dimension and hyper doesn't
    if(length(hyper) != 1 && !is.null(dim(x)) && is.null(dim(hyper))) {
        stop("hyperparameter does have a dimensions attribute (and the stochastic variable does).")
    }

    ## check dimension equality if both x and hyper have dims
    if(!is.null(dim(x)) && !is.null(dim(hyper))) {
        if(!all.equal(dim(x),dim(hyper))) {
            stop("dimension of hyperparameter must be the same as the stochastic variable")
        }
    }
}

mcmc.normal <- function(x,mu,tau,observed=FALSE) {
    if(missing(mu)) stop("required argument 'mu' missing.")
    if(missing(tau)) stop("required argument 'tau' missing.")
    if(length(mu) > length(x) || length(tau) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,mu)
    check.dim.eq(x,tau)

    attr(x,"distributed") <- "normal"
    attr(x,"mu") <- substitute(mu)
    attr(x,"tau") <- substitute(tau)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

mcmc.uniform <- function(x,lower,upper,observed=FALSE) {
    if(missing(lower)) stop("required argument 'lower' missing.")
    if(missing(upper)) stop("required argument 'upper' missing.")
    if(length(lower) > length(x) || length(upper) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,lower)
    check.dim.eq(x,upper)

    attr(x,"distributed") <- "uniform"
    attr(x,"lower") <- substitute(lower)
    attr(x,"upper") <- substitute(upper)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

mcmc.gamma <- function(x,alpha,beta,observed=FALSE) {
    if(missing(alpha)) stop("required argument 'alpha' missing.")
    if(missing(beta)) stop("required argument 'beta' missing.")

    if(length(alpha) > length(x) || length(beta) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,alpha)
    check.dim.eq(x,beta)

    attr(x,"distributed") <- "gamma"
    attr(x,"alpha") <- substitute(alpha)
    attr(x,"beta") <- substitute(beta)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

mcmc.beta <- function(x,alpha,beta,observed=FALSE) {
    if(missing(alpha)) stop("required argument 'alpha' missing.")
    if(missing(beta)) stop("required argument 'beta' missing.")

    if(length(alpha) > length(x) || length(beta) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,alpha)
    check.dim.eq(x,beta)

    attr(x,"distributed") <- "beta"
    attr(x,"alpha") <- substitute(alpha)
    attr(x,"beta") <- substitute(beta)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

mcmc.bernoulli <- function(x,p,observed=FALSE) {
    if(missing(p)) stop("required argument 'p' missing.")

    if(length(p) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,p)

    attr(x,"distributed") <- "bernoulli"
    attr(x,"p") <- substitute(p)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}

mcmc.binomial <- function(x,n,p,observed=FALSE) {
    if(missing(n)) stop("required argument 'n' missing.")
    if(missing(p)) stop("required argument 'p' missing.")

    if(length(n) > length(x) || length(p) > length(x)) {
        stop("dimensions of hyperparmeters are larger than the stochastic variable itself. ",
             "(is this really what you wanted to do?)")
    }
    check.dim.eq(x,n)
    check.dim.eq(x,p)

    attr(x,"distributed") <- "binomial"
    attr(x,"n") <- substitute(n)
    attr(x,"p") <- substitute(p)
    attr(x,"observed") <- observed
    attr(x,"env") <- parent.frame()
    class(x) <- "mcmc.object"
    x
}
