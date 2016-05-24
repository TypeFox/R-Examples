##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

sample0 <- function(x, ...) x[sample.int(length(x), ...)]
bsample <- function(x, ...) x[sample.int(length(x), replace = TRUE, ...)]

isMatrix <- function(x) length(dim(x)) == 2L

pboot <- function(data, statistic, R, sim = "ordinary",
                  stype = c("i", "f", "w"),
                  strata  =  rep(1, n), L = NULL, m = 0, weights = NULL,
                  ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...)
{  
#
# R replicates of bootstrap applied to  statistic(data)
# Possible sim values are: "ordinary", "balanced", "antithetic",
#                     "permutation", "parametric"
# Various auxilliary functions find the indices to be used for the
# bootstrap replicates and then this function loops over those replicates.
#
  call <- match.call()
  stype <- match.arg(stype)
  
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0, so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    n <- NROW(data)
  if ((n == 0) || is.null(n))
    stop("no data in call to boot")
  temp.str <- strata
  strata <- tapply(seq_len(n),as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L))
      L <- boot::empinf(data = data, statistic = statistic,
                        stype = stype, strata = strata, ...)
    if (sim != "ordinary") m <- 0
    else if (any(m < 0)) stop("negative value of m supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata))))
      stop("length of m incompatible with strata")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (isMatrix(weights) && (nrow(weights) != length(R)))
        stop("dimensions of R and weights do not match")}
    else weights <- NULL
    if (!is.null(weights))
      weights <- t(apply(matrix(weights, n, length(R), byrow = TRUE),
                         2L, normalize, strata))
    if (!simple) i <- index.array(n, R, sim, strata, m, L, weights)
    
    original <- if (stype == "f") rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    } else seq_len(n)
    
    t0 <- if (sum(m) > 0L) statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
        t0
  } else # "parametric"
  statistic(data, ...)
  
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ## force promises, so values get sent by parallel
    ran.gen; data; mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  } else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[ , (n+1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- boot::freq.array(i)
      rm(i)
      if (stype == "w") f <- f/ns
      if (sum(m) == 0L) function(r) statistic(data, f[r,  ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, ], ...)
    } else if (sum(m) > 0L)
      function(r) statistic(data, i[r, ], pred.i[r,], ...)
    else if (simple)
      function(r)
        statistic(data,
                  index.array(n, 1, sim, strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  
  res <- .Call("pboot",
               as.list(seq_len(RR)),
               fn)

  t.star <- matrix(, RR, length(t0))
  for(r in seq_len(RR)) t.star[r, ] <- res[[r]]

  if (is.null(weights)) weights <- 1/tabulate(strata)[strata]
  boot.return(sim, t0, t.star, temp.str, R, data, statistic, stype,
                     call, seed, L, m, pred.i, weights, ran.gen, mle)
}

# Copied from boot class bootfuns.q 
boot.return <- function(sim, t0, t, strata, R, data, stat, stype, call,
			seed, L, m, pred.i, weights, ran.gen, mle)
#
# Return the results of a bootstrap in the form of an object of class
# "boot".
#
{
    out <- list(t0=t0, t=t, R=R, data=data, seed=seed,
                statistic=stat, sim=sim, call=call)
    if (sim == "parametric")
        out <- c(out, list(ran.gen=ran.gen, mle=mle))
    else if (sim == "antithetic")
        out <- c(out, list(stype=stype, strata=strata, L=L))
    else if (sim == "ordinary") {
        if (sum(m) > 0)
            out <- c(out, list(stype=stype, strata=strata,
                               weights=weights, pred.i=pred.i))
        else 	out <- c(out, list(stype=stype, strata=strata,
                                   weights=weights))
    } else if (sim == "balanced")
        out <- c(out, list(stype=stype, strata=strata,
                           weights=weights ))
    else
        out <- c(out, list(stype=stype, strata=strata))
    class(out) <- "boot"
    out
}

# Copied from boot class bootfuns.q 
index.array <- function(n, R, sim, strata=rep(1,n), m=0, L=NULL, weights=NULL)
{
#
#  Driver function for generating a bootstrap index array.  This function
#  simply determines the type of sampling required and calls the appropriate
#  function.
#
    indices <- NULL
    if (is.null (weights)) {
        if (sim == "ordinary") {
            indices <- ordinary.array(n, R, strata)
            if (sum(m) > 0)
                indices <- cbind(indices, extra.array(n, R, m, strata))
        }
    else if (sim == "balanced")
        indices <- balanced.array(n, R, strata)
    else if (sim == "antithetic")
        indices <- antithetic.array(n, R, L, strata)
    else if (sim == "permutation")
        indices <- permutation.array(n, R, strata)
    } else {
        if (sim == "ordinary")
            indices <- importance.array(n, R, weights, strata)
        else if (sim == "balanced")
            indices <- importance.array.bal(n, R, weights, strata)
    }
    indices
}

# Copied from boot class bootfuns.q 
normalize <- function(wts, strata)
{
#
# Normalize a vector of weights to sum to 1 within each strata.
#
    n <- length(strata)
    out <- wts
    inds <- as.integer(names(table(strata)))
    for (is in inds) {
        gp <- seq_len(n)[strata == is]
        out[gp] <- wts[gp]/sum(wts[gp]) }
    out
}



# Copied from boot class bootfuns.q 
ordinary.array <- function(n, R, strata)
{
#
# R x n array of bootstrap indices, resampled within strata.
# This is the function which generates a regular bootstrap array
# using equal weights within each stratum.
#
    inds <- as.integer(names(table(strata)))
    if (length(inds) == 1L) {
        output <- sample.int(n, n*R, replace=TRUE)
        dim(output) <- c(R, n)
    } else {
        output <- matrix(as.integer(0L), R, n)
        for(is in inds) {
            gp <- seq_len(n)[strata == is]
            output[, gp] <- if (length(gp) == 1) rep(gp, R) else bsample(gp, R*length(gp))
        }
    }
    output
}


# Copied from boot class bootfuns.q 
extra.array <- function(n, R, m, strata=rep(1,n))
{
#
# Extra indices for predictions.  Can only be used with
# types "ordinary" and "stratified".  For type "ordinary"
# m is a positive integer.  For type "stratified" m can
# be a positive integer or a vector of the same length as
# strata.
#
    if (length(m) == 1L)
        output <- matrix(sample.int(n, m*R, replace=TRUE), R, m)
    else {
        inds <- as.integer(names(table(strata)))
        output <- matrix(NA, R, sum(m))
        st <- 0
        for (i in inds) {
            if (m[i] > 0) {
                gp <- seq_len(n)[strata == i]
                inds1 <- (st+1):(st+m[i])
                output[,inds1] <- matrix(bsample(gp, R*m[i]), R, m[i])
                st <- st+m[i]
            }
        }
    }
    output
}

# Copied from boot class bootfuns.q 
balanced.array <- function(n, R, strata)
{
#
# R x n array of bootstrap indices, sampled hypergeometrically
# within strata.
#
    output <- matrix(rep(seq_len(n), R), n, R)
    inds <- as.integer(names(table(strata)))
    for(is in inds) {
        group <- seq_len(n)[strata == is]
        if(length(group) > 1L) {
            g <- matrix(rperm(output[group,  ]), length(group), R)
            output[group,  ] <- g
        }
    }
    t(output)
}


# Copied from boot class bootfuns.q 
antithetic.array <- function(n, R, L, strata)
#
#  Create an array of indices by antithetic resampling using the
#  empirical influence values in L.  This function just calls anti.arr
#  to do the sampling within strata.
#
{
    inds <- as.integer(names(table(strata)))
    out <- matrix(0L, R, n)
    for (s in inds) {
	gp <- seq_len(n)[strata == s]
        out[, gp] <- anti.arr(length(gp), R, L[gp], gp)
    }
    out
}


# Copied from boot class bootfuns.q 
permutation.array <- function(n, R, strata)
{
#
# R x n array of bootstrap indices, permuted within strata.
# This is similar to ordinary array except that resampling is
# done without replacement in each row.
#
    output <- matrix(rep(seq_len(n), R), n, R)
    inds <- as.integer(names(table(strata)))
    for(is in inds) {
        group <- seq_len(n)[strata == is]
        if (length(group) > 1L) {
            g <- apply(output[group,  ], 2L, rperm)
            output[group,  ] <- g
        }
    }
    t(output)
}

# Copied from boot class bootfuns.q 
importance.array <- function(n, R, weights, strata){
#
#  Function to do importance resampling  within strata based
#  on the weights supplied.  If weights is a matrix with n columns
#  R must be a vector of length nrow(weights) otherwise weights
#  must be a vector of length n and R must be a scalar.
#
    imp.arr <- function(n, R, wts, inds=seq_len(n))
        matrix(bsample(inds, n*R, prob=wts), R, n)
    output <- NULL
    if (!isMatrix(weights))
        weights <- matrix(weights, nrow=1)
    inds <- as.integer(names(table(strata)))
    for (ir in seq_along(R)) {
        out <- matrix(rep(seq_len(n), R[ir]), R[ir], n, byrow=TRUE)
        for (is in inds) {
            gp <- seq_len(n)[strata == is]
            out[, gp] <- imp.arr(length(gp), R[ir],
                                 weights[ir,gp], gp)
        }
        output <- rbind(output, out)
    }
    output
}

# Copied from boot class bootfuns.q 
anti.arr <- function(n, R, L, inds=seq_len(n))
{
#  R x n array of bootstrap indices, generated antithetically
#  according to the empirical influence values in L.
    unique.rank <- function(x) {
# Assign unique ranks to a numeric vector
        ranks <- rank(x)
        if (any(duplicated(ranks))) {
            inds <- seq_along(x)
            uniq <- sort(unique(ranks))
            tab <- table(ranks)
            for (i in seq_along(uniq))
                if (tab[i] > 1L) {
                    gp <- inds[ranks == uniq[i]]
                    ranks[gp] <- rperm(inds[sort(ranks) == uniq[i]])
                }
        }
        ranks
    }
    R1 <- floor(R/2)
    mat1 <- matrix(bsample(inds, R1*n), R1, n)
    ranks <- unique.rank(L)
    rev <- inds
    for (i in seq_len(n)) rev[i] <- inds[ranks == (n+1-ranks[i])]
    mat1 <- rbind(mat1, matrix(rev[mat1], R1, n))
    if (R != 2*R1) mat1 <- rbind(mat1, bsample(inds, n))
    mat1
}

# Copied from boot class bootfuns.q 
## random permutation of x.
rperm <- function(x) sample0(x, length(x))



importance.array.bal <- function(n, R, weights, strata) {
#
#  Function to do balanced importance resampling within strata
#  based on the supplied weights.  Balancing is achieved in such
#  a way that each index appears in the array approximately in
#  proportion to its weight.
#
    imp.arr.bal <- function(n, R, wts, inds=seq_len(n)) {
        if (sum (wts) != 1) wts <- wts / sum(wts)
        nRw1 <- floor(n*R*wts)
        nRw2 <- n*R*wts - nRw1
        output <- rep(inds, nRw1)
        if (any (nRw2 != 0))
            output <- c(output,
                        sample0(inds, round(sum(nRw2)), prob=nRw2))
        matrix(rperm(output), R, n)
    }
    output <- NULL
    if (!isMatrix(weights))
        weights <- matrix(weights, nrow = 1L)
    inds <- as.integer(names(table(strata)))
    for (ir in seq_along(R)) {
        out <- matrix(rep(seq_len(n), R[ir]), R[ir], n, byrow=TRUE)
        for (is in inds) {
            gp <- seq_len(n)[strata == is]
            out[,gp] <- imp.arr.bal(length(gp), R[ir], weights[ir,gp], gp)
        }
        output <- rbind(output, out)
    }
    output
}

