######################################################################
#
# negenes.R
#
# copyright (c) 2002-2012, Karl W Broman
# last modified Mar, 2012
# first written June, 2002
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: negenes, sim.mutants
#
######################################################################

######################################################################
#
# negenes: "Number of essential genes"
#
# n.sites = (x_i) = no. tranposon sites in each gene (alone)
# counts  = (y_i) = no. mutants was observed for each gene (alone)
#
# n.sites2 = (w_i) = no. transposon sites shared by genes i and i+1
# counts2  = (z_i) = no. mutants shared by genes i and i+1
#            [in the above, take gene N+1 to be the same as gene 1
#
# n.mcmc    = number of Gibbs iterations to perform
# skip      = an integer; only save every skip+1st iteration
# burnin    = number of initial Gibbs steps to run (output discarded)
#
# startp    = Initial proportion of genes with no observed gene that will
#             be assumed essential for the Gibbs sampler.
#             If startp=0, we'll start with theta=1 for all genes
#             If startp=1, we'll start with theta=1 only for genes
#             for which a mutant was observed,
#             Otherwise, for genes for which a mutant was not
#             observed, theta ~ bernoulli(1-p), independently
#
# trace     = if TRUE, print iteration number occassionally
# calc.prob = if TRUE, return log posterior prob'y (up to scalar) for
#             each saved iteration
# return.output = if TRUE, include detailed Gibbs results in output
#
######################################################################

negenes <-
function(n.sites, counts, n.sites2, counts2,
         n.mcmc=5000, skip=49, burnin=500,
         startp=1, trace=TRUE,
         calc.prob=FALSE, return.output=FALSE)
{
  n.genes <- length(n.sites)

  # check for errors in the input
  if(length(counts) != n.genes)
    stop("n.sites and counts must be the same length")

  if(missing(n.sites2)) n.sites2 <- rep(0,n.genes)
  if(missing(counts2)) counts2 <- rep(0,n.genes)
  if(length(n.sites2) != n.genes)
    stop("n.sites2 and n.sites must be the same length")
  if(length(counts2) != n.genes)
    stop("counts2 and n.sites must be the same length")

  if(any(n.sites<0) || any(counts<0) || any(n.sites2<0) || any(counts2<0))
    stop("n.sites, counts, n.sites2, and counts2 must all be >= 0")

  # replace n.mcmc with number of iterations to be saved
  n.mcmc <- ceiling(n.mcmc/(skip+1))

  # update burnin, since using skip
  burnin <- ceiling(burnin/(skip+1))

  if(n.mcmc <= 0 || burnin < 0 || skip < 0)
    stop("n.mcmc, burnin or skip are incorrectly chosen")

  # genes known to be non-essential
  known <- as.numeric(((counts > 0) | (counts2 > 0) |
                       (c(counts2[n.genes],counts2[-n.genes])>0)))
  notknown <- ((1:n.genes)-1)[known==0]
  n.known <- sum(known)
  n.notknown <- n.genes-n.known

  if(startp < 0 || startp > 1)
    stop("startp must be between 0 and 1")

  temp <- as.numeric(return.output)

  n.mutants <- sum(counts) + sum(counts2)

  output <- .C("R_negenes",
               as.integer(n.genes),
               as.integer(n.mutants),
               as.integer(n.sites),
               # n.sites2 made so that first == last and last == first
               as.integer(c(n.sites2[length(n.sites2)],n.sites2,n.sites2[1])),
               as.integer(known),
               as.integer(n.mcmc),
               as.integer(burnin),
               as.integer(skip),
               output = as.integer(rep(0,(n.mcmc-1)*n.genes*temp+n.genes)),
               n.ess = as.integer(rep(0,n.mcmc)),
               geneprob = as.double(rep(0,n.genes)),
               as.integer(rep(0,n.genes+2)),
               as.integer(n.notknown),
               as.integer(notknown),
               as.integer(calc.prob),
               logprob = as.double(rep(0,n.mcmc)),
               as.integer(return.output),
               as.integer(trace),
               as.double(startp),
               PACKAGE="negenes")

  logprob <- output$logprob
  tot.ess <- output$n.ess
  geneprob <- output$geneprob
  if(return.output) {
    output <- matrix(output$output,ncol=n.genes)
    storage.mode(output) <- "integer"
  }
  summ <- c(mean=sum(geneprob),sd=sd(tot.ess),
            quantile(tot.ess,c(0.025,0.975)))

  if(return.output) {
    if(!calc.prob)
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,output=output))
    else
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,
                  logprob=logprob, output=output))
  }
  else {
    if(!calc.prob)
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob))
    else
      return(list(n.essential=tot.ess, summary=summ,
                  geneprob=geneprob,logprob=logprob))
  }
}





######################################################################
#
# sim.mutants
#
# Simulate mutant count data
#
######################################################################

sim.mutants <-
function(n.sites, essential, n.sites2, n.mutants)
{
  n.genes <- length(n.sites)
  if(length(essential) != n.genes)
    stop("n.sites and essential must be the same length")

  if(missing(n.sites2)) n.sites2 <- rep(0,n.genes)

  if(length(n.sites2) != n.genes)
    stop("n.sites and n.sites2 must be the same length")

  if(any(essential != 0 & essential != 1))
    stop("essential must contain only 0's and 1's.")

  if(n.mutants <= 0)
    stop("n.mutants must be positive")

  temp <- c(essential[-1],essential[1])
  p <- c(n.sites*(1-essential), n.sites2*(1-essential)*(1-temp))

  o <- table(factor(sample(1:(2*n.genes), n.mutants, replace=TRUE,
                           prob=p), levels=1:(2*n.genes)))
  names(o) <- NULL

  if(sum(n.sites2)==0) return(o[1:n.genes])
  else return(cbind(o[1:n.genes],o[-(1:n.genes)]))
}

# end of negenes.R
