# From SamplerCompare, (c) 2010 Madeleine Thompson

# This file contains glue for calling MCMC samplers written in C
# from R.  For more information, see the R help for each individual
# function and the vignette "R/C Glue in SamplerCompare" (doc/glue.pdf).

wrap.c.sampler <- function(sampler.symbol, sampler.context,
                           name, name.expression=NULL) {
  stopifnot(is.character(sampler.symbol) && length(sampler.symbol)==1)

  # Define the wrapper function.

  sampler <- function(target.dist, x0, sample.size, tuning=1) {
    stopifnot(target.dist$ndim==length(x0))
    stopifnot(is.numeric(tuning) && length(tuning)==1)

    # If the target distribution's log density function is written
    # in R, indicated by the absence of the c.log.density.and.grad
    # string in the distribution object, call sampler_glue_R_dist to
    # invoke the function named by sampler.symbol.  If it is written
    # in C, use sampler_glue_C_dist instead.

    if (is.null(target.dist$c.log.density.and.grad)) {
      stopifnot(!is.null(target.dist$log.density.and.grad))
      S <- .Call(sampler_glue_R_dist, sampler.symbol, sampler.context,
                 target.dist$log.density.and.grad,
                 x0, sample.size, tuning, new.env())
    } else {
      S <- .Call(sampler_glue_C_dist, sampler.symbol, sampler.context,
                 target.dist$c.log.density.and.grad,
                 target.dist$c.context, x0, sample.size, tuning)
    }
    return(S)
  }

  # Fill in "name" and possibly "name.expression" attributes on the
  # wrapper function.

  attr(sampler, 'name') <- name
  if (!is.null(name.expression))
    attr(sampler, 'name.expression') <- name.expression

  return(sampler)
}

# Returns a function pointer for the specified symbol, stored as
# an R raw vector.

raw.symbol <- function(symbol) {
  .Call(raw_symbol, symbol)
}

# Wrap the ARMS C implementation.

arms.sample <- wrap.c.sampler("arms_sample", NULL, "ARMS")

# Wrap the shrinking rank slice sampler.  It has to call wrap.c.sampler
# each time it is invoked because wrap.c.sampler does not provide a
# mechanism for mapping arguments to the sampler function into the
# sample context.

shrinking.rank.sample <- function(target.dist, x0, sample.size,
                                    tuning=1, downscale=0.95, min.dimension=1) {
  sub <- wrap.c.sampler("transition_sample",
                        list(raw.symbol("sr_draw"), downscale,
                             as.integer(min.dimension)),
                        "c-shrinking-rank-sub")
  sub(target.dist, x0, sample.size, tuning)

}
attr(shrinking.rank.sample, 'name') <- 'Shrinking Rank'

# The nonadaptive crumb sampler is implemented as a shrinking rank
# slice sampler that never shrinks rank.

nonadaptive.crumb.sample <- function(target.dist, x0, sample.size,
    tuning=1, downscale=0.95) {

  shrinking.rank.sample(target.dist, x0, sample.size, tuning,
                        downscale=downscale,
                        min.dimension=target.dist$ndim)
}

attr(nonadaptive.crumb.sample, 'name') <- 'Nonadaptive Crumb'
