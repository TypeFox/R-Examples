
as.rvsummary.bugs <- function (x, list.=TRUE, ...) {
  if (list.) {
    lapply(as.rv(x, list.=TRUE), as.rvsummary, ...)
  } else {
    as.rvsummary(as.rv(x, list.=FALSE), ...)
  }
}

as.rv.bugs <- function (x, list.=TRUE, ...) {
  # as.rv.bugs - transform a R2WinBUGS object into a random variable object
  # n.chains,
  # n.iter
  # n.burnin
  # n.thin
  # n.keep
  # n.sims
  # sims.array=sims.array[,,rank.long]
  # sims.list
  # sims.matrix=sims[,rank.long]
  # summary
  # mean=summary.mean
  # sd=summary.sd
  # median=summary.median
  # root.short
  # long.short
  # dimension.short
  # indexes.short
  # last.values
  #     sims <- c(bugs.sims(parameters.to.save, n.chains, n.iter, 
  #          n.burnin, n.thin, DIC), model.file = model.file, 
  #          is.DIC = DIC)
  #      class(sims) <- "bugs"
  #      return(sims)
  if (is.null(x['sims.array']))
    stop("Argument does not contain a 'sims.array'")
  A <- x$sims.array
  dnA <- dimnames(A)
  dim(A) <- c(prod(dim(A)[1:2]), dim(A)[3])
  dimnames(A) <- list(dnA[[1]], dnA[[3]])
  r <- rvsims(A)
  rvattr(r, 'Rhat')  <- x$summary[,'Rhat']
  rvattr(r, 'n.eff') <- x$summary[,'n.eff']
  if (list.) {
    return(splitbyname(r))
  } else {
    return(r)
  }
}

