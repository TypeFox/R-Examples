# Calculate divergence between two probability density functions using the
# Kullback-Leibler distance measure.
# h - Returns matrix
# window - number of points to sample (defaults to anylength(h))
# count - number of observations to create
# filter - function to apply
# replace - whether to use replacement in bootstrapping
# divergence(sp500.subset, 25, filter=RandomMatrixDenoiser())
# divergence(sp500.subset, 25, filter=ShrinkageDenoiser())
# Can measure information (the default) or stability. Measuring stability will
# resample twice to get two forms of the correlation matrix.

#### TODO - Rethink this API
deform(object, type) %::% AssetReturns : character : matrix
deform(object, type) %when% {
  type == 'matrix'
} %as% {
  col.names <- colnames(object)
  row.names <- format(index(object), '%Y-%m-%d')
  h <- matrix(object, ncol=ncol(object))
  colnames(h) <- col.names
  rownames(h) <- row.names
  h
}

deform(object, type) %::% zoo : character : matrix
deform(object, type) %when% {
  type == 'matrix'
} %as% {
  col.names <- colnames(object)
  row.names <- format(index(object), '%Y-%m-%d')
  h <- matrix(object, ncol=ncol(object))
  colnames(h) <- col.names
  rownames(h) <- row.names
  h
}

KullbackLeibler(...) %as% { list(...) }

divergence(ret, count, filter) %when% {
  ret %isa% zoo
} %as% {
  p <- TawnyPortfolio(ret, as.numeric(nrow(ret)))
  divergence(p, count, filter)
}

divergence(p, count, filter) %when% {
  p %hasa% returns
  is.function(filter)
} %as% {
  divergence(p, count, filter, KullbackLeibler(measure='information'))
}

divergence(p, count, filter, algo) %when% {
  p %hasa% returns
  algo %isa% KullbackLeibler
  algo$measure=='information'
} %as% {
  flog.debug("Row names: %s", rownames(p$returns))
  # Convert to matrix to allow duplicates
  h <- deform(p$returns,'matrix')
  if (is.null(p$window)) { p$window <- anylength(h) }

  # This should just use rollapply on a TawnyPortfolio
  div <- function(junk, h.full)
  {
    h.window <- h.full[sample(index(h.full), p$window, replace=TRUE), ]
    c.sample <- cov2cor(cov.sample(h.window))
    p <- TawnyPortfolio(zoo(h.window, rownames(h.window)), p$window)
    c.model <- filter(p)

    divergence <- divergence.kl(c.sample, c.model)
    return(divergence)
  }
  ds <- sapply(1:count, div, h)

  theory <- divergence_lim(ncol(h), p$window, algo)
  flog.trace("Theoretical divergence is %s",theory)

  return(c(mean=mean(ds, na.rm=TRUE), sd=sd(ds, na.rm=TRUE), limit=theory))
}

#divergence <- function(h, count, window=NULL, filter=getCorFilter.RMT(), 
#  measure='information')
#{
#  fn <- paste('divergence', measure, sep='.')
#  do.call(fn, list(h, count, window, filter))
#}

#divergence.information <- function(h, count, window, filter)
#{
#  if (is.null(window)) { window <- anylength(h) }
#  # Convert to matrix to allow duplicates
#  col.names <- colnames(h)
#  row.names <- format(index(h), '%Y-%m-%d')
#  h <- matrix(h, ncol=ncol(h))
#  colnames(h) <- col.names
#  rownames(h) <- row.names
#
#  div <- function(junk, h.full)
#  {
#    h.window <- h.full[sample(index(h.full), window, replace=TRUE), ]
#    c.sample <- cov2cor(cov.sample(h.window))
#    c.model <- filter(h.window)
#
#    divergence <- divergence.kl(c.sample, c.model)
#    return(divergence)
#  }
#  ds <- sapply(1:count, div, h)
#
#  theory <- divergenceLimit.kl(ncol(h), window)
#  #cat("Theoretical divergence is",theory,"\n")
#
#  return(c(mean=mean(ds, na.rm=TRUE), sd=sd(ds, na.rm=TRUE), limit=theory))
#}

# Measuring information compares sample correlation matrix with filtered
# correlation matrix
# Measuring stability averages all permutations of the KL divergence of two
# instances of the filtered correlation matrix
divergence.kl(sigma.1, sigma.2) %as%
{
  term.1 <- log(det(sigma.2) / det(sigma.1))
  term.2 <- sum(diag(solve(sigma.2) %*% sigma.1))
  0.5 * (term.1 + term.2 - nrow(sigma.1))
}

# The expected value of the divergence for random matrices (sample versus 
# true correlation matrix)
divergence_lim(ps, model) %when% {
  model %isa% KullbackLeibler
} %as% {
  divergence_lim(ps[1], ps[2], model)
}

divergence_lim(m, t, model) %when% {
  model %isa% KullbackLeibler
} %as% {
  l <- t - m + 1
  0.5 * ( m * log(t/2) - sum(digamma((l:t)/2)) )
}

# plotDivergenceLimit.kl(100, 80:499, col='green', ylim=c(0,55))
# plotDivergenceLimit.kl(80, 80:499, col='orange', overlay=TRUE)
# plotDivergenceLimit.kl(40, 80:499, col='red', overlay=TRUE)
plotDivergenceLimit.kl(m, t.range, ..., overlay=FALSE) %as%
{
  model <- KullbackLeibler()
  ns <- rep(m,length(t.range))
  limit <- apply(matrix(c(ns, t.range), ncol=2), 1, 
    function(m) divergence_lim(m, model))
  if (! overlay)
  {
    ylab <- 'Expected KL divergence'
    plot(t.range/m, limit, ylab=ylab, xlab='Q', type='l', ...)
  }
  else
  {
    lines(t.range/m, limit, type='l', ...)
  }

  invisible(limit)
}

# Limit for stability (distance between two sample correlation matrices)
stability_lim(m, model) %::% a : KullbackLeibler : a
stability_lim(m, model) %as%
{
  t <- m[2]
  m <- m[1]
  stability_lim(m, t, model)
}

stability_lim(m, t, model) %::% a : a : KullbackLeibler : a
stability_lim(m, t, model) %as%
{ 
  0.5 * m * (m+1) / (t - m - 1)
}

# Determine the stability of the filter.
divergence.stability(h, count, window, filter) %as%
{
  if (is.null(window)) { window <- anylength(h) }
  # Convert to matrix to allow duplicates
  h <- matrix(h, ncol=ncol(h))

  div <- function(junk, h.full)
  {
    h.window.1 <- h.full[sample(index(h.full), window, replace=TRUE), ]
    c.sample.1 <- cov2cor(cov.sample(h.window.1))
    h.window.2 <- h.full[sample(index(h.full), window, replace=TRUE), ]
    c.sample.2 <- cov2cor(cov.sample(h.window.2))

    divergence <- divergence.kl(c.sample.1, c.sample.2)
    return(divergence)
  }
  ds <- sapply(1:count, div, h)

  theory <- stability_lim(ncol(h), window)
  #cat("Theoretical divergence is",theory,"\n")

  return(c(mean=mean(ds, na.rm=TRUE), sd=sd(ds, na.rm=TRUE), limit=theory))
}

