#!/usr/bin/Rscript

infer.clonality <- function(
  read.count.matrix, 
  variance.method = 'fpc.max',
  estimate.abundances = F,
  num.iterations = 1,
  internal.parameters = list()
  ) {

take.pos <- function(x)
{
  x[x < 0] <- 0
  return(x)
}

if (length(internal.parameters$replicates > 0)) {
  replicates <- internal.parameters$replicates
} else {
  # normalize replicates, and get rep.gram.matrix once
  replicates <- as.matrix(read.count.matrix)
  replicates <- t(t(replicates) / colSums(replicates))
}

if (length(internal.parameters$rep.gram.matrix) > 0) {
  rep.gram.matrix <- internal.parameters$rep.gram.matrix
} else {
  rep.gram.matrix <- t(replicates) %*% replicates
}

n <- nrow(replicates)
num.replicates <- ncol(replicates)
num.pairs <- num.replicates * (num.replicates - 1) / 2

if (length(internal.parameters$simple.precision.clonality) > 0) {
  simple.precision.clonality <- internal.parameters$simple.precision.clonality
} else {
  ## simple model where each read is independent
  reads.per.replicate <- as.numeric(apply(read.count.matrix, 2, sum))
  simple.precision.weights <- matrix(reads.per.replicate, nrow = num.replicates) %*% 
    matrix(reads.per.replicate, ncol = num.replicates)
  simple.precision.weights <- lower.tri(simple.precision.weights) * simple.precision.weights
  simple.precision.clonality <- sum(simple.precision.weights * rep.gram.matrix) / 
    sum(simple.precision.weights)
}
if (num.replicates < 3) # not enough replicates
{
  return(list(
  internal.parameters = internal.parameters,
  variance.method = variance.method,
  simple.precision.clonality = simple.precision.clonality, 
  lymphclon.clonality = "Too few replicates: at least 3 are needed. 6 is recommended."
    ))
}

if (length(internal.parameters$num.clones.est) > 0) {
  num.clones.est <- internal.parameters$num.clones.est
} else {
  int.chao <- function (x) { # taken from fossil package
    so <- length(x[x > 0])
    s1 <- length(x[x == 1])
    s2 <- length(x[x == 2])
    if ((s1 - s2)^2 == (s1 + s2)^2) 
        return(so + s1 * (s1 - 1)/((s2 + 1) * 2))
    else return(so + s1^2/(s2 * 2))
  }

  num.clones.est <- int.chao(apply(replicates > 0, 1, sum))
}

if (length(internal.parameters$use.squared.err.est) > 0) {
  use.squared.err.est <- as.matrix(internal.parameters$use.squared.err.est)
} else {
  use.squared.err.est <- c()
}

if (length(internal.parameters$use.replicate.var.est) > 0) {
  use.replicate.var.est <- as.matrix(internal.parameters$use.replicate.var.est)
} else {
  use.replicate.var.est <- c()
}

compute.variances.d1jkn <- (num.replicates >= 4) # We need at least 4 to mix
if (length(internal.parameters$compute.variances.d1jkn > 0)) {
  compute.variances.d1jkn <- internal.parameters$compute.variances.d1jkn
}


internal.parameters <- list(
  replicates = replicates,
  rep.gram.matrix = rep.gram.matrix,
  simple.precision.clonality = simple.precision.clonality,
  num.clones.est = num.clones.est, 
  use.squared.err.est = use.squared.err.est,
  use.replicate.var.est = use.replicate.var.est,
  compute.variances.d1jkn = compute.variances.d1jkn
)

curr.clonality.score.estimate <- simple.precision.clonality
fpc.iter.estimates <- c()

for (curr.iter.number in 1:num.iterations) {

replicates.cov.off.diagonal.value <-
  take.pos(curr.clonality.score.estimate - (1 / num.clones.est)) 

replicates.cov.off.diagonals <- 
    matrix(replicates.cov.off.diagonal.value, num.replicates, num.replicates)
    - diag(rep(curr.clonality.score.estimate, num.replicates))

if (variance.method %in% c('fpc.add'))
{ # diagonal is the clonality score plus the abs difference of the self-inner products from it
  replicates.cov.diagonals <- 
    ifelse(
      diag(rep.gram.matrix) > replicates.cov.off.diagonal.value,
      diag(rep.gram.matrix), 
      2 * replicates.cov.off.diagonal.value - diag(rep.gram.matrix))
    + rep((1 / num.clones.est), num.replicates) 
    # regularize by smallest possible clonality, given number of clones seen  

    replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
} else if (variance.method %in% c('fpc.max')) 
{ # diagonal is the max of clonality score, or the self-inner products
  replicates.cov.diagonals <- 
    ifelse(
      diag(rep.gram.matrix) > replicates.cov.off.diagonal.value,
      diag(rep.gram.matrix), 
      replicates.cov.off.diagonal.value)
    + rep((1 / num.clones.est), num.replicates)  

  # regularize by smallest possible clonality, given number of clones seen
  replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
} else if (variance.method %in% c('mle.cov')) {
  replicates.cov <- cov(replicates) * n
} else if (variance.method %in% c('usr.var')) {
  replicates.cov.diagonals <- use.replicate.var.est
  replicates.cov <- diag(replicates.cov.diagonals) + replicates.cov.off.diagonals      
}

# usr.rer: internal.parameters$use.squared.err.est specifies conditional variances of replicates
# usr.var internal.parameters$use.replicate.var.est specifies variances of replicates
# fpc.max: fixed point covariance: average off diagonals: variances are set to their lower bound when empirically seen to be lower (labeled positive expectation in the paper)
# fpc.add: fixed point covariance: average off diagonals: variances are set to their lower bound plus the empirical difference from that bound, when empirically seen to be lower (labeled concentric shells in the paper)
# mle.cov: use the empirical variance on the read count matrix abundances, without any projections
# corpcor: corpcor covariance, from the corpcor package

Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
ptinv.Lambda.matrix <- matrix(data = NA, nrow = num.replicates, ncol = num.replicates)
if (variance.method == 'usr.rer') {
  epsilon.vec <- use.squared.err.est
  inv.eps.vec <- 1 / epsilon.vec
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else if (variance.method %in% c('fpc.add', 'fpc.max', 'mle.cov', 'corpcor', 'usr.var')) {
  if (variance.method %in% c('fpc.add', 'mle.cov', 'fpc.max', 'usr.var')) {
    inv.eps.vec <- diag(ginv(replicates.cov))
  } else { 
    # variance.method == 'corpcor'
    capture.output(inv.eps.vec <- diag(invcov.shrink(replicates)))
  }

  epsilon.vec <- 1 / inv.eps.vec # the estimated conditional errors associated with each replicate
  diag(Lambda.matrix) <- inv.eps.vec
  diag(ptinv.Lambda.matrix) <- epsilon.vec
} else {
  write(sprintf('unknown variance method: %s\n', variance.method), stderr())
  write('list of valid methods: usr.rer, fpc.add, fpc.max, corpcor', stderr())
}

contributions.to.replicate.cov.matrix <- -Lambda.matrix # negative conditional covariances
diag(contributions.to.replicate.cov.matrix) <- diag(ptinv.Lambda.matrix) # inverse conditional variances

num.estimators <- num.replicates * (num.replicates - 1) / 2

estimators.rownums <- diag(c(1:num.replicates)) %*% matrix(1, num.replicates, num.replicates)
estimators.colnums <- t(estimators.rownums)
lowtri.indx <- as.vector(lower.tri(rep.gram.matrix))
estimators.vec <- as.vector(rep.gram.matrix)[lowtri.indx]
rownums.vec <- as.vector(estimators.rownums)[lowtri.indx]
colnums.vec <- as.vector(estimators.colnums)[lowtri.indx]

row.col.to.indx <- function(r, c) {return (((r - 1) * (r - 2)) / 2 + c)} 
# this is for the lower half of cov mat, ie where r<c  r-1 choose 2 + c

# the reverse will be pre computed
indx.to.row <- rep(0, num.replicates * num.replicates)
indx.to.col <- rep(0, num.replicates * num.replicates)
curr.indx <- 0

# indx.to.row maps each integer indexing the n-choose-2 estimators to an integer row
for (r in 2:num.replicates) {
  for (c in 1:(r - 1)) {
    curr.indx <- row.col.to.indx(r, c) # could have been "++1" but just to be safe...
    indx.to.row[curr.indx] <- r
    indx.to.col[curr.indx] <- c 
  }
}

# fill out the actual estimators
estimator.vec.forcov <- rep(0, num.pairs)
for (r in 2:num.replicates) {
  for (c in 1:(r-1)) {
    curr.indx <- row.col.to.indx(r, c) 
    estimator.vec.forcov[curr.indx] <- rep.gram.matrix[r, c]
  }
}

# debug.data <- c()
full.cov <- matrix(data = 0, nrow = num.estimators, ncol = num.estimators)
# (n choose 2) by (n choose 2)
for (r1 in 2:num.replicates) {
  for (c1 in 1:(r1 - 1)) {
    for (r2 in 2:num.replicates) {
      for (c2 in 1:(r2 - 1)) {
        curr.indx1 <- row.col.to.indx(r1, c1) # could have been "++1" but just to be safe...
        curr.indx2 <- row.col.to.indx(r2, c2)
        curr.cov.components <- c(
          contributions.to.replicate.cov.matrix[r1, r2],
          contributions.to.replicate.cov.matrix[r1, c2],
          contributions.to.replicate.cov.matrix[c1, r2],
          contributions.to.replicate.cov.matrix[c1, c2])
        set.count <- sum(!is.na(curr.cov.components))
        curr.cov.components[is.na(curr.cov.components)] <- 0
        # debug.data <- rbind(debug.data, c(curr.indx1, curr.indx2, r1, c1, r2, c2, set.count, sum(curr.cov.components)))
        full.cov[curr.indx1, curr.indx2] <- sum(curr.cov.components)
      }
    }
  }
}

unreg.matrix <- full.cov

# use the covariance matrix to compute the mvg MLE estimate, given the n-choose-2 estimators
cov.to.clonality <- function(target.matrix, reg.coefs) {
  target.term <- (target.matrix * reg.coefs)
  unreg.term <- (unreg.matrix * (1 - reg.coefs))
  linear.combo.matrix <- target.term + unreg.term
  prec.matrix <- ginv(linear.combo.matrix)
  root.prec.matrix <- sqrtm(prec.matrix)
  numerator <- rep(1, num.pairs) %*% prec.matrix %*% estimator.vec.forcov
  denominator <- rep(1, num.pairs) %*% prec.matrix %*% rep(1, num.pairs)
  return(abs(numerator / denominator)) # the clonality score
} # cov.to.clonality

# regularize the n-choose-2 by n-choose-2 covariance matrix
mean.eps2 <- mean(epsilon.vec)
min.eps2 <- min(epsilon.vec)

regularization.method.names <- 
  c('unregularized', 'ue.zr.full', 
  'eq.zr.half', 'ue.zr.half', 'eq.eq.half', 'ue.eq.half',
  'ue.mn.half', 'ue.mn.full', 'ue.mn.js1')
num.reg.methods <- length(regularization.method.names)
regularized.estimates <- rep(NA, num.reg.methods)
names(regularized.estimates) <- regularization.method.names

# initialize default values for regularization
reg.coefs <- rep(0.5, num.reg.methods)
names(reg.coefs) <- regularization.method.names

# populate the list with default values to be overwritten
target.matrices <- lapply(regularization.method.names, function(x){unreg.matrix}) 
names(target.matrices) <- regularization.method.names

reg.coefs['unregularized'] <- 0 # doesn't matter
target.matrices[['unregularized']] <- unreg.matrix

reg.coefs['ue.zr.full'] <- 1
target.matrices[['ue.zr.full']] <- diag(diag(unreg.matrix))

reg.coefs['eq.zr.half'] <- 0.5
target.matrices[['eq.zr.half']] <- diag(rep(2 * mean.eps2, nrow(unreg.matrix)))

reg.coefs['ue.zr.half'] <- 0.5
target.matrices[['ue.zr.half']] <- diag(diag(unreg.matrix))

reg.coefs['eq.eq.half'] <- 0.5
target.matrices[['eq.eq.half']][unreg.matrix > 0] <- mean.eps2
diag(target.matrices[['eq.eq.half']]) <- 2 * mean.eps2

reg.coefs['ue.eq.half'] <- 0.5
target.matrices[['ue.eq.half']][unreg.matrix > 0] <- mean.eps2
diag(target.matrices[['ue.eq.half']]) <- 2 * mean.eps2

ue.mn.coef.denom <- sum((epsilon.vec - min.eps2) ^ 2) + (1e-14)
ue.mn.coef.numer <- sum((epsilon.vec - mean(epsilon.vec)) ^ 2) + (2e-14)
reg.coefs['ue.mn.half'] <- 0.5
reg.coefs['ue.mn.full'] <- 1
reg.coefs['ue.mn.js1'] <- ue.mn.coef.numer / ue.mn.coef.denom

ue.mn.matrix <- unreg.matrix
ue.mn.matrix[unreg.matrix > 0] <- min.eps2
diag(ue.mn.matrix) <- diag(unreg.matrix)
target.matrices[['ue.mn.half']] <- ue.mn.matrix
target.matrices[['ue.mn.full']] <- ue.mn.matrix
target.matrices[['ue.mn.js1']]  <- ue.mn.matrix

for (curr.reg in regularization.method.names)
{
  regularized.estimates[curr.reg] <-
    cov.to.clonality(target.matrices[[curr.reg]], reg.coefs[curr.reg])
}

curr.clonality.score.estimate <- as.numeric(regularized.estimates['ue.zr.half']) 
# fpc: fixed point iterations. Usually 1 is best

fpc.iter.estimates <- append(fpc.iter.estimates, curr.clonality.score.estimate)
} # for (curr.iter.number in 1: num.iterations)


# jackknife related default values
mixture.clonality <- NA
mixture.estimates <- NA # the vector of different ways to average
d1jkn.covariance <- NA
if (compute.variances.d1jkn) { # use jackknife to estimate variance
mix.regnames <- c('unregularized', 'ue.zr.half', 'ue.mn.half', 'eq.zr.half')
mix.values <- c(regularized.estimates[mix.regnames], simple.precision.clonality)

# include an extra component in the mixture for simple clonality.
unnormalized.d1jkn.cov.sqrt <- 
  matrix(NA, length(mix.regnames) + 1, num.replicates)
for (i in c(1:num.replicates)) {

  curr.subset.internal.parameters <- list(
    replicates = replicates[, -i],
    rep.gram.matrix = rep.gram.matrix[-i, -i],
    use.squared.err.est = use.squared.err.est[-i],
    use.replicate.var.est = use.replicate.var.est[-i],
    compute.variances.d1jkn = F)
  curr.subset.clonality <- infer.clonality(
    read.count.matrix = read.count.matrix[, -i], 
    variance.method = variance.method,
    internal.parameters = curr.subset.internal.parameters)

  curr.d1jk.terms <-
    c((curr.subset.clonality$regularized.estimates)[mix.regnames],
      curr.subset.clonality$simple.precision.clonality)

  unnormalized.d1jkn.cov.sqrt[, i] <- curr.d1jk.terms - mix.values
} # for (i in c(1:num.replicates))

n.var.coef <- (num.replicates - 1) / num.replicates
d1jkn.covariance <- n.var.coef * unnormalized.d1jkn.cov.sqrt %*% t(unnormalized.d1jkn.cov.sqrt)
d1jkn.ev <- eigen(d1jkn.covariance)$values
d1jkn.adjustment <-  # regularize the diagonal of the jackknifed covariance matrix
  (-max(min(d1jkn.ev), 0)) + # negate the smallest eigenvalue, if it is negative
  (1e-4) * max(d1jkn.ev) +   # 1e-4 of the largest eigenvalue
  all(d1jkn.ev == 0) * .Machine$double.eps * 32 # if matrix was identically 0

  diag(d1jkn.covariance) <- diag(d1jkn.covariance) + d1jkn.adjustment

  get.reweighted.clonality.given.covariance <- 
    function(covariance.of.regularized.estimates, mix.values) {
    precision <- ginv(covariance.of.regularized.estimates)
    numerator <- rep(1, nrow(precision)) %*% 
      precision %*% mix.values
    denominator <- rep(1, nrow(precision)) %*% 
      precision %*% rep(1, nrow(precision))
    return(numerator / denominator) # the clonality score
  }
  
  mixture.matrix.weight.clonality <- get.reweighted.clonality.given.covariance(d1jkn.covariance, mix.values)
  mixture.scalar.weight.clonality <- get.reweighted.clonality.given.covariance(diag(diag(d1jkn.covariance)), mix.values)
  mixture.equal.weight.clonality  <- get.reweighted.clonality.given.covariance(diag(rep(1, nrow(d1jkn.covariance))), mix.values)
  mixture.estimates <- c(
    matrix.mix = mixture.matrix.weight.clonality, 
    scalar.mix = mixture.scalar.weight.clonality, 
    equalw.mix = mixture.equal.weight.clonality)
  
  mixture.clonality <- mixture.matrix.weight.clonality
  if (mixture.clonality < min(mix.values) || mixture.clonality > max(mix.values))
  {
    mixture.clonality <- mixture.scalar.weight.clonality
  }
  
} # if (compute.variances.d1jkn)

return.results <- list(
  estimated.abundances = 
    'The estimate.abundances parameter was set to false when infer.clonality was called.',
  internal.parameters = internal.parameters,
  d1jkn.covariance = d1jkn.covariance,
  estimated.squared.errs = epsilon.vec,
  estimated.precisions = inv.eps.vec,
  variance.method = variance.method,
  fpc.iter.estimates = fpc.iter.estimates,
  simple.precision.clonality = simple.precision.clonality, 
  regularized.estimates = regularized.estimates,
  mixture.estimates = mixture.estimates,
  mixture.clonality = mixture.clonality,
  lymphclon.clonality = ifelse(is.na(mixture.clonality), 
    regularized.estimates['ue.zr.half'], 
    mixture.clonality)
    )

if (estimate.abundances) {
  return.results[['estimated.abundances']] = (replicates %*% inv.eps.vec) / sum(inv.eps.vec)
}

return (return.results)
}
