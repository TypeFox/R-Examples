predict_clusterwise <-
function(object, newdata, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  probs <- rcpp_calculate_haplotype_probabilities_clusterwise(newdata, object$y, object$disclap_parameters, object$tau)
  return(probs)
}

happrobsum_within <- function(object, alleles, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  if (ncol(alleles) != ncol(object$y)) stop("Please specify alleles for exactly the number of loci the model was fitted for")
  prob_sum <- rcpp_calculate_haplotype_probabilities_sum(alleles, object$y, object$disclap_parameters, object$tau)
  return(prob_sum)
}

happrobsum_within_between <- function(objects, alleles, ...) {
  if (!is.list(objects) || length(objects) <= 1L) {
    stop("Please specify at least two disclapmixfit fits")
  }

  for (o in objects) {
    if (!is(o, "disclapmixfit")) stop("all objects must be a disclapmixfit")
    if (ncol(alleles) != ncol(o$y)) stop("Please specify alleles for exactly the number of loci that all the objects was fitted for")
  }

  fits <- lapply(objects, function(o) {
    return(list(y = o$y, p = o$disclap_parameters, tau = o$tau))
  })

  res <- rcpp_calculate_haplotype_probabilities_sum_CLASS(fits, alleles)

  class(res) <- c("happrobsum_within_between", class(res))
  
  return(res)
}

happrobsum_within_between_cache <- function(objects, alleles, normalise = FALSE, ...) {
  if (!is.list(objects) || length(objects) <= 1L) {
    stop("Please specify at least two disclapmixfit fits")
  }

  for (o in objects) {
    if (!is(o, "disclapmixfit")) stop("all objects must be a disclapmixfit")
    if (ncol(alleles) != ncol(o$y)) stop("Please specify alleles for exactly the number of loci that all the objects was fitted for")
  }

  fits <- lapply(objects, function(o) {
    return(list(y = o$y, p = o$disclap_parameters, tau = o$tau))
  })

  res <- rcpp_calculate_haplotype_probabilities_sum_CLASS_Cache(fits, alleles)
  
  class(res) <- c("happrobsum_within_between", class(res))

  if (!is.null(normalise) && !is.na(normalise) && is.logical(normalise) && length(normalise) == 1L && normalise == TRUE) {
    res <- happrobsum_normalise(res)
  }
  
  return(res)
}

happrobsum_normalise <- function(hapsum, change_sums = TRUE, ...) {
  if (!is(hapsum, "happrobsum_within_between")) stop("object must be a happrobsum_within_between")

  hapsum$match_within <- hapsum$match_within / (hapsum$hap_sum * hapsum$hap_sum)
  
  n <- length(hapsum$match_within)
  for (i1 in 1L:(n - 1L)) {
    for (i2 in (i1 + 1L):n) {
      hapsum$match_between[i1, i2] <- hapsum$match_between[i1, i2] / (hapsum$hap_sum[i1] * hapsum$hap_sum[i2])
    }    
  }
  
  if (!is.null(change_sums) && !is.na(change_sums) && is.logical(change_sums) && length(change_sums) == 1L && change_sums == TRUE) {
    hapsum$hap_sum <- rep(1, length(hapsum$hap_sum))
  }
  
  return(hapsum)    
}

happrobsum_within_between_binomial <- function(dbs, ...) {
  if (length(dbs) <= 1L) {
    stop("At least two databases are required")
  }
  
  for (db in dbs) {
    if (!is.matrix(db) | !is.integer(db)) {
      stop("haplotypes must be an integer matrix")
    }
  }

  loci <- ncol(dbs[[1L]])
  
  for (db in dbs) {
    if (ncol(db) != loci) {
      stop("All databases must have the same number of loci")
    }
  }
  
  res <- rcpp_calculate_haplotype_probabilities_sum_binomial(dbs)
  class(res) <- c("happrobsum_within_between", class(res))  
  return(res)
}

# Expect count column to be the last column
happrobsum_within_between_binomial_compact_dbs <- function(compact_dbs, ...) {
  if (length(compact_dbs) <= 1L) {
    stop("At least two databases are required")
  }
  
  for (db in compact_dbs) {
    if (!is.matrix(db) | !is.integer(db)) {
      stop("haplotypes must be an integer matrix")
    }
  }
  
  for (db in compact_dbs) {
    if (ncol(db) != ncol(compact_dbs[[1L]])) {
      stop("All databases must have the same number of loci")
    }
  }
  
  counts <- lapply(compact_dbs, function(db) as.integer(db[, ncol(db)]))
  dbs <- lapply(compact_dbs, function(db) (db[, 1L:(ncol(db)-1L), drop = FALSE]))
  
  res <- rcpp_calculate_haplotype_probabilities_sum_binomial_compact_dbs(dbs, counts)
  class(res) <- c("happrobsum_within_between", class(res))  
  return(res)
}

# Normalises
happrobsum_within_between_normalised <- function(objects, dbs_to_eval, ...) {
  if (!is.list(objects) || length(objects) <= 1L) {
    stop("Please specify at least two disclapmixfit fits")
  }
  
  for (db in dbs_to_eval) {  
    if (!is.matrix(db) | !is.integer(db)) {
      stop("All db in dbs_to_eval must be an integer matrix")
    }    
  }
  
  loci <- ncol(dbs_to_eval[[1L]])

  for (o in objects) {
    if (!is(o, "disclapmixfit")) stop("all objects must be a disclapmixfit")
    if (ncol(o$y) != loci) stop("Please specify dbs_to_eval for exactly the number of loci that all the objects was fitted for")
  }

  fits <- lapply(objects, function(o) {
    return(list(y = o$y, p = o$disclap_parameters, tau = o$tau))
  })

  # FIXME: Faster? C++ hash?
  db_unique <- do.call(rbind, dbs_to_eval)
  db_unique <- db_unique[!duplicated(db_unique), ]

  res <- rcpp_hapsums_disclap_normalised(fits, db_unique)

  class(res) <- c("happrobsum_within_between", class(res))
  
  return(res)
}

theta_weir <- function(hapsum, ...) {
  if (!is(hapsum, "happrobsum_within_between")) stop("object must be a happrobsum_within_between")

  mi <- hapsum$match_within
  mij <- hapsum$match_between

  mW <- mean(mi)
  mA <- mean(mij[upper.tri(mij)])
  
  r <- length(mi)
  disclap_theta_approx <- (mW - mA) / (1 - mA)
  disclap_theta <- (((r-1)/r) * disclap_theta_approx) / (1 - (1/r)*disclap_theta_approx)
  
  #return(weir_theta)
  return(list(theta = disclap_theta, theta_approx = disclap_theta_approx))
}

#match_prob_quantities <- function(object, alleles, ...) {
#  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
#  if (ncol(alleles) != ncol(object$y)) stop("Please specify alleles for exactly the number of loci the model was fitted for")
#  ms <- rcpp_match_quantities(alleles, object$y, object$disclap_parameters, object$tau)
#  return(ms)
#}

#  * New function to calculate theta
#estimate_theta <- function(object, alleles, ...) {
#  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
#  if (ncol(alleles) != ncol(object$y)) stop("Please specify alleles for exactly the number of loci the model was fitted for")
#  if (nrow(object$y) <= 1L) stop("Need a fit with at least two (2) clusters")
#  
##  #prob_sums <- happrobsum(object, alleles)
#  ms <- match_prob_quantities(object, alleles)
#
#  mis <- unlist(lapply(1L:nrow(object$y), function(j) {
#    rcpp_calculate_haplotype_probabilities_sum(alleles, 
#      object$y[j, , drop = FALSE], object$disclap_parameters[j, , drop = FALSE], 1)[2]
#  }))
#  
#  mW <- mean(mis)
#  mA <- mean(ms[upper.tri(ms)])
#  betaW <- (mW - mA) / (1 - mA)
#  
#  return(betaW)
##  #Mw <- betaW + (1-betaW)*as.numeric(prob_sums[2L])
##  #theta <- (Mw - prob_sums[2L]) / (1 - prob_sums[2L])
##
##  #return(list(theta = betaW, prob_sums = prob_sums))
##
##  #return(list(theta = theta, prob_sums = prob_sums))
##}

#  * New helper function: convert_to_compact_db
convert_to_compact_db <- function(x) { 
  #if (!is.matrix(x) || !is.data.frame(x)) {
  #  stop("x must be a matrix or data.frame")
  #}

  if (nrow(x) <= 1L) {
    ret <- data.frame(x, Ndb = 1L)
    ret$ind <- list("1" = 1L)
    return(ret)
  }

  #if (is.matrix(x) && (dim(x)[2L] == 1L)) {
  #  x <- as.vector(x) 
  #}
 
  x_ord <- do.call(order, as.data.frame(x))
   
  #if (is.vector(x)) {
  #  same_as_previous <- x[tail(x_ord, -1L)] == x[head(x_ord, -1L)]
  #} else {
    same_as_previous <- rowSums(x[tail(x_ord, -1L), , drop = FALSE] != x[head(x_ord, -1L), , drop = FALSE]) == 0L
  #}	
 
  indices <- split(x_ord, cumsum(c(TRUE, !same_as_previous)))
 
  #if (is.vector(x)) {
  #  x <- x[sapply(indices, function (x) x[[1L]]), drop = FALSE]
  #} else {
    x <- x[sapply(indices, function (x) x[[1L]]), , drop = FALSE]
  #}
  
  return(data.frame(x, Ndb = sapply (indices, length), ind = I(indices)))
}

#  * New helper function: find_haplotype_in_matrix
find_haplotype_in_matrix <- function(mat, haplotype) {
  if (!is.matrix(mat) || !is.integer(mat) || !is.integer(haplotype)) {
    stop("mat must be an integer matrix and haplotype an integer vector")
  }
  if (length(haplotype) != ncol(mat)) stop("Wrong dimensions")

  i <- rcpp_find_haplotype_in_matrix(mat, haplotype)

  if (i <= 0L) {
    i <- NULL
  }

  return(i)
}



haplotype_diversity <- function(object, nsim = 1e4L) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (is.null(nsim) || length(nsim) != 1L || !is.integer(nsim) || nsim <= 0L) {
    stop("nsim must be >= 1L (note the L postfix for integer)")
  }
  
  get_db_counts <- function(x) { 
    order_x <- do.call(order, as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order_x, -1),] != x[head(order_x, -1),]) == 0
    indices <- split(order_x, cumsum(c(TRUE, !equal.to.previous)))
    Ns <- unlist(lapply(indices, length))
    return(Ns)
  }
  
  db <- simulate.disclapmixfit(object, nsim = nsim)
  Ns <- get_db_counts(db)
  freqs <- Ns / nsim
  #D <- 1 - sum(freqs^2)
  D <- (nsim / (nsim - 1)) * (1 - sum(freqs^2))
  
  return(D)
}


