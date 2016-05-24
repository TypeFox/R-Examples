#' Distribution evaluation.
#' 
#' @aliases inverse.simpson diversity gini chao1 gini.simpson
#' 
#' @description
#' Functions for evaluating the diversity of species or objects in the given distribution.
#' See the \code{repOverlap} function for working with clonesets and a general interface to
#' all of this functions.
#' 
#' Warning!
#' Functions will check if .data is a distribution of a random variable (sum == 1) or not.
#' To force normalisation and / or to prevent this, set .do.norm to TRUE (do normalisation)
#' or FALSE (don't do normalisation), respectively.
#' 
#' - True diversity, or the effective number of types, refers to the number
#' of equally-abundant types needed for the average proportional abundance
#' of the types to equal that observed in the dataset of interest 
#' where all types may not be equally abundant.
#' 
#' - Inverse Simpson index is the effective number of types that is obtained when
#' the weighted arithmetic mean is used to quantify average 
#' proportional abundance of types in the dataset of interest.
#' 
#' - The Gini coefficient measures the inequality among values
#' of a frequency distribution (for example levels of income). A Gini coefficient of zero
#' expresses perfect equality, where all values are the same (for example, where everyone
#' has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality
#' among values (for example where only one person has all the income).
#' 
#' - The Gini-Simpson index is the probability of interspecific encounter, i.e., probability that two entities
#' represent different types.
#' 
#' - Chao1 estimator is a nonparameteric asymptotic estimator of species richness (number of species in a population).
#' 
#' @usage
#' inverse.simpson(.data, .do.norm = NA, .laplace = 0)
#' 
#' diversity(.data, .q = 5, .do.norm = NA, .laplace = 0)
#' 
#' gini(.data, .do.norm = NA, .laplace = 0)
#' 
#' gini.simpson(.data, .do.norm = NA, .laplace = 0)
#' 
#' chao1(.data)
#' 
#' @param .data Numeric vector of values for proportions or for numbers of individuals.
#' @param .q q-parameter for the Diversity index.
#' @param .do.norm One of the three values - NA, T or F. If NA than check for distrubution (sum(.data) == 1)
#' and normalise if needed with the given laplace correction value. if T then do normalisation and laplace
#' correction. If F than don't do normalisaton and laplace correction.
#' @param .laplace Value for Laplace correction which will be added to every value in the .data.
#' 
#' @return Numeric vector of length 1 with value for all functions except \code{chao1}, which returns 4 values:
#' estimated number of species, standart deviation of this number and two 95% confidence intervals for the species number.
#' 
#' @seealso \link{repOverlap}, \link{entropy}, \link{similarity}
#' 
#' @examples
#' data(twb)
#' # Next two are equal calls:
#' stopifnot(gini(twb[[1]]$Read.count, TRUE, 0) - 0.7609971 < 1e-07)
#' stopifnot(gini(twb[[1]]$Read.proportion, FALSE) - 0.7609971 < 1e-07)
#' stopifnot(chao1(twb[[1]]$Read.count)[1] == 1e+04)
inverse.simpson <- function (.data, .do.norm = NA, .laplace = 0) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  1 / sum(.data ^ 2)
}

diversity <- function (.data, .q = 5, .do.norm = NA, .laplace = 0) {
  .data <- check.distribution(.data, .do.norm, .laplace)
  if (.q == 0) {
    length(.data)
  } else if (.q == 1) {
    1 / prod(.data ^ .data)
  } else if (.q > 1) {
    1 / (sum(.data ^ .q) ^ (1 / (.q - 1)))
  } else {
    NA
  }
}

gini <- function (.data, .do.norm = NA, .laplace = 0) {
  .data <- sort(check.distribution(.data, .do.norm, .laplace, .warn.sum = F))
  n <- length(.data)
  1 / n * (n + 1 - 2 * sum((n + 1 - 1:n) * .data) / sum(.data))
}

gini.simpson <- function (.data, .do.norm = NA, .laplace = 0) {
  1 - 1 / inverse.simpson(.data, .do.norm, .laplace)
}

chao1 <- function (.data) {
  counts <- table(.data)
  e <- NA
  v <- NA
  lo <- NA
  hi <- NA
  n <- sum(.data)
  D <- length(.data)
  f1 <- counts['1']
  f2 <- counts['2']
  # f1 == 0 && f2 == 0
  if (is.na(f1) && is.na(f2)) {
    e <- D
    i <- 1:max(.data)
    i <- i[unique(.data)]
    v <- sum(sapply(i, function(i) sum(.data == i) * (exp(-i) - exp(-2 * i)))) - (sum(sapply(i, function(i) i * exp(-i) * sum(.data == i))))^2/n
    P <- sum(sapply(i, function(i) sum(.data == i) * exp(-i)/D))
    lo <- max(D, D/(1 - P) - qnorm(1 - .05/2) * sqrt(v)/(1 - P))
    hi <- D/(1 - P) + qnorm(1 - .05/2) * sqrt(v)/(1 - P)
  }
  # f1 != 0 && f2 == 0
  else if (is.na(f2)) {
    e <- D + f1 * (f1 - 1) / 2 * (n - 1) / n
    v <- (n-1)/n * f1 * (f1 - 1) /2  + ((n -1)/n)^2 * f1 * (2 * f1 - 1)^2/4 - ((n-1)/n)^2*f1^4/4/e
    t <- e - D
    K <- exp(qnorm(1 - .05/2) * sqrt(log(1 + v/t^2)))
    lo <- D + t/K
    hi <- D + t*K
  }
  # f1 != && f2 != 0
  else {
    const <- (n - 1) / n
    e <- D + f1^2 / (2 * f2) * const
    f12 <- f1 / f2
    v <- f2 * (const * f12^2 / 2 + const^2 * f12^3 + const^2 * f12^4 / 4)
    t <- e - D
    K <- exp(qnorm(.975) * sqrt(log(1 + v/t^2)))
    lo <- D + t/K
    hi <- D + t*K
  }
  c(Estimator = e, SD = sqrt(v), Conf.95.lo = lo, Conf.95.hi = hi)
}

# hill.numbers <- function (.data, .max.q = 6, .min.q = 1) {
#   print(head(.data))
#   print(sum(.data))
#   # .data <- check.distribution(.data)
#   print(head(.data))
#   
#   if (.min.q < 0) { .min.q <- 0 }
#   res <- c()
#   
#   if (.min.q == 0) {
#     res <- length(.data)
#     .min.q <- 1
#   }
#   
#   if (.min.q == 1) {
#     res <- c(res, exp(-sum(.data * log(.data))))
#     .min.q <- 2
#   }
#   
#   if (.max.q >= 2) {
#     for (q in .min.q:.max.q) {
#       res <- c(res, sum(.data ^ q)^(1 / (1 - q)))
#     }
#   }
#   
#   res
# }


#' Diversity evaluation using rarefaction.
#' 
#' @description
#' Sequentially resample the given data with growing sample size the given data and compute mean number of unique clones.
#' For more details on the procedure see "Details".
#' 
#' @param .data Data frame or a list with data frames.
#' @param .step Step's size.
#' @param .quantile Numeric vector of length 2 with quantiles for confidence intervals.
#' @param .extrapolation If N > 0 than perform extrapolation of all samples to the size of the max one +N reads or UMIs.
#' @param .col Column's name from which choose frequency of each clone.
#' @param .verbose if T then print progress bar.
#' 
#' @return
#' Data frame with first column for sizes, second columns for the first quantile,
#' third column for the mean, fourth columns for the second quantile, fifth columns
#' for the name of subject.
#' 
#' @details
#' This subroutine is designed for diversity evaluation of repertoires. On each step it computes a
#' mean unique clones from sample of fixed size using bootstrapping. Unique clones for each sample from bootstrap computed
#' as a number of non-zero elements in a vector from multinomial distribution with input vector of probabilities from the \code{.col} column
#' using function \code{rmultinom} with parameters n = .n, size = i * .step, prob = .data[, .col] (i is an index of current iteration)
#' and choosing for lower and upper bound \code{quantile} bounds of the computed distribution of unique clones.
#' 
#' @seealso \link{vis.rarefaction} \link{rmultinom}
#' 
#' @examples
#' \dontrun{
#' rarefaction(immdata, .col = "Read.count")
#' }
rarefaction <- function (.data, .step = 30000, .quantile = c(.025, .975), .extrapolation = 200000, .col = 'Umi.count', .verbose = T) {
  if (has.class(.data, 'data.frame')) {
    .data <- list(Sample = .data)
  }
  
  # multinom
#   .alpha <- function (n, Xi, m) {
#     k <- Xi
#     if (k <= n - m) {
#       prod((n - k):(n - m - k + 1) / n:(n - m + 1))
#     } else {
#       0
#     }
#   }

  # poisson
  .alpha <- function (n, Xi, m) {
    k <- Xi
    return((1 - m / n)^Xi)
  }
  
  if (.verbose) {
    pb <- set.pb(sum(sapply(1:length(.data), function (i) {
      bc.vec <- .data[[i]][, .col]
      bc.sum <- sum(.data[[i]][, .col])
      sizes <- seq(.step, bc.sum, .step)
      if (sizes[length(sizes)] != bc.sum) {
        sizes <- c(sizes, bc.sum)
      }
      length(sizes)
    } )))
  }
  
  muc.list <- lapply(1:length(.data), function (i) {
    Sobs <- nrow(.data[[i]])
    bc.vec <- .data[[i]][, .col]
    Sest <- chao1(bc.vec)
    n <- sum(bc.vec)
    sizes <- seq(.step, n, .step)
    if (sizes[length(sizes)] != n) {
      sizes <- c(sizes, n)
    }
    counts <- table(bc.vec)
    muc.res <- t(sapply(sizes, function (sz) {
      freqs <- as.numeric(names(counts))
      
      # multinom
      alphas <- sapply(freqs, function (k) .alpha(n, k, sz))
#       Sind <- Sobs - sum(sapply(freqs, function (k) .alpha(n, k, sz) * counts[as.character(freqs)]))
#       SD <- sqrt(sum(sapply(freqs, function (k) (1 - .alpha(n, k, sz))^2 * counts[as.character(freqs)])) - Sind^2/Sest[1])
      
      # poisson
      Sind <- sum(sapply(1:length(freqs), function (k) (1 - alphas[k]) * counts[k]))
      if (Sest[1] == Sobs) {
        SD <- 0
      } else {
        SD <- sqrt(sum(sapply(1:length(freqs), function (k) (1 - alphas[k])^2 * counts[k])) - Sind^2/Sest[1])
      }
      t <- Sind - Sobs
      K <- exp(qnorm(.975) * sqrt(log(1 + (SD / t)^2)))
      lo <- Sobs + t*K
      hi <- Sobs + t/K
      res <- c(sz, Sind, Sind, Sind)
      names(res) <- c('Size', paste0('Q', .quantile[1]), 'Mean', paste0('Q', .quantile[2]))
      if (.verbose) add.pb(pb)
      res
    }))

    if (.extrapolation > 0) {
      sizes <- seq(sum(.data[[i]][, .col]), .extrapolation + max(sapply(.data, function (x) sum(x[, .col]))), .step)
      if (length(sizes) != 1) {        
        ex.res <- t(sapply(sizes, function (sz) {
          f0 <- Sest[1] - Sobs
          f1 <- counts['1']
          if (is.na(f1) || f0 == 0) {
            Sind <- Sobs
          } else {
            Sind <- Sobs + f0 * (1 - exp(-(sz - n)/n * f1 / f0))
          }
          res <- c(sz, Sind, Sind, Sind)
          names(res) <- c('Size', paste0('Q', .quantile[1]), 'Mean', paste0('Q', .quantile[2]))
          if (.verbose) add.pb(pb)
          res
        }))
        df1 <- data.frame(muc.res, People = names(.data)[i], Type = 'interpolation', stringsAsFactors = F)
        df2 <- data.frame(ex.res, People = names(.data)[i], Type = 'extrapolation', stringsAsFactors = F)
        rbind(df1, df2)
      } else {
        df1 <- data.frame(muc.res, People = names(.data)[i], Type = 'interpolation', stringsAsFactors = F)
      }
    } else {
      data.frame(muc.res, People = names(.data)[i], stringsAsFactors = F)
    }
  })
  if (.verbose) close(pb)
  
  do.call(rbind, muc.list)
}