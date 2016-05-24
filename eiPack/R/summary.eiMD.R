summary.eiMD <- function(object, short = TRUE, ...) {
  "%w/o%" <- function(x,y) x[!x %in% y]
  get2 <- function(x) x[2]
  getm1 <- function(x) x[2:length(x)]

  if (is.mcmc(object$draws$Cell.counts)) { 
    tnames <- strsplit(colnames(object$draws$Cell.counts), "ccount.")
    idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
    idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE,
                                        nrow = length(idx), ncol = length(idx[[1]]))))
    idx <- lapply(idx, as.character)
    idx <- lapply(idx, unique)
  } else {
    idx <- dimnames(object$draws$Cell.counts)[1:2]
  }
  names(idx) <- c("rows", "columns")
  cnames <- apply(expand.grid(idx), 1, paste, collapse = ".")
  
  cells <- prod(sapply(idx, length))
  for (ii in names(object$acc.ratios) %w/o% c("beta.acc")) {
    ll <- length(object$acc.ratios[[ii]])
    if (ll == length(idx[[1]])) {
      names(object$acc.ratios[[ii]]) <- idx[[1]]
    }
    else if (ll < cells) {
      cc <- ll / length(idx[[1]])
      object$acc.ratios[[ii]] <- matrix(object$acc.ratios[[ii]],
                                        nrow = length(idx[[1]]),
                                        ncol = cc,
                                        dimnames = list(idx[[1]], idx[[2]][1:cc]))
    }
    else if (ll == cells) {
      object$acc.ratios[[ii]] <- matrix(object$acc.ratios[[ii]],
                                        nrow = length(idx[[1]]),
                                        ncol = length(idx[[2]]),
                                        dimnames = idx[1:2])
    }
  }
  if (short) {
    # old code created r by c array, not r by (c-1) by l array
    #tmp <- array(object$acc.ratios$beta.acc,
    #               dim = sapply(idx, length),
    #               dimnames = idx)
    rr <- length(idx[[1]])#
    cc <- length(idx[[2]]) - 1#
    ll <- length(object$acc.ratios$beta.acc)/(rr*cc)#
    tmp <- array(object$acc.ratios$beta.acc,#
                   dim = c(rr,cc,ll))#
    object$acc.ratios$beta.acc <- apply(tmp, c(1,2), mean)
    dimnames(object$acc.ratios$beta.acc) <- list(idx[[1]], idx[[2]][1:cc])#
  }  else  {
    #old code filled matrix in wrong direction
    bacc <- object$acc.ratios$beta.acc
    rr <- length(idx[[1]])#
    cc <- length(idx[[2]]) - 1#
    ll <- length(object$acc.ratios$beta.acc)/(rr*cc)#
    object$acc.ratios$beta.acc <-  matrix(bacc,
                                          nrow = ll,
                                          ncol = rr*cc,
                                          dimnames =
                                          list(as.character(1:ll),cnames[1:(rr*cc)]), byrow=TRUE)#
  }

  for (ii in names(object$draws) %w/o% c("Beta")) {
    aa <- object$draws[[ii]]
    if (!is.mcmc(aa)) {
      if (length(dim(aa)) > 2) {
        nc <- prod(dim(aa)[1:2])
        aa <- matrix(c(aa), nrow = dim(aa)[3], ncol = nc,
                     byrow = TRUE, dimnames = list(NULL, cnames[1:nc]))
      }
      else
        aa <- t(aa)
    }
    object$draws[[ii]] <- cbind(apply(aa, 2, mean), apply(aa, 2, sd),
                                t(apply(aa, 2, quantile, c(0.025,0.975))))
    colnames(object$draws[[ii]])[1:2] <- c("Mean", "Std. Error")
    if (ncol(aa) == length(idx[[1]]))
      rownames(object$draws[[ii]]) <- idx[[1]]
    else if (ncol(aa) <= cells)
      rownames(object$draws[[ii]]) <- cnames[1:ncol(aa)]
  }
    
  object$short <- short
  class(object) <- "eiMDsum"
  object
}
