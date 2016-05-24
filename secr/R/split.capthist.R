############################################################################################
## package 'secr'
## split.capthist.R
## last changed 2009 06 11 2009 07 10 2009 10 05 2012 07 26 2012 09 04 2015-10-11
############################################################################################

# split.capthist <- function (x, f, drop = FALSE, prefix='S', bytrap = FALSE, ...) {
#   if (!inherits(x, 'capthist'))
#       stop ("argument to 'split.capthist' should have class 'capthist'")
#   if (inherits(x, 'list'))
#       stop ("split not suitable for multi-session 'capthist'")
#   options(warn=-1)
#
#   f <- as.factor(f)  # retains unused levels
#   if (any(!is.na(as.numeric(levels(f))))) {
#       ## leadingzero added 2012-09-04
#    #   f <- factor(paste (prefix,leadingzero(f),sep=''))
#          levels(f) <- paste (prefix,leadingzero(levels(f)),sep='')
#   }
#   options(warn=0)
#
#   if (bytrap) {
#       if (length(f)!=nrow(traps(x)))
#           stop ("length of f should match number of traps")
#   }
#   else {
#       if (length(f)!=nrow(x))
#           stop ("length of f should match number of rows in capthist")
#   }
#
#   out <- list()
#   for (i in levels(f)) {
#     if (bytrap) {
#         temp <- subset (x, traps = f == i, ...)
#     }
#     else {
#         temp <- subset (x, subset = f == i, ...)
#     }
#     session(temp) <- i
#     if (!drop | (nrow(temp)>0))
#       out[[i]] <- temp
#   }
#   class (out) <- c('list', 'capthist')
#   out
# }
############################################################################################

split.capthist <- function (x, f, drop = FALSE, prefix='S', bytrap = FALSE,
                            byoccasion = FALSE, ...) {
  if (!inherits(x, 'capthist'))
      stop ("argument to 'split.capthist' should have class 'capthist'")
  if (inherits(x, 'list'))
      stop ("split not suitable for multi-session 'capthist'")
  options(warn=-1)

  f <- as.factor(f)  # retains unused levels
  if (any(!is.na(as.numeric(levels(f))))) {
      ## leadingzero added 2012-09-04
   #   f <- factor(paste (prefix,leadingzero(f),sep=''))
         levels(f) <- paste (prefix,leadingzero(levels(f)),sep='')
  }
  options(warn=0)

  if (bytrap) {
      ## 2015-10-11
      ## if (length(f)!= nrow(traps(x)))
      if (length(f)!= ndetector(traps(x)))
          stop ("length of f should match number of detectors")
  }
  else if (byoccasion) {
      if (length(f)!=ncol(x))
          stop ("length of f should match number of columns in capthist")
  }
  else {
      if (length(f)!=nrow(x))
          stop ("length of f should match number of rows in capthist")
  }
  if (bytrap & byoccasion)
      stop("specify only one of bytrap and byoccasion")

  out <- list()
  for (i in levels(f)) {
    if (bytrap) {
        temp <- subset (x, traps = f == i, ...)
    }
    else if (byoccasion) {
        temp <- subset (x, occasions = f == i, ...)
    }
    else {
        temp <- subset (x, subset = f == i, ...)
    }
    session(temp) <- i
    if (!drop | (nrow(temp)>0))
      out[[i]] <- temp
  }
  class (out) <- c('list', 'capthist')
  out
}
############################################################################################

extract.estimates <- function (x, simplify = FALSE) {
## compile a dataframe (simplify = T) or list of data.frames of session-specific real parameter estimates
## from a list of separate secr model fits
   if (!is.list(x) | !inherits(x[[1]], 'secr'))
       stop ("requires list of fitted secr models")
   temp <- lapply(x, predict)
   temp <- lapply(temp, function(x) x[,-1])  ## drop unwanted 'link' column
   temp <- lapply(temp, function(x) {x$Parameter <- row.names(x); x})
   sessions <- names(temp)
   nsessions <- length(temp)
   parnames <- row.names(temp[[1]])
   nrealpar <- nrow(temp[[1]])
   temp2 <- data.frame(abind(temp, along = 1), row.names = NULL, stringsAsFactors=F)
   temp2[,1:4] <- sapply(temp2[,1:4], as.numeric)
   temp2$Session <- rep(sessions, rep(nrealpar, nsessions))
   if (simplify) {
       temp3 <- temp2[order(temp2$Parameter, temp2$Session), c('Parameter','Session','estimate', 'SE.estimate', 'lcl', 'ucl')]
       row.names(temp3) <- NULL
   }
   else {
      temp3 <- split(temp2[order(temp2$Session), c('Session','estimate', 'SE.estimate', 'lcl', 'ucl')],
          temp2$Parameter)
      temp3 <- lapply(temp3, function(x) {row.names(x) <- NULL; x})
   }
   temp3
}
############################################################################################

# extract.estimates(revi)$density


# f1 <- c(rep('A', 30), rep('B', nrow(captdata)-30))
# spltcpt <- split.capthist(captdata, f= f1)
# lapply(spltcpt, summary)
# sessionfits <- lapply(spltcpt, secr.fit)
# extract.estimates (sessionfits)

