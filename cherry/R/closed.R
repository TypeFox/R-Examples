setClass("closure",
  representation(
    hypotheses = "character",   # names of hypotheses
    alpha = "numeric",          # stores chosen alpha for testing
    adjusted = "numeric",       # stores adjusted p-values
                                #  (has value 1 if adjusted p > alpha)
                                #  (numeric(0) if no adjusted p-values calculated)
    max.alpha = "numeric",      # highest alpha at which adjusted p-values are calculated
                                #  value NA if no adjusted p-values present
    defining = "integer"        # defining hypotheses at chosen alpha
                                #  NA if alpha = NA
  )
)

###
# This function performs closed testing
###
closed <- function(test, hypotheses, alpha = 0.05, adjust = FALSE) {

  # reverse argument order if mistaken
  if (is(test, "character") && is(hypotheses, "function")) {
    tmp <- test
    test <- hypotheses
    hypotheses <- tmp
  }

  # alpha=NA means adjusted p-values
  if ((!missing(alpha)) && is.na(alpha)) {
    alpha <- 1
    adjust <- TRUE
  }
  
  # default of alpha is 1 if adjust = TRUE
  if (adjust && missing(alpha))
    alpha <- 1

  # preparation of closure
  N <- length(hypotheses)
  Nmax <- log2(.Machine$integer.max+1)
  if (N > Nmax)
    stop("no more than ", Nmax, " hypotheses supported in full closed testing.\n Use a shortcut-based test.")
  closure <- 1:(2^N-1)
  base <- 2^(1:N-1)

  # finds offspring hypotheses of a hypothesis (NB including self)
  offspring <- function(x) {
    res <- bitAnd(x, closure)
    res[res != 0]
  }

  # sort the closure to decreasing number of participating hypotheses
  lengths <- rowSums(sapply(base, function(bs) bitAnd(closure, bs) != 0))
  closure <- closure[sort.list(lengths, decreasing = TRUE)]

  # perform closed testing (adjusted p variant)
  if (adjust) {
    adjusted <- numeric(2^N-1)
    for (i in closure) {
      if (adjusted[i] < 1) {
        localtest <- test(hypotheses[.bit2boolean(i, N)])
        if (localtest > alpha)    # p-values over alpha are set to 1 to save calculations
          localtest <- 1
        if (localtest > adjusted[i]) {
          offs <- offspring(i)
          adjusted[offs] <- pmax(adjusted[offs], localtest)
        }
      }
    }
    def <- as.integer(NA)
  } else {
  # perform closed testing (rejection variant)
    rejected <- !logical(2^N-1)
    for (i in closure) {
      if (rejected[i]) {
        localtest <- test(hypotheses[.bit2boolean(i, N)])
        if (localtest > alpha) {
          offs <- offspring(i)
          rejected[offs] <- FALSE
        }
      }
    }
    adjusted <- numeric(0)
    def <- .defining(which(rejected), N)
  }
  
  # reduce: find defining rejections

  # return
  max.alpha <- ifelse(adjust, alpha, as.numeric(NA))
  alpha <- ifelse(adjust, as.numeric(NA), alpha)
  out <- new("closure", hypotheses = hypotheses,
                        alpha = alpha,
                        adjusted = adjusted,
                        max.alpha = max.alpha,
                        defining = def)
  return(out)
}

###
# show method for closure object
###
setMethod("show", "closure", function(object) {
  cat("Closed testing result on", length(object@hypotheses), "elementary hypotheses.\n")
  if (!is.na(object@alpha)) {
    cat("At confidence level ", 1-object@alpha, ": ", sep="")
    res <- pick(object, object@hypotheses, silent=TRUE)
    cat("False hypotheses >= ", res, "; ", sep="")
    cat("True hypotheses <= ", length(object@hypotheses) - res, ".\n", sep="")
  }
  object
})

###
# slot accession methods for closure object
###
setGeneric("hypotheses", function(object, ...) standardGeneric("hypotheses"))
setMethod("hypotheses", "closure", function(object, ...) {
  object@hypotheses
})

setGeneric("alpha", function(object, ...) standardGeneric("alpha"))
setMethod("alpha", "closure", function(object, ...) {
  object@alpha
})

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))
setMethod("alpha<-", "closure", function(object, value) {
  if (length(object@adjusted)==0)
    stop("Only closure objects with adjusted p-values can reset alpha.")
  if (is.null(value) || is.na(value)) {
    object@alpha <- as.numeric(NA)
    object@defining <- integer(0)
  } else {
    if (value > object@max.alpha)
      stop("Adjusted p-values only available up to alpha = ", object@max.alpha)
    object@alpha <- value
    rejected <- (object@adjusted <= value)
    object@defining <- .defining(which(rejected), length(object@hypotheses))
  }
  object
})

setGeneric("defining", function(object, alpha, ...) standardGeneric("defining"))
setMethod("defining", "closure", function(object, alpha, ...) {
  if (missing(alpha)) {
    if (is.na(object@alpha)) {
      if (object@max.alpha < 1)
        alpha <- object@max.alpha
      else
        stop("Please specify alpha.")
    }
  } else
    alpha(object) <- alpha
  .num2names(object@defining, object@hypotheses)
})


setGeneric("shortlist", function(object, alpha, ...) standardGeneric("shortlist"))
setMethod("shortlist", "closure", function(object, alpha, ...) {
  if (missing(alpha)) {
    if (is.na(object@alpha)) {
      if (object@max.alpha < 1)
        alpha <- object@max.alpha
      else
        stop("Please specify alpha.")
    }
  } else
    alpha(object) <- alpha
  .num2names(.shortlist(object), object@hypotheses)
})


adjusted <- function(closure, reject, n = 0) {
  
  # check if adjusted p available
  if (is.na(closure@max.alpha))
    stop("no adjusted p-values in this closure object")
  
  # transform to number
  N <- length(closure@hypotheses)
  reject <- which(closure@hypotheses %in% reject)
  M <- length(reject)
  reject <- sum(2^(reject-1))
  
  if (n>M)
      stop("value of n larger than number of rejected hypotheses")
  
  clos <- 1:(2^N-1)
  base <- 2^(1:N-1)
  interest <- unique(bitAnd(reject, clos)) 
  lengths <- lapply(base, function(bs) bitAnd(interest, bs) != 0)
  lengths <- rowSums(do.call(cbind, lengths))
  interest <- interest[lengths==M-n]
  return(max(closure@adjusted[interest]))
}


pick <- function(closure, reject, alpha, silent = FALSE, plot = FALSE) {

  # preparation: create closure
  N <- length(closure@hypotheses)
  if (missing(reject)) reject <- closure@hypotheses 
  reject <- which(closure@hypotheses %in% reject)
  M <- length(reject)
  clos <- 1:(2^N-1)
  base <- 2^(1:N-1)

  # the part of the closure that actually matters for this pick
  interest <- unique(bitAnd(sum(2^(reject-1)), clos))
  interest <- interest[interest>0]

  # should adjusted p variant be used?
  if (!missing(alpha))
    alpha(closure) <- alpha

  if (!is.na(closure@alpha)) {
    # fixed alpha variant
    isAncestor <- function(x,y) { # is x an ancestor of y?
      bitOr(x,y) == x
    }
    interest <- interest[!apply(outer(interest, closure@defining, isAncestor), 1, any)]
    lengths <- lapply(base, function(bs) bitAnd(interest, bs) != 0)
    lengths <- rowSums(do.call(cbind, lengths))
    out <- max(lengths,0)
    if (!silent) {
      cat(length(reject), " hypotheses selected. At confidence level ", 1-closure@alpha, ":\n", sep="")
      cat("False null-hypotheses >= ", length(reject)-out, "; ", sep="")
      cat("True null-hypotheses  <= ", out, ".\n", sep="")
      return(invisible(length(reject)-out))
    } else
      return(length(reject)-out)
  } else {
    # adjusted p variant
    lengths <- lapply(base, function(bs) bitAnd(interest, bs) != 0)
    lengths <- rowSums(do.call(cbind, lengths))
    cumulative <- tapply(closure@adjusted[interest], lengths, max)
    diffs <- rev(diff(rev(c("0"=1, cumulative, 0))))
    out <- data.frame(alpha = rev(cumulative),
                        confidence = 1-rev(cumulative),
                        "true<=" = M:1-1,
                        "false>=" = 1:M, check.names=F, row.names=1:M)
    if (plot) {
      bp <- barplot(diffs, xlab="True hypotheses", ylab = "Confidence probability mass function", ylim=c(0,max(diffs)*1.1))
      mids <- (bp[-1] + bp[-length(bp)])/2
      text(mids, max(diffs), round(1-cumsum(diffs)[-length(diffs)],3), pos=3)
    }
    out
  }
}






###
# helper functions
###

# converts integer to set notation
.num2names <- function(rejected, vars) {
  N <- length(vars)
  bools <- lapply(rejected, .bit2boolean, N=N)
  lapply(bools, function(b) vars[b])
}

# converts from integer to boolean (as binary)
.bit2boolean <- function(x, N) {
  base <- 2^(1:N-1)
  bitAnd(x, base) != 0
}

# gets the defining rejections from all rejections
.defining <- function(rejected, N) {
  closure <- 1:(2^N-1)
  ancestors <- function(x) {
    bitOr(x, closure)
  }
  isdone <- integer(0)
  todo <- rejected
  while (length(todo) > 0) {
    isdone <- c(setdiff(isdone, ancestors(todo[1])), todo[1])
    todo <- setdiff(todo, ancestors(todo[1]))
  }
  isdone
}

# reverses defining hypotheses from (A or B or ...) and (C or D or ...) and ...
# to (A and C and ...) or (B and C and ...) or ...
# result is still in bit-form, to be transformed with .num2names()
# reverses defining hypotheses from (A or B or ...) and (C or D or ...) and ...
# to (A and C and ...) or (B and C and ...) or ...
# result is still in bit-form, to be transformed with .num2names()
.shortlist <- function(cl) {
  N <- length(cl@hypotheses)
  M <- length(cl@defining)
  base <- 2^(1:N-1)
  res <- 0
  for (i in 1:M) {
    whichs <- which(bitAnd(cl@defining[i], base) != 0)
    comb <- outer(res, 2^(whichs-1), bitOr)
    res <- unique(c(comb))
  }
  isAncestor <- function(x,y) { # is x an ancestor of y?
    bitOr(x,y) == x
  }
  ancs <- outer(res, res, isAncestor)
  diag(ancs) <- FALSE
  res[!apply(ancs, 1, any)]
}

.shortlist_old <- function(cl) {
  N <- length(cl@hypotheses)
  M <- length(cl@defining)
  base <- 2^(1:N-1)
  lengths <- sapply(cl@defining, function(x) sum(bitAnd(x, base) != 0))
  total <- prod(lengths)
  whichs <- lapply(cl@defining, function(x) which(bitAnd(x, base) != 0))
  ands <- unique(sapply(1:total, function(k) {
    ix <- (k-1) %/% (total/cumprod(lengths)) %% lengths + 1
    choice <- sapply(1:M, function(i) whichs[[i]][ix[i]])
    sum(2^(unique(choice)-1))
  }))
  isAncestor <- function(x,y) { # is x an ancestor of y?
    bitOr(x,y) == x
  }
  ancs <- outer(ands, ands, isAncestor)
  diag(ancs) <- FALSE
  ands[!apply(ancs, 1, any)]
}


