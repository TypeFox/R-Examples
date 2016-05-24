"bounds" <-
function(y, type = "MRC", ...)
{
  switch(type,
    MRC = bounds.MRC(y, ...),
    stop("type ", type, " not known!")
  )
}

# compute two-sided bounds based on MRC, allows to combine these for a subset of right indices (candidates)
"bounds.MRC" <-
function(y, q, alpha = 0.05, r = ceiling(50 / min(alpha, 1 - alpha)), lengths = if(family == "gaussKern") 2^(floor(log2(length(y))):ceiling(log2(length(param$kern)))) else 2^(floor(log2(length(y))):0), penalty = c("none", "len", "var", "sqrt"), name = if(family == "gaussKern") ".MRC.ktable" else ".MRC.table", pos = .GlobalEnv, family = c("gauss", "gaussvar", "poisson", "binomial","gaussKern"), param = NULL, subset, max.iter = 1e2, eps = 1e-3)
{
  family <- match.arg(family)
  if(family == "binomial") {
    if(is.null(param)) stop("param must be size of binomial distribution")
    if(any(y < 0) | any(y > param)) stop("only values in {0, ..., param} allowed for family binomial")
    y <- as.integer(y)
  }
  if(family == "poisson") {
    if(any(y < 0)) stop("only values in the natural numbers allowed for family poisson")
    y <- as.integer(y)
  }
  if(family == "gaussKern" & is.null(param)) stop("param must be dfilter")
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
    pen <- switch(penalty,
      none = rep(0, length(lengths)),
      len = log(lengths / length(y)),
      sqrt = sqrt(2 * (1 - log(lengths / length(y)))),
      var = log(lengths / length(y)),
      stop("unknown penalty")
    )
    if(family == "binomial" & penalty == "var") {
      lengths <- lengths[lengths * param > 1] # length * size <= 1 always results in trivial bounds, so can be ignored
    }
  } else {
    stop("arbitrary penalty not implemented")
  }
  algo <- switch(family,
    gaussKern = "gauss",
    family
  )
  if(family == "gauss" & is.null(param)) {
    stdev <- sdrobnorm(y)
  } else {
    stdev <- param
  }
  if(family == "gaussKern") {
    kl <- length(param$kern)
    lengths <- pmin(lengths, length(y) - kl)
    stdev <- sdrobnorm(y, lag = kl)
  }
  # sort lengths increasingly
  pen <- pen[order(lengths)]
  lengths <- sort(lengths)
  # compute quantile
  if(missing(q)) {
    if(is.null(r)) stop("q or r need to be specified!")
    q <- if(family == "gaussKern") {
      kMRC.quant(1 - alpha, length(y), r, param$kern, lengths, name, pos)
    } else {
      MRC.quant(1 - alpha, length(y), r, lengths, penalty, name, pos)
    }
  }
  # compute cumulative sums
  if(algo == "gaussvar") {
    S <- t(MRCoeff(y^2, lengths = lengths, norm = 1, signed = TRUE))
    keep <- is.finite(S)
    len <- rep(lengths, length(y))[keep]
    meanY2 <- S[keep] / len
    if(penalty=="sqrt") {
      rhs <- ( q + pen )^2 / lengths  + 1 # = z - log(z), z = meanY2 / sigma2
    }else
      rhs <- ( q - pen ) / lengths * 2 + 1 # = z - log(z), z = meanY2 / sigma2
  } else {
    S <- t(MRCoeff(y, lengths = lengths, norm = 1, signed = TRUE))
    keep <- is.finite(S)
    S <- S[keep]
    len <- rep(lengths, length(y))[keep]
    pen <- rep(pen, length(y))[keep]
  }
  # compute bounds and corresponding indices
  bounds <- data.frame(li = as.integer( rep(1:length(y), each = length(lengths))[keep] ))
  bounds$ri <- as.integer( bounds$li + len - 1 )
  if(algo == "gauss") {
    if (penalty=="sqrt"){
      b <- (q + pen) / sqrt(len) * stdev
    } else{
      b <- sqrt(2 * ( q - pen ) / len) * stdev
    }
    bounds$lower <- as.numeric( S / len - b )
    bounds$upper <- as.numeric( S / len + b )
  } else if(algo == "gaussvar") {
    # solve rhs = z - log(z), z = meanY2 / sigma2, set z = exp(r), i.e. exp(r) - r - rhs = 0
    a <- -rhs
    b <- -1
    #   c <- 1
    upp <- rhs
    for(i in 1:max.iter) {
      cer <- exp(upp)
      cer.inv <- 1 / cer
      if(all(abs( cer + b * upp + a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      p2 <- b * cer.inv + 1
      # always disregard the left half of the parabola, go linear on the left side
      upp <- ifelse(upp >= 0, sqrt(pmax(p2^2 - ( a + b * upp ) * cer.inv * 2 - 2, 0)) - p2 + upp, -(a + (1 - upp) * cer) / (b + cer) )
    }
    bounds$lower <- meanY2 / rep(exp(upp), length(y))[keep]
    low <- exp(-rhs)
    for(i in 1:max.iter) {
      cer <- exp(low)
      cer.inv <- 1 / cer
      if(all(abs( cer + b * low + a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      p2 <- b * cer.inv + 1
      # always disregard the left half of the parabola, go linear on the left side
      low <- ifelse(low >= 0, -sqrt(pmax(p2^2 - ( a + b * low ) * cer.inv * 2 - 2, 0)) - p2 + low, -(a + (1 - low) * cer) / (b + cer) )
    }
    bounds$upper <- meanY2 / rep(exp(low), length(y))[keep]
  } else if(algo == "poisson") {
    # compute cumulative sums
    S <- round(S)
    nonzeros <- S != 0
    nontriv <- if(penalty != "var") nonzeros else S != 1 # solution is trivial for S = 0 and no penalty, or S = 1 and log-penalty
    # compute bounds and corresponding indices, we use esolve above to solve for the log of length * mu
    Snt <- S[nontriv] 
    if(penalty == "var") a <- -q - log(sum(y))
    if(penalty == "sqrt") a <- -(q+pen)^ 2 / 2
    if(penalty=="len" || penalty=="none") a <- -q + pen
    
    a <- a + ifelse(nonzeros, -S * ( 1 - log(S) ), 0)
    bounds$lower <- rep(0, length(S))
#     bounds$upper <- ifelse(nontriv, log(ifelse(nonzeros, S + q, q)), -a / len)
    bounds$upper <- ifelse(nontriv, log(ifelse(nonzeros, max(S + q,0.001), max(q,0.001))), -a / len) # allow for negative q , see previous line
    uppnt <- bounds$upper[nontriv]
    a <- a[nontriv]
    b <- -ifelse(nonzeros, S, ifelse(penalty == "var", 0, 1)) + ifelse(penalty == "var", 1, 0)
    b <- b[nontriv]
    #   c <- 1
    for(i in 1:max.iter) {
      cer <- exp(uppnt)
      cer.inv <- 1 / cer
      if(all(abs( cer + b * uppnt + a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      p2 <- b * cer.inv + 1
      # always disregard the left half of the parabola, go linear on the left side
      uppnt <- ifelse(uppnt >= 0, sqrt(pmax(p2^2 - ( a + b * uppnt ) * cer.inv * 2 - 2, 0)) - p2 + uppnt, -(a + (1 - uppnt) * cer) / (b + cer) )
    }
    if(!is.na(i)) warning("upper did not converge")
    bounds$upper[nontriv] <- exp(uppnt) / len[nontriv]
    nontriv <- nontriv & nonzeros
    Snt <- S[nontriv] 
    a <- - Snt * ( 1 - log(Snt) )
    if(penalty == "var") a <- a - log(sum(y))
    if(penalty=="len" || penalty=="none") a <- a -q  + pen[nontriv]
    if(penalty=="sqrt") a<- a -(q+pen[nontriv])^ 2 / 2
    lownt <- a / Snt
    b <- -Snt + ifelse(penalty == "var", 1, 0)
    #   c <- 1
    for(i in 1:max.iter) {
      cer <- exp(lownt)
      cer.inv <- 1 / cer
      if(all(abs( cer + b * lownt + a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      p2 <- b * cer.inv + 1
      # always disregard the left half of the parabola, go linear on the left side
      lownt <- ifelse(lownt >= 0, -sqrt(pmax(p2^2 - ( a + b * lownt ) * cer.inv * 2 - 2, 0)) - p2 + lownt, -(a + (1 - lownt) * cer) / (b + cer) )
    }
    if(!is.na(i)) warning("lower did not converge")
    bounds$lower[nontriv] <- exp(lownt) / len[nontriv]
  } else if(algo == "binomial") {
    # compute cumulative sums
    S <- round(S)
    sizelen <- param * len
    NS <- sizelen - S
    zeros <- if(penalty == "var") S <= 1 else S == 0
    sizes <- if(penalty == "var") NS == 1 else NS == 0
    nonzs <- !(zeros | sizes)
    sizelens <- sizelen[nonzs]
    Ss <- S[nonzs]
    NSs <- NS[nonzs]
    # compute lower bounds
    if(penalty == "var") {
      probs <- S / sizelen
      # estimate total variance
      totvar <- ( sum(y[-length(y)] * (param - y[-1])) + sum(y[-1] * (param - y[-length(y)])) ) / param / 2
#       if(totvar == 0) totvar <- sum(x) * mean(param - x) / param # estimate by constant
      A <- -q + ifelse(S == 0, 0, S * log(probs)) + ifelse(NS == 0, 0, NS * log(1 - probs)) + log(sizelen) - log(totvar)
      a <- A[nonzs]
#       print(totvar)
      Ss <- Ss - 1
      NSs <- NSs - 1
      sizelens <- sizelens - 2
      probs <- (Ss + 2) / (sizelens + 4)
#       print(rbind(Ss, NSs, sizelens, probs, a))
    } else {
      probs <- Ss / sizelens
      if(penalty=="sqrt"){
        a <- Ss * log(probs) + NSs * log(1 - probs) - 1/2 * (q + pen[nonzs])^2
      } else
        a <- Ss * log(probs) + NSs * log(1 - probs) - q + pen[nonzs]
#       print(a)
    }
    ps <- probs * 0.9
    yp <- atanh( 2 * ps - 1 )
    yprob <- atanh( 2 * probs - 1 )
    for(i in 1:max.iter) {
#       print(i);print(ps);print(Ss * log(ps) + NSs * log(1 - ps) - a)
      if(all(abs( Ss * log(ps) + NSs * log(1 - ps) - a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      a2 <- -0.5 * sizelens * ( 1 - tanh(yp)^2 )
      a1 <- Ss * (1 - tanh(yp)) - NSs * (1 + tanh(yp))
      a0 <- a - Ss * log(ps) - NSs * log(1 - ps)
      p2 <- a1 / a2 / 2
      root <- pmax(p2^2 + a0 / a2, 0)
#       print(rbind(a2, a1, a0, root))
      yp <- yp + ifelse(root > 0 & sizelens != 0, -p2 - sqrt(root), a0 / a1)
      ps <- (tanh(yp) + 1) / 2
    }
    if(!is.na(i)) warning("lower did not converge")
    bounds$lower <- rep(0, length(S))
    bounds$lower[nonzs] <- ps
    if(penalty == "var") bounds$lower[sizes & !zeros] <- exp(A / (S - 1))[sizes & !zeros] 
    if(penalty == "len" || penalty== "none") bounds$lower[sizes] <- exp((-q + pen[sizes]) / param / len[sizes])
    if(penalty == "sqrt")  bounds$lower[sizes] <- exp(-(q+pen[sizes])^2/(2*param*len[sizes]))
    
     # compute upper bounds
    if(penalty == "var") {
      zeros <- S == 1
      sizes <- NS <= 1
      nonzs <- !(zeros | sizes)
      sizelens <- sizelen[nonzs]
      Ss <- S[nonzs]
      NSs <- NS[nonzs]
      probs <- S / sizelen
      a <- A[nonzs]
      Ss <- Ss - 1
      NSs <- NSs - 1
      sizelens <- sizelens - 2
      probs <- (Ss + 2) / (sizelens + 4)
#       print(rbind(Ss, NSs, sizelens, probs, a))
    }
    ps <- 1 - ( 1 - probs ) * 0.9
    yp <- atanh( 2 * ps - 1 )
    for(i in 1:max.iter) {
      if(all(abs( Ss * log(ps) + NSs * log(1 - ps) - a ) < eps)) {
#         print(i)
        i <- NA
        break
      }
      a2 <- -0.5 * sizelens * ( 1 - tanh(yp)^2 )
      a1 <- Ss * (1 - tanh(yp)) - NSs * (1 + tanh(yp))
      a0 <- a - Ss * log(ps) - NSs * log(1 - ps)
      p2 <- a1 / a2 / 2
      root <- pmax(p2^2 + a0 / a2, 0)
      yp <- yp + ifelse(root > 0 & sizelens != 0, -p2 + sqrt(root), a0 / a1)
      ps <- (tanh(yp) + 1) / 2
    }
    if(!is.na(i)) warning("upper did not converge")
    bounds$upper <- rep(1, length(S))
    bounds$upper[nonzs] <- ps
    if(penalty == "var")  bounds$upper[zeros & !sizes] <- 1 - exp(A / (NS - 1))[zeros & !sizes]
    if(penalty == "len" || penalty== "none") bounds$upper[zeros] <- 1 - exp((-q + pen[zeros]) / param / len[zeros])
    if(penalty == "sqrt") bounds$upper[zeros] <- 1 - exp((-1/2*(q+pen[zeros])^2) / param / len[zeros])
  } else {
    stop("unknown family")
    }
  
  start <- cumsum(keep)[length(lengths) * (0:(length(y)-1)) + 1] - 1 # C-style!
  # lengthen intervals if kernel is present
  if(family == "gaussKern") {
    kj <- param$jump
    bounds$li <- bounds$li - kj + 1
    bounds$ri <- bounds$ri + kl - kj
    start <- c(start[kj:(length(y) - kl + kj)], rep(NA, kl - 1))
    # remove impossible bounds
    rem <- bounds$li < 1 | bounds$ri > length(y)
    bounds <- bounds[!rem,]
    start <- start - cumsum(rem)[start + 1]
  }
  # combine bounds for subset of right indices: each block needs to fulfill all of its bounds
  if(!missing(subset)) {
    bounds$li <- sapply(bounds$li, function(i) sum(subset < i)) + 1
    bounds$ri <- sapply(bounds$ri, function(i) sum(subset < i)) + 1
    # sort
    bounds <- bounds[order(bounds$li, bounds$ri),]
    # maybe work with diff(li) instead for recomputing start?
    start <- cumsum(sapply(tapply(bounds$li, ordered(bounds$li, levels = 1:max(bounds$ri)), identity), length))
    start <- c(0, start[-length(start)]) # C-style
    start[is.na(tapply(bounds$li, ordered(bounds$li, levels = 1:max(bounds$ri)), length))] <- NA
  }
  # check feasibility
#   feas <- sapply(which(!is.na(start)), function(i) with(bounds[bounds$li == i & bounds$ri == i,], max(lower) <= min(upper)))
  if(missing(subset)) {
    feas <- TRUE
  } else {
    si <- c(start[!is.na(start)], nrow(bounds)) + 1
    feas <- sapply(1:(length(si)-1), function(i) with(bounds[si[i]:(si[i+1]-1),], {wi <- ri == li[1]; if(any(wi)) max(lower[wi]) <= min(upper[wi]) else TRUE}))
  }
  ret <- list(bounds = bounds, start = as.integer(start), feasible = all(feas))
  class(ret) <- c("bounds", class(ret))
  ret
}

# allow subsetting by right indices (candidates)
"[.bounds" <-
function(x, subset)
{
  bounds <- x$bounds
  # combine bounds for subset of right indices: each block needs to fulfill all of its bounds
  if(!is.null(subset)) {
    bounds$li <- sapply(bounds$li, function(i) sum(subset < i)) + 1
    bounds$ri <- sapply(bounds$ri, function(i) sum(subset < i)) + 1
  }
  # sort
  bounds <- bounds[order(bounds$li, bounds$ri),]
  start <- cumsum(sapply(tapply(bounds$li, ordered(bounds$li, levels = 1:max(bounds$ri)), identity), length))
  start <- c(0, start[-length(start)]) # C-style
  start[is.na(tapply(bounds$li, ordered(bounds$li, levels = 1:max(bounds$ri)), length))] <- NA
  # check feasibility
#   feas <- sapply(which(!is.na(start)), function(i) with(bounds[bounds$li == i & bounds$ri == i,], max(lower) <= min(upper)))
  si <- c(start[!is.na(start)], nrow(bounds)) + 1
  feas <- sapply(1:(length(si)-1), function(i) with(bounds[si[i]:(si[i+1]-1),], {wi <- ri == li[1]; if(any(wi)) max(lower[wi]) <= min(upper[wi]) else TRUE}))
  list(bounds = bounds, start = as.integer(start), feasible = all(feas))
}
