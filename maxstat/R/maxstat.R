# $Id: maxstat.R 417 2015-10-05 18:45:38Z hothorn $


print.maxtest <- function(x, digits = getOption("digits"), ...) {
  x$stats <- NULL
  x$cuts <- NULL
  x$quant <- NULL
  x$method <- paste("Maximally selected", x$smethod,
                    "statistics using",
                    x$pmethod, collapse=" ")

  class(x) <- "htest"
  print(x, digits = digits, quote = TRUE, prefix = "", ...)
} 

print.mmaxtest <- function(x, digits = getOption("digits"), ...) {
  cat("\n\t Optimally Selected Prognostic Factors \n\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  cat("\n Selected: \n")
  p.value <- x$p.value
  sx <- x$maxstats[[x$whichmin]]
  sx$method <- x$method
  sx$stats <- NULL
  sx$cuts <- NULL
  sx$quant <- NULL
  sx$method <- paste("Maximally selected", sx$smethod, "statistics using",
                      sx$pmethod, collapse=" ")
  class(sx) <- "htest"
  print(sx, digits = digits, quote = TRUE, prefix = "", ...)
  cat("Adjusted p.value: \n")
  cat(x$p.value, ", error: ", attr(x$p.value, "error"), "\n\n")
} 

plot.maxtest <- function(x, xlab=NULL, ylab=NULL, ...) {
  xname <- unlist(strsplit(x$data.name, "by"))[2]
  if (is.na(x$quant)) {
    smethod <- x$smethod
    if (smethod == "LogRank") smethod <- "log-rank"
    if (is.null(ylab)) ylab <- paste("Standardized", smethod, "statistic")
    if (is.null(xlab)) xlab <- xname
    plot(x$cuts, x$stats, type="b", xlab=xlab, ylab=ylab, ...)
    lines(c(x$estimate, x$estimate), c(0, x$statistic), lty=2)
  } else {
    smethod <- gsub("LogRank", "log-rank", x$smethod)
    ylim <- c(min(x$quant, min(x$stats)), max(x$quant, max(x$stats)))
    ylim <- c(ylim[1]*0.95, ylim[2]*1.05)
    xlength <- range(x$cuts)
    if (is.null(ylab)) ylab <- paste("Standardized", smethod,
                                     "statistic using", x$pmethod)
    if (is.null(xlab)) xlab <- xname
    plot(x$cuts, x$stats, type="b", xlab=xlab, ylab=ylab, ylim=ylim, ...)
    lines(c(x$estimate, x$estimate), c(0, x$statistic), lty=2)
    lines(xlength, c(x$quant, x$quant), col="red")
  }
}

plot.mmaxtest <- function(x, xlab=NULL, ylab=NULL, nrow=2, ...) {
  old.par <- par(no.readonly=TRUE)
  np <- length(x$maxstats)
  ncol <- ceiling(np/nrow)
  par(mfrow=c(nrow, ncol))
  ylim <- c()
  for (i in 1:np) ylim <- c(ylim, x$maxstats[[i]]$stats) 
  ylim <- range(ylim)
  for (i in 1:np)
    plot(x$maxstats[[i]], xlab=xlab, ylab=ylab, ylim=ylim, ...)
  par(old.par)
}

pLausen92 <- function(b, minprop=0.1, maxprop=0.9)
{
  ### correction suggested by "Schell, Michael J." <Michael.Schell@moffitt.org>
  if (b < 1) return(1)
  db <- dnorm(b)
  p <- 4*db/b + db*(b - 1/b)*log((maxprop*(1 - minprop))/((1-maxprop)*minprop))
  max(p,0)
}

qLausen92 <- function(p, minprop=0.1, maxprop=0.9)
{
  test <- function(x)
    abs(pLausen92(x, minprop, maxprop) - p)

  return(optimize(test, interval=c(0,10))$minimum)
}

pLausen94 <- function(b, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  if(is.null(m))
    m <- floor(N*minprop):floor(N*maxprop)
  if (length(m) <= 1L) {
      D <- 0
  } else {
      m1 <- m[1:(length(m)-1)]
      m2 <- m[2:length(m)]
      t <- sqrt(1 - m1*(N-m2)/((N-m1)*m2))
      D <- sum(1/pi*exp(-b^2/2)*(t - (b^2/4 -1)*(t^3)/6))
  }
  1 - (pnorm(b) - pnorm(-b)) + D
}

qLausen94 <- function(p, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  test <- function(x)
    abs(pLausen94(x, N, minprop, maxprop, m) - p)

  return(optimize(test, interval=c(0,10))$minimum)
}

index <- function(x) {
  ux <- unique(x)
  ux <- ux[ux < max(ux)]
  lapply(ux, wh, x = x)
}

wh <- function(cut, x)
  which(x <= cut)

cmatrix <- function(X) {
  N <- nrow(X)
  lindx <- unlist(test <- apply(X, 2, index), recursive=FALSE)
  a <- .Call("corr", as.list(lindx), as.integer(N), PACKAGE="maxstat")
  a
}
  
pexactgauss <- function(b, x, minprop=0.1, maxprop=0.9, ...)
{
  if (length(x) > 1) {
    cm <- corrmsrs(x, minprop, maxprop) 
    if (is.null(dim(cm))) return(pnorm(-b)*2)
    p <- pmvnorm(mean=rep(0, nrow(cm)),
                 corr=cm, lower=rep(-b, nrow(cm)),
               upper=rep(b, nrow(cm)), ...)
    msg <- attr(p, "msg")
    if (msg != "Normal Completion") warning(msg)
    return(1 - p)
  }
  if (length(x) == 1) {
    return(pnorm(-b)*2)
  }
}



qexactgauss <- function(p, x, minprop=0.1, maxprop=0.9,...)
{
  test <- function(a)
    abs(pexactgauss(a, x, minprop=minprop, maxprop=maxprop, ...) - p)

  return(optimize(test, interval=c(0,10))$minimum)
}


pmaxstat <- function(b, scores,  msample, quant=FALSE)
{

  # for integers only

  if (!all(round(scores) == scores))
    stop("scores are not integers in pmaxstat")
  if (length(scores) < length(msample))
    stop("incorrect number of cutpoints in pmaxstat")
  if (length(b) != 1)
    stop("b is not a single number in pmaxstat")

  N <- length(scores)

  scores <- scores - min(scores)

  # sample sizes

  m <- 1:(N-1)

  # Expectation and Variance of a Linear Rank Test

  E <- m/N*sum(scores)
  V <- m*(N-m)/(N^2*(N-1))*(N*sum(scores^2) - sum(scores)^2)
  if (length(msample) == 1) {
    b <- b * sqrt(V[msample]) + E[msample]
    return(pperm(b, scores, m=msample, alternative="two.sided"))
  }

  #if(sum(scores) > sum(1:200)) { 
  #  warning("Cannot compute SR p-value. Sum of scores > 20100")
  #  p <- list(1, 1)
  #  names(p) <- c("upper", "lower")
  #  return(p)
  #}

  H <- rep(0, sum(scores)*N)

  totsum <- sum(scores)
  sc <- rep(1, N)

  # Streitberg / Roehmel in C, package "exactRankTest"

  H <- exactRankTests:::cpermdist2(as.integer(N),
                as.integer(totsum), as.integer(sc),
                as.integer(scores), as.logical(FALSE))

  # add last row, column for compatibility

  H <- matrix(H, nrow=N+1, byrow=TRUE)

  S <- rep(1:(ncol(H)-1), nrow(H) -2)
  S <- matrix(S, nrow(H) -2, ncol(H)-1, byrow=TRUE)
  EM <- matrix(rep(E, ncol(H) -1), nrow(H) -2, ncol(H) - 1)
  VM <- matrix(rep(V, ncol(H) -1), nrow(H) -2, ncol(H) - 1 )
  S <- (S- E)/sqrt(V)

  # remove technical parts

  H <- H[2:(nrow(H)-1), ]
  H <- H[, 2:(ncol(H))]

  # S is the matrix of the standardized values

  S <- abs(S)
  S[H == 0] <- 0

  # extend to number of permutations

  H <- H*gamma(m+1)*gamma(N -m +1)

  if (quant)
    return(list(scores=scores, H=H, E=E, S=S, msample=msample, N=N))

  # those are in general not needed

  H[S <= b] <- 0


  # delete those, which are in m+1 and + max(scores) still > b
  # well, that's the trick

  sH <- apply(H, 1, sum)

  for (i in min(msample):(nrow(H)-1)) {
    indx <- which(H[i,] > 0)
    if (length(indx) > 0) {
      indxmax <- indx[indx < E[i]]
      indxmax <- indxmax[S[i+1, indxmax + max(scores)] > b]
      if (length(indxmax) > 0 & all(!is.na(indxmax)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmax]) 
      indxmin <- indx[indx > E[i]]
      indxmin <- indxmin[S[i+1, indxmin + min(scores)] > b]
      if (length(indxmin) > 0 & all(!is.na(indxmin)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmin])
    }
  }

  # only meaningful sample sizes

  sH <- sH[msample]

  gaN <- gamma(N+1)   
  lower <- min(sum(sH)/gaN, 1)  # <- this is a better approx.
  #  upper <- max(apply(H, 1, sum))/gaN
  #  p <- list(upper, lower)
  #  names(p) <- c("upper", "lower")
  #  cat("hl working: ", all.equal(hl(scores, H, E, S, msample, N, b), p), "\n")
  lower
}

qmaxstat <- function(p, scores, msample)
{
  if (p > 1 | p < 0) stop("p not in [0,1]")
  sr <- pmaxstat(3, scores, msample, quant=TRUE)
  tr <- rev(sort(unique(round(sr$S,5))))
  i <- 1
  pw <- 0
  while (pw < p) {
    pw <- hl(sr$scores, sr$H, sr$E, sr$S, sr$msample, sr$N, tr[i])$lower
    i <- i+1
  }
  return(tr[i-1])
}


hl <- function(scores, H, E, S, msample, N, b)
{
  H[S <= b] <- 0

  # delete those, which are in m+1 and + max(scores) still > b
  # well, that's the trick

  sH <- apply(H, 1, sum)

  for (i in min(msample):(nrow(H)-1)) {
    indx <- which(H[i,] > 0)
    if (length(indx) > 0) {
      indxmax <- indx[indx < E[i]]
      indxmax <- indxmax[S[i+1, indxmax + max(scores)] > b]
      if (length(indxmax) > 0 & all(!is.na(indxmax)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmax]) 
      indxmin <- indx[indx > E[i]]
      indxmin <- indxmin[S[i+1, indxmin + min(scores)] > b]
      if (length(indxmin) > 0 & all(!is.na(indxmin)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmin])
    }
  }

  # only meaningful sample sizes

  sH <- sH[msample]

  gaN <- gamma(N+1)   
  lower <- min(sum(sH)/gaN, 1)  # <- this is a better approx.
  upper <- max(apply(H, 1, sum))/gaN
  p <- list(upper, lower)
  names(p) <- c("upper", "lower")
  p
}

### just internal functions, not exported (yet), so less sanity checking...

pmaxperm <- function(b, scores, msample, expect,
                     variance, B = 10000, ...) {
  N <- length(scores)
  if (any(msample > N)) stop("invalid split points in msample")
  p <- .Call("maxstatpermdist", scores = as.double(scores),
                           msample = as.integer(msample),
                           expect = as.double(expect),
                           variance = as.double(variance),
                           Nsim = as.integer(B),    
                           pvalonly = as.logical(TRUE),
                           ostat = as.double(b), PACKAGE = "maxstat")
  p
}

qmaxperm <- function(p, scores, msample, expect,
                     variance, B = 10000, ...) {
  N <- length(scores)
  if (length(p) > 1 || p > 1 || p < 0) stop("p must be in [0,1]")
  if (any(msample > N)) stop("invalid split points in msample")
  cp <- .Call("maxstatpermdist", scores = as.double(scores),
                           msample = as.integer(msample),
                           expect = as.double(expect),
                           variance = as.double(variance),
                           Nsim = as.integer(B),    
                           pvalonly = as.logical(FALSE),
                           ostat = NULL, PACKAGE = "maxstat")
  names(cp) <- c("T", "Prob")
  # class(cp) <- c("data.frame", "excondens")
  cs <- cumsum(cp$Prob) 
  quant <- max(cp$T[which(cs <= p)])
  RET <- list(quant = quant, exdens = cp)
  return(RET)
}
