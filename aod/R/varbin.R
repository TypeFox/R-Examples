varbin <- function(n, y, data, alpha = 0.05, R = 5000){
  CALL <- match.call()
  resp <- c(deparse(substitute(n)), deparse(substitute(y)))
  datan <- data[ , resp]
  if(any(datan[ , 2] > datan[ , 1]))
    stop("Some ", deparse(substitute(y)), " were > ", deparse(substitute(n)), ".")
  if(any(datan[ , 1] <= 0))
    warning("Data with ", deparse(substitute(n)), " <= 0 were discarded.")
  names(datan) <- c("n", "y")
  datan <- datan[datan$n > 0, ]
  n <- datan$n; y <- datan$y; N <- length(n); ntot <- sum(n); ytot <- sum(y)
  # binomial
  p <- ytot / ntot
  varp <- p * (1 - p) / (ntot - 1)
  # ratio
  pratio <- p
  varpratio <- N * (N - 1)^(-1) * ntot^(-2) * sum((y - n * p)^2)
  # arithmetic
  parithm <- mean(y / n)
  varparithm <- var(y / n) / N
  # jackknife
  z <- (y - n * p) / (ntot - n)
  pseudovalue <- p + (N - 1) * z
  #pjack <- p + (N - 1) * mean(z)
  #varpjack <- (N^(-1)) * (N - 1) * sum((z - mean(z))^2)
  pjack <- mean(pseudovalue)
  varpjack <- var(pseudovalue) / N
  # bootstrap
  if(!require(boot, quietly = TRUE))
    stop("This function requires the recommended package dQuote(boot).")
  foo <- function(d, f) p <- sum(d$y * f) / sum(d$n * f)
  res <- boot(data = datan, statistic = foo, stype = "w", R = R)
  pboot <- mean(res$t)
  varpboot <- var(res$t)
  # results
  tab <- data.frame(p = c(p, pratio, parithm, pjack, pboot),
                    varp = c(varp, varpratio, varparithm, varpjack, varpboot))
  tab$lower <- pmax(0, tab$p - qnorm(1 - alpha/2) * sqrt(tab$varp))
  tab$upper <- pmin(1, tab$p + qnorm(1 - alpha/2) * sqrt(tab$varp))
  rownames(tab) <- c("Binomial", "Ratio", "Arithmetic", "Jackknife", "Bootstrap")
  features <- c(N = N, n = ntot, y = ytot)
  # outputs
  new(Class = "varbin", CALL = CALL, tab = tab, pboot = res$t, alpha = alpha, features = features)
  }    

# show
setMethod(f = "show", signature = "varbin", definition = function(object){                   
  summry <- object@tab
  nam <- rownames(summry)
  summry$varp <- sqrt(summry$varp)
  names(summry)[2] <- "se"
  print(object@CALL)
  feat <- object@features
  cat("N = ", feat["N"], " clusters, n = ", feat["n"], " subjects, y = ", feat["y"], " cases.\n", sep = "")
  cat("\nProportion", "\n----------\n", sep = "")
  List <- lapply(summry, function(x) ifelse(is.na(x), "", format(round(x, digits = 4), nsmall = 3)))
  summ <- as.data.frame(t(do.call("rbind", List)))
  rownames(summ) <- nam
  print(summ)
  cat("\nalpha = ", object@alpha, " for the CIs; R = ", length(object@pboot),
      " samples for the bootstrap estimates.\n", sep = "")
  invisible(summry)
  })

