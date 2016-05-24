initsPrint <- function(x) {
  x$B <- NULL
  var <- initsVar(x)
  paste("list(", paste(var, collapse = ",\n "), ")\n")
}

initsVar <- function(x) {
  n <- length(x)
  var <- vector(length = n)
  for (i in seq_len(n)) var[i] <- paste(names(x)[i], " = c(", paste(x[[i]], collapse = ", "),")", sep = "")
  var
}

initsFS <- function(x, k, restrict, initialValues = list()) {
  if (missing(restrict)) restrict <- ""
  x <- as.matrix(x)
  if (any(names(initialValues) %in% c("mu", "eta", "tau"))) stop("initialValues are not specified correctly")
  eta <- rep(1/k, k)
  eta[k] = 1-sum(eta[-k])
  if (restrict == "mu")  mu <- mean(x)
  else mu <- stats::quantile(x, probs = seq(1/(k+1),k/(k+1), length.out = k))
  names(mu) <- NULL
  R <- stats::IQR(x)
  sigma2 <- (R/1.34)^2
  if (restrict == "tau") tau <- 1/sigma2
  else tau <- rep(1/sigma2,k)
  z <- c(list(eta = eta, mu = mu, tau = tau), initialValues)
  z
}

