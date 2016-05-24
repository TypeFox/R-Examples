# AUXILIARY FUNCTIONS

# ----------------------------
# pooled covariance matrix
# in: D2.disc(), boxM()
pooledCov <- function(data, classes)
{
    n <- nrow(data)
    p <- ncol(data)
    stopifnot(n == length(classes))
    classes <- as.factor(classes)
    lev <- levels(classes)
    dfs <- tapply(classes, classes, length) - 1
    if (any(dfs < p)) 
       warning("such a few observations for many variables!")
    covs <- aux <- list()
    for (i in 1:nlevels(classes)) {
       covs[[i]] <- cov(data[classes == lev[i], ])
       aux[[i]] <- covs[[i]] * dfs[i]
    }
    names(covs) <- lev
    pooled <- Reduce("+", aux)/sum(dfs)
    return(pooled)
}

# -----------------------------
# insert significance symbols
# in: multicor.test()
indicate.signif <-
function(x)
{
   symbol <- NULL
   if (x <= 0.1 & x > 0.05) {
      symbol <- "."
      } else if (x <= 0.05 & x > 0.01) {
      symbol <- "*"
      } else if (x <= 0.01 & x > 0.001) {
      symbol <- "**"
      } else if (x <= 0.001) {
      symbol <- "***"
      } else {
      symbol <- " "
   }
   return(symbol)
}

# ---------------------------
# simulated p-value
# in: mantelTest()
# Simulated p-value
simpval <- 
function(null, obs, alternative)
{
    stopifnot(is.atomic(null))
    stopifnot(is.numeric(null))
    if(length(obs) != 1)
	stop("'obs' must be a vector of length 1!")
    stopifnot(is.numeric(obs))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if(alternative == "two.sided") {
          count <- 2 * min(sum(null <= -abs(obs), na.rm = TRUE),
                sum(null >= abs(obs), na.rm = TRUE))
    } else if(alternative == "less") {
          count <- sum(null <= obs, na.rm = TRUE)
    } else if(alternative == "greater") {
          count <- sum(null >= obs, na.rm=TRUE)
    }
    p <- (count + 1) / (sum(!is.na(null)) + 1)
    return(p)
}

# -------------------------------------
# in: tocher()
# aux function to find the two farthest objects
maxmat <- function(mat)
{
   n <- ncol(mat)
   v1 <- v2 <- NULL
   aux <- data.frame(v1 = rep(colnames(mat), each = n),
      v2 = rep(colnames(mat), times = n),
      val = as.vector(mat))
   aux2 <- subset(aux, v1 != v2)
   ind <- which.max(aux2[, "val"])
   ma <- aux2[ind, c("v1", "v2")]
   return(c(as.matrix(ma)))
}

# ------------------------------------
# in: gencovtest()
# collinearity analysis
conditionNumber <- function(m)
{
    eigval <- svd(m)$d
    cn <- max(eigval) / min(eigval)
    meaning <- NULL
    if (cn < 100) {
       meaning <- "weak collinearity"
    } else if (cn > 1000) {
       meaning <- "severe collinearity"
    } else {
       meaning <- "moderate to severe collinearity"
    }
    attr(cn, "meaning") <- meaning
    return(cn)
}

# --------------------------
# in: gencovtest()
# winsorized data
windata <- 
function(x, p)
{
    if(length(p) != 1 || p < 0 || p > 0.5)
       stop('"p" deve ser um valor entre 0 e 0.5!')
    qx <- quantile(x, c(p, 1-p))
    x[x < qx[1]] <- qx[1]
    x[x > qx[2]] <- qx[2]
    return(x)
}


# --------------------------------------
# on loading biotools
.welcome <- function(text = NULL)
   {
   if(is.null(text))
      text <- "Welcome to biotools!"
   if(!inherits(text, "character") || length(text) != 1)
      stop("'text' must be a character vector of length 1!")
   vec <- strsplit(text, "")[[1]]
   lab <- c(vec, "\n")
   for(i in 1:length(lab)) {
      setTxtProgressBar(txtProgressBar(char = lab[i]), 0.01)
      Sys.sleep(0.05)
   }
}

.onAttach <- function(lib, pkg)
{
   vers <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
   packageStartupMessage(.welcome(paste("---\nbiotools version", vers)))
}