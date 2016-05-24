##
## trend.R - trend extraction functions
##
## Authors:
##   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>
##   Dipl. Stat. Karen Schettlinger <schettlinger@statistik.uni-dortmund.de>
##
## Bugs:
##

match.trend.extraction <- function(name) {
  if(!is.character(name))
    name <- as.character(name)
  trend <- switch(name,
                  MED=trendMED,
                  RM=trendRM,                  
                  LTS=trendLTS,
                  LMS=trendLMS)
  if (is.null(trend)) {
    ## No match found. Throw error message.
    stop("Value of 'trend.extraction' must be one of 'MED', 'RM', 'LTS' or 'LMS'.")
  }
  return(trend)
}

trendMED <- function(x) {
  c(median(x), 0)
}

trendRM <- function(x) {
  n <- length(x) 
  l <- 1:n
  medquotient.t <- numeric(0)
  m <- ceiling((n-1)/2)
  for (k in l)
    medquotient.t <- c(medquotient.t, median((x[k] - x[-k]) / (k - l[-k]), na.rm=TRUE))
  
  if (all(is.na(medquotient.t))) {
    beta.RM <- NA
    if (!all(is.na(x))) {
      mu.RM <- median(x, na.rm=TRUE)
    } else {
      mu.RM <- NA
    }
  } else {
    beta.RM <- median(medquotient.t, na.rm=TRUE)
    ## intercept for x = med(x) - meaning the center of the time window!
    mu.RM <- median(x - beta.RM * (l - (m+1)), na.rm=TRUE) 
  }
  return(c(mu.RM, beta.RM))
}

trendLMS <- function(y) {
  m <- (length(y) - 1)/2
  x <- -m:m
  mdl <- lqs(x, y, method="lms")
  return(coef(mdl))
}

trendLTS <- function(y) {
  m <- (length(y) - 1)/2
  x <- -m:m
  mdl <- lqs(x, y, method="lts")
  return(coef(mdl))
}
