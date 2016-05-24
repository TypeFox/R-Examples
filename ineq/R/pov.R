pov <- function(x, k, parameter = NULL, type = c("Watts", "Sen", "SST", "Foster"), na.rm = TRUE)
{
  switch(match.arg(type),
  Watts = Watts(x, k, na.rm = na.rm),
  Sen = Sen(x, k, na.rm = na.rm),
  SST = SST(x, k, na.rm = na.rm),
  Foster = Foster(x, k, parameter = parameter, na.rm = na.rm))
}

Sen <- function(x, k, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
  {
    H <- sum(x<k)/length(x)
    I <- sum((k-x2)/k)/length(x2)
    G <- Gini(x2)
    H*(I+(1-I)*G)
  }
}

SST <- function(x, k, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  x2 <- sort(x[x < k])
  n <- length(x)
  q <- length(x2)
  if(q < 1) 0 else {
    sum((2 * n - 2 * 1:q + 1) * (k - x2)/k)/n^2
  }
}

Watts <- function(x, k, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(log(k/x2))/length(x)
}

Foster <- function(x, k, parameter = 1, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  if(is.null(parameter)) parameter <- 1
  x2 <- x[x<k]
  if(length(x2)<1)
    0
  else
    sum(((k-x2)/k)^(parameter-1))/length(x)
}
