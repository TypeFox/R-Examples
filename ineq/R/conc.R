conc <- function(x, parameter = NULL, type = c("Herfindahl", "Rosenbluth"), na.rm = TRUE)
{
  switch(match.arg(type),
  Herfindahl = Herfindahl(x, parameter = parameter, na.rm = na.rm),
  Rosenbluth = Rosenbluth(x, na.rm = na.rm))
}

Herfindahl <- function(x, parameter = 1, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  m <- if(is.null(parameter)) 1 else parameter
  Herf <- x/sum(x)
  Herf <- Herf^(m+1)
  Herf <- sum(Herf)^(1/m)
  Herf
}

Rosenbluth <- function(x, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  x <- sort(x)
  HT <- (n:1)*x
  HT <- 2*sum(HT/sum(x))
  HT <- 1/(HT-1)
  HT
}

