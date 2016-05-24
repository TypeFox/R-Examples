

# ---------------------------------------------------------------------------------- #

moments <- function(x) {

  camp <- sort(x)
  n <- length(camp)

  m <- mean(camp)
  s <- sd(camp)
  v <- var(camp)
  cv <- s/m
  sk <- sum((camp - m)^3)/(n * s^3)
  ku <- sum((camp - m)^4)/(n * v^2) - 3

  mom <- c(m,s,cv,sk,ku)
  names(mom) <- c("mean","sd","cv","skew","kur")

  return(mom)
}


# -------------------------------------------------------------------------------- #

CV <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  m <- mean(x)
  s <- sd(x)
  cv <- s/m

  return(cv)
}


# -------------------------------------------------------------------------------- #

skew <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  m <- mean(x)
  s <- sd(x)
  #v <- var(x)
  #cv <- s/m
  sk <- sum((x - m)^3)/(n * s^3)
  #ku <- sum((x - m)^4)/(n * v^2) - 3

  return(sk)
}


# -------------------------------------------------------------------------------- #

kurt <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  m <- mean(x)
  #s <- sd(x)
  v <- var(x)
  #cv <- s/m
  #sk <- sum((x - m)^3)/(n * s^3)
  ku <- sum((x - m)^4)/(n * v^2) - 3

  return(ku)
}


