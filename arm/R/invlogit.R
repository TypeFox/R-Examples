#R function for the logistic function
logit <- function (x) {
  log(x/(1-x))
}

invlogit <- function (x) {
    1/(1+exp(-x))
}
