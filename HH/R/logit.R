### logit and antilogit

logit <- function(p) log(p/(1-p))

##antilogit <- function(x) exp(x)/(1+exp(x))  #boundary problem when x==Inf

antilogit <- function(x) {
  tmp <- (x != Inf)
  result <- x
  result[tmp] <- exp(x[tmp])/(1+exp(x[tmp]))
  result[!tmp] <- 1
  result
}


#### odds and antiodds

odds <- function(p) p/(1-p)

antiodds <- function(o) {
  p <- o/(o+1)
  p[is.nan(p)] <- 1
  p
}
