sim.rasch <-function(persons, items, seed = NULL, cutpoint = "randomized")
{
#produces rasch homogeneous data
#cutpoint... probability or "randomized"

if (length(items) == 1) {
  if (!is.null(seed)) set.seed(seed)
  schwierig <- rnorm(items)      #standard normal distributed
  n.items <- items
} else {
  schwierig <- items
  n.items <- length(items)
}

if (length(persons) == 1) {
  if (!is.null(seed)) set.seed(seed)
  faehig <- rnorm(persons)
  n.persons <- persons
} else {
  faehig <- persons
  n.persons <- length(persons)
}

fsmat <- outer(faehig, schwierig, "-")
psolve <- exp(fsmat)/(1+exp(fsmat))

if (cutpoint == "randomized") {
  if (!is.null(seed)) set.seed(seed)
    R <-(matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
} else {
   R <- (cutpoint < psolve)*1
 }

return(R)
}



