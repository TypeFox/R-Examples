sim.2pl <- function(persons, items, discrim = 0.25, seed = NULL, cutpoint = "randomized")
{

# simulation of Birnbaum's 2-PL (non-parallel ICCs)
# violation is steered by the standard deviation sdlog.
# meanlog in rlnorm is 0 which implies that the random numbers lie asymmetrically around 1. If sdlog = 0 the data are Rasch homogeneous.
# For IRT applications, values up to 0.5 should be considered. 

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


if (length(discrim) > 1) {
 alpha <- discrim 
} else {
 if (!is.null(seed)) set.seed(seed) 
 alpha <- rlnorm(n.items, 0, sdlog = discrim)         #discrimination parameter
}

psolve <- matrix(0, n.persons, n.items)

for (i in 1:n.persons)
	for (j in 1:n.items)
	psolve[i,j]<-exp(alpha[j]*(faehig[i]-schwierig[j]))/(1+exp(alpha[j]*(faehig[i]-schwierig[j])))

if (cutpoint == "randomized") {
  if (!is.null(seed)) set.seed(seed)
    R <-(matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
} else {
    R <- (cutpoint < psolve)*1
}




return(R)
}


