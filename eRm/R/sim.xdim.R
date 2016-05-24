sim.xdim <- function(persons, items, Sigma, weightmat, seed = NULL, cutpoint = "randomized")
{

# Sigma ... VC matrix for multinormal distribution
# weightmat ... matrix of dimension k times D with weights. If omitted, equal weights are used.

if (missing(Sigma)) {
  ndim <- ncol(persons)
} else {
  ndim <- nrow(Sigma)                      #number of dimensions
}

if (length(persons) == 1) {                #simulating
  if (!is.null(seed)) set.seed(seed)
  faehig <- mvrnorm(persons, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
} else {
  faehig <- persons
}
if (length(items) == 1) {
  if (!is.null(seed)) set.seed(seed)
  schwierig <- rnorm(items)
} else {
  schwierig <- items
}


n.persons <- nrow(faehig)
n.items <- length(schwierig)

if (missing(weightmat)) {                      #specifying the weight matrix
  weightmat <- matrix(0, ncol = ndim, nrow = n.items)
  if (!is.null(seed)) set.seed(seed)
  indvec <- sample(1:ndim, n.items, replace = TRUE)
  for (i in 1:n.items) weightmat[i,indvec[i]] <- 1
}

Wp <- apply(weightmat, 1, function(wi) {      #n.persons times n.items matrix
                     Xw <- t(wi) %*% t(faehig)})

psolve <- matrix(0,n.persons,n.items)

#class<-rep(1,n.persons)
#class[sample(n.persons)[1:round(n.persons/2,0)]]<-2

for (j in 1:n.items)
  for (i in 1:n.persons)
    psolve[i,j] <- exp(Wp[i,j]-schwierig[j])/(1+ exp(Wp[i,j]-schwierig[j]))

if (cutpoint == "randomized") {
  if (!is.null(seed)) set.seed(seed)
  R <-(matrix(runif(n.items*n.persons),n.persons,n.items) < psolve)*1
} else {
  R <- (cutpoint < psolve)*1
}

return(R)
}


