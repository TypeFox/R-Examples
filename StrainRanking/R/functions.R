
## Parameter functions used in the simulation of the regression model
generation.alpha.3strains <- function(x) {
	alpha=NULL
	alpha <- cbind(alpha,(cos(x[,2]) + (3/2))*100)
	alpha <- cbind(alpha,(sin(x[,1]) + (3/2))*100)
	alpha <- cbind(alpha,(sin(x[,2]) + (3/2))*100)
	return(alpha)
}

## Function generating draws from the Dirichlet distribution
## code taken from the R-package "gregmisc"
.rdirichlet=function (n, alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

## Function providing the intensity of the risk of infection
.infection.potential=function(par,r,y){
  colSums(y*par[1]*exp(-r/par[2]))
}

## Epanechnikov smoothing kernel
.noyau <- function(u) {
	if (u >= 0  &  u <= 1) {
    	K <- 1 - u^2
  	} else {
    	K <- 0
  	}
  	return(K)
}

## Distance function
.distance <- function(x1i,x2i,jeuj) {
	d <- ((jeuj[1] - x1i) ^ 2) + ((jeuj[2] - x2i) ^ 2)
	return(d)
}

## Function computing the weights of the kernel smoothing
.calcul.wij <- function(jeu,j,x1i,x2i,b){
  	if (jeu[j,1] == x1i & jeu[j,2] == x2i) {
    	wij <- 1
  	} else {
    	wij <- .noyau(.distance(x1i, x2i, jeu[j,]) / b)
  	}
  	return(wij)
}

## Function estimating proportions of strains
.estimation.piS <- function(jeu,x1i,x2i,s,b) {
  	sommeWij <- 0
  	for (j in 1:nrow(jeu)) {
    	sommeWij <- sommeWij + .calcul.wij(jeu,j,x1i,x2i,b)
  	}
  	sommeWijS <- 0
  	for (j in 1:nrow(jeu)) {
    	sommeWijS <- sommeWijS + .calcul.wij(jeu,j,x1i,x2i,b) * jeu[j,2+s] /	sum(jeu[j,-(1:2)])
  	}
  	return(sommeWijS / sommeWij)
}

