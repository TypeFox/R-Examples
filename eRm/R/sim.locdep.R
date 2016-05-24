sim.locdep <- function(persons, items, it.cor = 0.25, seed = NULL, cutpoint = "randomized")
{
# simulating data according to the local dependence model by Jannarone (1986)
# it.cor represents the pairwise item correlation. If it is a single value, it is constant over all items,
# otherwise a symmetric matrix of dimension n.items x n.items
# it.cor = 1 reflects strong violation, it.cor = 0 corresponds to the Rasch model.


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

if (is.matrix(it.cor)) {
  #if (dim(it.cor)!= c(n.items, n.items)) stop("it.cor must be symmetric and of dimension number of items")
  delta <- it.cor
} else {
  delta <- matrix(it.cor, ncol = n.items, nrow = n.items)
}

Loesprob<-matrix(0,n.persons,n.items)

if (!is.null(seed)) set.seed(seed)
Random.numbers<-matrix(runif(n.items*n.persons),n.persons,n.items)
R<-matrix(-5,n.persons,n.items)

for (j in 1:n.items)
	{
	for (i in 1:n.persons)
		{
                if ((j %% 2) == 0)
                {
                  Loesprob[i,j]<-exp(faehig[i]-schwierig[j]+(R[i,j-1]-0.5)*delta[j,j-1])/(1+exp(faehig[i]-schwierig[j]+(R[i,j-1]-0.5)*delta[j,j-1]))
		} else {
                  Loesprob[i,j]<-exp(faehig[i]-schwierig[j])/(1+exp(faehig[i]-schwierig[j]))
		}}
	R[,j]<-(Random.numbers[,j]<Loesprob[,j])*1
	}

return(R)
}

