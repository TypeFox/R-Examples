##load the data
data(mesa.model)

##Compute dimensions for the data structure
dim <- loglikeSTdim(mesa.model)

##Find out in which order parameters should be given
loglikeST(NULL, mesa.model)

##Let's create random vectors of values
x <- runif( dim$nparam.cov )
x.all <- runif( dim$nparam )

##Evaluate the log-likelihood for these values
loglikeST(x.all, mesa.model, "f")
loglikeST(x, mesa.model, "p")

