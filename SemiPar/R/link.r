######## S-function: link ########

# For computing the link function
# for different families of distribution
# Currently supported families are :
#   Bernoulli           
#   Poisson             
#   Gaussian

# Last changed: 07/15/98 

link <- function(x,family)
{
   if ( family == "binomial") 
      return(log(x/(1-x)))

   if ( family == "poisson") 
      return(log(x))

   if ( family == "gaussian")
      return(x)
}

########## End of link ##########

