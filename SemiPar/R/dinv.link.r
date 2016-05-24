######## S-function: dinv.link ########

# For computing the derivative of the inverse 
# link function for different families of distribution
# Currently supported families are :
#   Bernoulli           
#   Poisson             

# Last changed: 01/06/98 

dinv.link <- function(x,family)
{
   if ( family == "binomial") 
   { 
     ex <- exp(x)
     return(ex/((1+ex)^2))
   }

   if ( family == "poisson") 
      return(exp(x))
}

########## End of dinv.link ##########

