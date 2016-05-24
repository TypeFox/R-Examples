######## S-function: inv.link ########

# For computing the inverse link function
# for different families of distribution
# Currently supported families are :
#   Bernoulli           
#   Poisson   
#   Gaussian          

# Last changed: 07/15/98 

inv.link <- function(x,family)
{
   if ( family == "binomial") 
      return(1/(1+exp(-x)))

   if ( family == "poisson") 
      return(exp(x))

   if ( family == "gaussian")
      return(x)
}

########## End of inv.link ##########

