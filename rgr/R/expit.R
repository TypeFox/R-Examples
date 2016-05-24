expit <-
function(z)
{
     # Function to compute the inverse-logit transformation.
     #
     p <- exp(z) / (1 + exp(z))
     return(p = p)
}
