gambin_exp <-
function(alpha, maxoctave, total_species)
{
  exp <- dgambin(alpha, maxoctave)
  exp <- exp * total_species 
  return(exp)
}
