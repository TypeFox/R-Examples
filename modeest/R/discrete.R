# Author : P. Poncet

discrete <-
function(x,   # sample (the data, assumed to come from a discrete law)
         ...) # 
{
##############################################
# Estimation of the mode(s), discrete case
# The most frequent value(s) is (are) returned
##############################################

  ## Frequency table ('tabulate' is similar to 'table', but simpler)
  f <- factor(x)
  tf <- tabulate(f)

  ## Output
  return(as.numeric(levels(f)[tf==max(tf)]))
}

# Alias
mfv <- discrete # 'mfv' = 'most frequent value'


