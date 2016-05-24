# calculation of the entropy
entropy <- function(x)
{
-sum(x*log(x, base=length(x)))
} 
