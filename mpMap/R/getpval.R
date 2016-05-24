getpval <- 
function(model)
{
  summ <- summary(model)$fstatistic
  return(1-pf(summ[1], summ[2], summ[3]))
}
