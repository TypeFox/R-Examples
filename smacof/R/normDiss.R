`normDiss` <-
function(diss,wghts)
{
  return(diss/sqrt(sum(wghts*diss^2)))
}

