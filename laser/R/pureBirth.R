`pureBirth` <-
function(x)
{
  # calculates ML estimates of speciation rate
  #under yule model & associated LH of full set of branching times
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  res <- vector("list", 3)
  x <- rev(sort(x))
  temp <- yuleint2(x, x[1], 0)
  res$LH <- temp[2]
  res$aic <- (-2* temp[2]) + 2
  res$r1 <- temp[1]
  return(res)
}



#`IpureBirth` <-
#function(x)
#{
#  res <- list()
#  temp <- yuleint2(x, x[1], 0)
#  res$LH <- temp[2]
#  res$r1 <- temp[1]
#  return(res)
#}

