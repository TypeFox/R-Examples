dST <-
function(y, df=10, sd=1, skew=1, log=FALSE)
{
if(log){
  lnumerator <- log(2) - log(skew + 1/skew)
  ldenom1 <- lbeta(0.5, df/2) + log(sd) + log(df)/2
  ldenom2 <- (df+1) * log( 1+y^2/(skew^(2*sign(y))*df*sd^2) )/2
  result <-  lnumerator - ldenom1 - ldenom2
}else{
  numerator <- (2/(skew + 1/skew))
  denom1 <- beta(0.5, df/2) * sd * sqrt(df)
  denom2 <- (1 + y^2/(skew^(2*sign(y)) * df * sd^2 ) )^( (df+1)/2 )
  result <-  numerator/(denom1 * denom2)
}
return(result)
}
