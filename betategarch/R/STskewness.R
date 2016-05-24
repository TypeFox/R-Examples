STskewness <-
function(df, skew=1)
{
mueps <- STmean(df,skew=skew)
vareps <- STvar(df,skew=skew)
sdeps <- sqrt(vareps)
#Eabseps3 <- df^(1.5)*gamma(2)*gamma((df-3)/2)/(gamma(1/2)*gamma(df/2))
Eabseps3 <- df^(1.5) * beta((df-3)/2,1.5)*2/pi
Eeps3 <- (skew^4-1/skew^4)*Eabseps3/(skew+1/skew)
skewness <- (Eeps3-3*mueps*sdeps^2*mueps^3)/sdeps^3
return(as.numeric(skewness))
}
