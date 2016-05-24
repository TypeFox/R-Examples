STkurtosis <-
function(df, skew=1)
{
mueps <- STmean(df,skew=skew)
vareps <- STvar(df,skew=skew)
sdeps <- sqrt(vareps)
Eabseps2 <- df * gamma(3/2) * beta((df-2)/2,1)/sqrt(pi)
Eeps2 <- (skew^3+1/skew^3)*Eabseps2/(skew+1/skew)
Eabseps3 <- df^(1.5) * beta((df-3)/2,1.5)*2/pi
Eeps3 <- (skew^4-1/skew^4)*Eabseps3/(skew+1/skew)
#Eabseps4 <- df^2*gamma(5/2)*gamma((df-4)/2)/(gamma(1/2)*gamma(df/2))
Eabseps4 <- df^2 * 3/4 * beta((df-4)/2,2)
Eeps4 <- (skew^5+1/skew^5)*Eabseps4/(skew+1/skew)
kurtosis <- (Eeps4-4*mueps*Eeps3+6*mueps^2-3*mueps^4)/vareps^2
return(as.numeric(kurtosis))
}
