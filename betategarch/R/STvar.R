STvar <-
function(df, skew=1)
{
mueps <- STmean(df,skew=skew)
#Eabseps2 <- df*gamma(3/2)*gamma((df-2)/2)/(gamma(1/2)*gamma(df/2))
Eabseps2 <- df * gamma(3/2) * beta((df-2)/2,1)/sqrt(pi)
Eeps2 <- (skew^3+1/skew^3)*Eabseps2/(skew+1/skew)
vareps <- Eeps2 - mueps^2
return(as.numeric(vareps))
}
