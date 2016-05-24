STmean <-
function(df, skew=1)
{
#(skew-1/skew) * sqrt(df) * gamma((df-1)/2)/(gamma(df/2) * sqrt(pi))
as.numeric( (skew - 1/skew) * sqrt(df) * beta((df - 1)/2, 0.5)/pi )
}
