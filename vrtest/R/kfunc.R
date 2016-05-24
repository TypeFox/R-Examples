kfunc <-
function(X)
{
WS <- (25/(12*pi^2*X^2))*( sin(6*pi*X/5)/(6*pi*X/5) - cos(6*pi*X/5) )
return(WS)
}
