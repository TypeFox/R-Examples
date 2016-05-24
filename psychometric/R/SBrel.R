"SBrel" <-
function(Nlength, rxx)
{
rxxp <- Nlength*rxx/(1+(Nlength-1)*rxx)
return(rxxp)
}

