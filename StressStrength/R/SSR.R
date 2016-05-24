SSR <-
function(parx, pary, family="normal")
{
pnorm((parx[1]-pary[1])/sqrt(parx[2]^2+pary[2]^2))
}

