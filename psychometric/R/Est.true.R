"Est.true" <-
function (obs, mx, rxx)
{
that <- mx*(1-rxx)+rxx*obs
return(that)
}

