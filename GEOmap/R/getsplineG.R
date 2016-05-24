`getsplineG` <-
function(x,y, kdiv)
{

nlen = length(x)

n = (nlen - 1) * kdiv + 1;

ex = rep(0, n)
ey = rep(0, n)

spl = .C("CALL_jspline", PACKAGE = "GEOmap", as.double(x), as.double(y), as.integer(nlen), as.integer(kdiv),
as.double(ex), as.double(ey))



nex = as.numeric(spl[[5]])
ney = as.numeric(spl[[6]])

return(list(x=nex, y=ney))


}

