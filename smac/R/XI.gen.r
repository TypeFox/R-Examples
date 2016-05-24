XI.gen = function(k)
{
XI = matrix(0,k-1,k)

XI[,1]=rep((k-1)^(-1/2),k-1)
for (ii in 2:k)
{
XI[,ii]=rep( -(1+sqrt(k))/((k-1)^(1.5)), k-1)
XI[ii-1,ii]=XI[ii-1,ii]+sqrt(k/(k-1))
}

return(XI)
}
