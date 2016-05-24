VAR.modul <-
function(b,p)
{k <- nrow(b)
if(p == 1) bmat <- b[,1:(k*p)]
if(p > 1) bmat <- rbind( b[,1:(k*p)], cbind(diag(k*(p-1)),matrix(0,nrow=k*(p-1),ncol=k)) )
modul <- Mod(eigen(bmat)$value)
return(modul)
}
