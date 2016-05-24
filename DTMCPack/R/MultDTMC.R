MultDTMC <-
function(nchains,tmat,io,n)
{
 chains<-replicate(nchains,DTMC(tmat,io,n,trace=FALSE),simplify=FALSE) 
 return(chains)
}

