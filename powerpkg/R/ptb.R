"ptb" <-
function(g,m,ld,nfam,alpha)
{
# Create a power table for different values of q
#  ptb(g=1.5,m=0.5,ld=0.75,nfam=2000,alpha=0.00000005)
 qv <- seq(0.2,0.8,0.1)
 tb <- NULL
 for (q in qv)
   {
   res <- ptdt(g,q,m,ld,nfam,alpha)
   tb <- rbind(tb,unlist(res))
   }
 return(tb)
}

