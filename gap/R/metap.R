# 7-3-2008 MRC-Epid JHZ

metap <- function(data,N,verbose="Y",prefixp="p",prefixn="n")
{
   M <- dim(data)[1]
   rawp <- rawn <- w <- matrix("numeric",M,N)
   x2 <- p <- metaz <- vector("numeric",M)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           rawp[j,i] <- data[paste(prefixp,i,sep="")][j,]
           rawn[j,i] <- data[paste(prefixn,i,sep="")][j,]
       }
       s1 <- rawn[j,]
       s2 <- as.numeric(s1)
       w[j,] <- sqrt(s2)/sqrt(sum(s2))
       t1 <- rawp[j,]
       t2 <- as.numeric(t1)
       x2[j] <- 2*sum(-log(t2))
       metaz[j] <- sum(as.numeric(w[j,])*qnorm(1-t2/2))
   }
   p <- pchisq(x2,2*N,lower.tail=FALSE)
   metap1 <- pnorm(-abs(metaz))
   metap2 <- 2*metap1
   if (toupper(verbose)=="Y")
   cat("\n\nFisher's method x2=",x2,"\ndf=",2*N,"\ntwo-sided p=",p,"\n\n")
   cat("\n\nCombined z=",metaz,"\none-sided p=",metap1,"\ntwo-sided p=",metap2,"\n\n")
   invisible(data.frame(x2=x2,p=p,z=metaz,p1=metap1,p2=metap2))
}
