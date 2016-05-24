fun.oldp2 <-
function(b,m,Data)
########################################################
#update for fun.ntitr(b,Z,Delta,n,m)
#ntitr.m
#######################################################
# version 0.2
# May 6, 2012
# Junlong Sun
# return [s,ru,u,gama,p,pl,deni,sm]
#######################################################
# May 6, 2012 - v0.1 Create
# May 7, 2012 - v0.2 Change multi-output interface
#######################################################
{
#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
	Z <- Data$Z
	Delta <- Data$Delta
	n <- Data$length

#-----------------------------------------------------------------#
## oldp.m
#-----------------------------------------------------------------#

	 sm <- matrix(n:1, nrow = n, ncol = 1)
       dlam <-(Delta /sm) %*% matrix(1,nrow = 1, ncol = m)
	 ww <- matrix(1,nrow = n, ncol = m)

       gama1 <- exp(-Z %*% b[1,])
       gama2 <- exp(-Z %*% b[2,])
       
       w2 <- ww*(dlam[1:n,])*gama2
       w1 <- ww*(dlam[1:n,])*gama1
		 
       lam2 <- matrix(fun.cumsum(w2) ,nrow = n, ncol = m)

       p <- exp(-lam2) 
       pl <- p
       pl[1,] <- matrix(1,nrow = 1, ncol = m)
       pl[2:n,] <- p[1:(n-1),]

	ru <- fun.cumsum(pl*w1)/p
       
#-----------------------------------------------------------------#
## ntitr.m
#-----------------------------------------------------------------#
  deni <- gama1+gama2*ru
	denil <- deni
	rul <- ru

	ow <- matrix(1,nrow = n, ncol = 1)
      u1 <- t(ow*Z*Delta)%*%(-gama1/denil)+t(ow*Z)%*%(ru/deni)
      u2 <- -t(ow*Z*Delta)%*%(gama2*rul/denil)-t(ow*Z)%*%(ru/deni)+t(ow*Z)%*%(log(deni/gama1)/gama2)
	
	u <- rbind(u1,u2)

      s <- u1^2 + u2^2

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    gama <- list(Num1=gama1,Num2=gama2)

    output <- list(s=s,ru=ru,u=u,gama=gama,p=p,pl=pl,deni=deni,sm=sm)
    return(output)
#-----------------------------------------------------------------#

}
