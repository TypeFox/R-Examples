fun.newBest <-
function(Data,u0,...)
{

n <- Data$length
oz <- Data$Z
Z <- Data$Z
od <- Data$Delta

#fun.hesscvx0Function(u0,oz,od,n)
nb <- matrix(c(0,0), nrow = 1, ncol = 2)
ob <- nb+.1

bc <- t(nb)
jh <- .01

#fun.hcvx0(bc,jh,u0,oz,od,n)
#------------------------------------------------------------------------------------------------------------------------#
	m <- 5
	b <- bc %*% matrix(1, nrow = 1, ncol = m) + matrix(c(0,0,jh,0,-jh,-0,0,jh,-0,-jh),nrow = 2, ncol = m)

  ## ntitr0.m
  data11 <- fun.ntitr0(u0,m,b,oz,od,Z,n)
  u <- data11$u
  s <- data11$s
  ru <- data11$ru

	h <- t(cbind((u[,2]-u[,3])/jh/2,(u[,4]-u[,5])/jh/2))
      h[1,2] <- h[1,2]/2 + h[2,1]/2
      h[2,1] <- h[1,2]

    ep <- 0.001
    cr <- 0

    #-----------------------------------------------------------------#
    while  (((det(h) > 1.e-4) * (h[1,1] >0) * (h[2,2]>0)<1) & (cr<100)) {
        ep <- 2*ep
        h <- h + ep * matrix(diag(rep(1,2)),nrow=2,ncol=2)
    }
#-----------------------------------------------------------------#

    po <- -solve(h) %*% u[,1]
#------------------------------------------------------------------------------------------------------------------------#

ou <- u[,1]
of <- s[1]

nm <- 4
ittm <- 50
cn <- 0

#while (sum(abs(round(1000*ob)-round(1000*nb)))>0)&(cn<ittm)
while ((sum(abs(round(1000*ob) - round(1000*nb)))>0) & (cn<ittm) ) {

    cn <- cn +1
    t <- 1
    bc <- t(nb) + t*po
    idb <- sum(abs(bc)>nm)
    tb <- 0
        
   #-----------------------------------------------------------------#
   while ((idb>0)&(tb<20)){
        tb <- tb+1
        t <- t*0.5
        bc <- t(nb)+ t*po
        idb <- sum(abs(bc)>nm)
    }
   #-----------------------------------------------------------------#

    m <- 1
    b <- bc
    
   ##return [s,ru]
   #data4 <- fun.ntitr(b,Z,Delta,n,m1)
   #s <- data4$s
   #ru <- data4$ru

	## ntitr0.m
	#------------------------------------#
  data12 <- fun.ntitr0(u0,m,b,oz,od,Z,n)
  u <- data12$u
  s <- data12$s
  ru <- data12$ru
	#------------------------------------#

    if (of==0){
        break
    } else if (abs((s/of-1)*1.e+8)<1){
        break
    }
 

    ob <- nb
    nb <- t(bc)
    of <- s

   ##return [po,s]
   #data3 <- fun.hcvxitr(bc,jh,Z,Delta,n)
   #po <- data3$po
   #s <- data3$s

   #data3 <- fun.hcvxitr1(bc,Data,Parameters)
   #po <- data3$po
   #s <- data3$s
   #fun.hcvx0(bc,jh,u0,oz,od,n)
#------------------------------------------------------------------------------------------------------------------------#
	m <- 5
	b <- bc %*% matrix(1,nrow = 1, ncol = m) + matrix(c(0,0,jh,0,-jh,-0,0,jh,-0,-jh),nrow = 2, ncol = m)

	## ntitr0.m
  data13 <- fun.ntitr0(u0,m,b,oz,od,Z,n)
  u <- data13$u
  s <- data13$s
  ru <- data13$ru
	#------------------------------------#

	h <- t(cbind((u[,2]-u[,3])/jh/2,(u[,4]-u[,5])/jh/2))
      h[1,2] <- h[1,2]/2 + h[2,1]/2
      h[2,1] <- h[1,2]

    ep <- 0.001
    cr <- 0

    #-----------------------------------------------------------------#
    while  (((det(h) > 1.e-4) * (h[1,1] >0) * (h[2,2]>0)<1) & (cr<100)) {
        ep <- 2*ep
        h <- h + ep * matrix(diag(rep(1,2)),nrow=2,ncol=2)
    }
#-----------------------------------------------------------------#

    po <- -solve(h) %*% u[,1]
#------------------------------------------------------------------------------------------------------------------------#


}

b <- cbind(t(nb), t(ob))
m <- 2

##return [s,ru]
#data5 <- fun.ntitr(b,Z,Delta,n,m)
#s <- data5$s
#ru <- data5$ru

## ntitr0.m
  data14 <- fun.ntitr0(u0,m,b,oz,od,Z,n)
  u <- data14$u
  s <- data14$s
  ru <- data14$ru
	#------------------------------------#

#data5 <- fun.oldp2(b,m,Data)
#s <- data5$s
#ru <- data5$ru
        
ots <- min(s)
jt <- which.min(s)
best <- t(b[,jt])
r <- ru[,jt]

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(best=best,r=r)
    return(output)
#-----------------------------------------------------------------#


#fun.hcvx0Function(bc,jh,u0,oz,od,n)

}
