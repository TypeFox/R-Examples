
# Quantile function
qcgb2 <- function(prob, shape1, scale,shape2,shape3,pl0,pl,decomp="r",tol=1e-08,ff=1.5,debug=FALSE,maxiter=50)
  {
	decomp1 <- decomp
	decomp2 <- "l" 
	if (decomp == "l") {decomp2 = "r"}
	
	pr <- sort(prob)
	ord <- order(prob)
	quant <- rep(NA,length(prob))
	ltot <- length(pr)
	lp0 <- length(pr[pr==0])
	lp1 <- length(pr[pr==1])
	lp <- lp0+lp1 # new 11.07.12
        quant[pr==1] <- Inf
        quant[pr==0] <- 0
        if (lp==length(pr)) return(quant)  # new 11.07.12



	  L1 <- seq_len(length(pr[(0<pr) & (pr<=0.5)]))
	  if(length(L1) > 0)
		{
		x0 <- moment.cgb2(1,shape1,1,shape2,shape3,pl0,pl,decomp1)  #new 11.07.12; scale=1
		p0 <- prcgb2(0,x0,shape1,1,shape2,shape3,pl0,pl,decomp1)      #new 11.07.12: scale=1
		pr1 <- 0
		q1 <- 0	
		pr2 <- pr[pr>0][1]

		cc <- 0
		while ((pr2 < p0) & (cc < maxiter))
			{
			cc <- cc+1
			x0 <- x0/ff        # new 11.07.12
			p0 <-prcgb2(q1,x0,shape1,1,shape2,shape3,pl0,pl,decomp1,tol,debug=FALSE)  #new 11.07.12
#			print(c(2,cc,x0,p0))
			}
		for (kk in L1)
			{
			i <- 0
			k <- kk + lp0
			pr2 <- pr[k]
			while (i < maxiter)
				{
                		delx <- (pr2-p0)/dcgb2(x0,shape1,1,shape2,shape3,pl0,pl,decomp1)
                		x0 <- max(tol,x0 + delx)
                		p0 <- pr1 + prcgb2(q1,x0,shape1,1,shape2,shape3,pl0,pl,decomp1,tol,debug=FALSE)

                		if (isgood(delx, tol)) 
                        		{ 
                        			if (debug)  print(as.vector(c("iterations=",i,"obs=",k)),quote=F)
                        			i <- maxiter+1
                               		}
                               else i <- i+1
				}


                	if (i==maxiter) warning("series not converged")
                        
			quant[k] <- scale*x0                     #new 11.07.12
			pr1 <- pr[k] 
			q1 <- x0
			}
		}
		
	  L2 <- seq_len(length(pr[(0.5<pr) & (pr<1)]))
	  if (length(L2) >0)
		{
		x0 <- moment.cgb2(1,shape1,1,shape3,shape2,pl0,pl,decomp2)  #new 11.07.12; scale=1
		p0 <- prcgb2(0,x0,shape1,1,shape3,shape2,pl0,pl,decomp2)      #new 11.07.12: scale=1
		pr1 <- 0
		q1 <- 0	
		pr2 <- 1-pr[ltot-lp1]

		cc <-0
		while ((pr2 < p0) & (cc < maxiter)) 
			{
			cc <- cc+1
			x0 <- x0/ff         # new 11.07.12
			p0 <-prcgb2(q1,x0,shape1,1,shape3,shape2,pl0,pl,decomp2,tol,debug=FALSE)  #new 11.07.12
#			print(c(4,cc,x0,p0))
			}
		 for (kk in L2)
			{
			i <- 0
			k <- -kk - lp1+ltot +1
			pr2 <- 1- pr[k]
			while (i < maxiter)
				{
                		delx <- (pr2-p0)/dcgb2(x0,shape1,1,shape3,shape2,pl0,pl,decomp2)
                		x0 <- max(tol,x0 + delx)
                		p0 <- pr1 + prcgb2(q1,x0,shape1,1,shape3,shape2,pl0,pl,decomp2,tol,debug=FALSE)

                		if (isgood(delx, tol)) 
				  { 
				  if (debug)  print(as.vector(c("iterations=",i,"obs=",k)),quote=F)
				  i <- maxiter+1
				  }
				else i <- i+1
				}
                	if (i==maxiter) warning("series not converged")
                        
			quant[k] <- scale/x0                     #new 11.07.12
			pr1 <-1- pr[k] 
			q1 <- x0
			}
		}
		
  return(quant[order(prob)])
  }
    
# Random generation
rcgb2 <- function(n,shape1,scale,shape2,shape3,pl0,pl,decomp="r",tol=1e-02,maxiter=100,debug = FALSE){
ranu <- runif(n,0,1)
rand <- rep(NA,n)
rand <-qcgb2(ranu,shape1,scale,shape2,shape3,pl0,pl,decomp=decomp,tol=tol,maxiter=maxiter, debug=debug)
return(rand)
}