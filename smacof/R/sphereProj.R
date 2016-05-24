#performs sphere projection for smacofSphere.primal() and smacofSphere.dual()

sphereProj <- function(y, v, init = FALSE, immediate = FALSE, ktmax = 100,
                       itmax = 5, oeps = 1e-10, ieps = 1e-6, iverbose = FALSE,
                       overbose = FALSE) 
{
  n <- nrow(y)
  nn <- 1:n
  itel <- 1
	
  u <- v%*%y                                       #matrix U = VY
  syy <- sum(y*u)
  tau <- rowSums(u^2)                              #tau_i = u_i'u_i
	if (is.matrix(init)) z<-init else z<-y/sqrt(rowSums(y^2))
	
  repeat {
		szy <- sum(z*u)
    vz <- v%*%z
    szz <- sum(z*vz)
    rho <- (szy^2)/szz                             #rho(Z) for each single column
    
		if (overbose) cat("Iteration: ",formatC(itel,digits=3,width=3),
							"Rho: ",formatC(rho,digits=6,width=10,format="f"),
							"\n")
		kmax<-1
		
    #--------- establish and solve polynomial --> (z-vector)
    repeat {		
		  for (i in nn) {
			    ui <- u[i,]                          #init values
          zi <- z[i,]
          vz <- v%*%z
          vi <- vz[i,]
          tti <- tau[i]
          szy <- sum(z*u) 
					
          if (immediate) {
						szz <- sum(z*vz)
						rho <- (szy^2)/szz
						}
						
					hi <-ui*(szy-sum(ui*zi))-rho*(vi-v[i,i]*zi) #compute h_i 
					sold <-(sum(ui*zi)^2)+2*sum(hi*zi)          #eta block relaxation
					ppi <-(sum(ui*hi)^2)/tti                    #projector p_i
          qqi <- sum(hi^2)-ppi                        #orthogonal complement
					
          p0 <- qqi*tti^2                             #define polynomial
          p1 <- -2*tti*qqi
          p2 <- ppi+qqi-tti^2
          p3 <- 2*tti
          p4<--1
					pol<-polynomial(c(p0,p1,p2,p3,p4))          #quartic equation
          the<-max(Re(solve(pol)))                    #solve polynomial --> theta
					zi<-(hi/the)-((1/the)+(1/(tti-the)))*(sum(ui*hi)/tti)*ui
          z[i,]<-zi                                   #compute z-value [i]
					
          snew <- (sum(ui*zi)^2)+2*sum(hi*zi)         #new eta block relaxation
					if (iverbose) cat("Iteration: ",formatC(itel,digits=3,width=3),
										"sol: ",formatC(sold,digits=6,width=10,format="f"),
										"sne: ",formatC(snew,digits=6,width=10,format="f"),
										"\n")
						
			}			#end for
			if ((kmax == ktmax) || (ieps<0)) break()
			kmax<-kmax+1
		} #------------ end polynomial loop 
		
    if ((itel == itmax) || (oeps<0)) return((szy/szz)*z)    #return rho(Z,Z)
		itel <- itel+1
	}
}
