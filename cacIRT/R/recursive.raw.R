recursive.raw <-
function(ip, theta, D = 1.7){
	
	ut <- theta
       
	ni <- dim(ip)[1]
	nn <- length(ut)
	sc <- ni+1
	
	 	   
		Pjt <-  irf(ip, ut, D)$f
		Qjt <-  1-Pjt
		H <- apply(Qjt,1,prod)
		Zjt	<- Pjt/Qjt

	zm <- array(0, dim = c(nn, sc, ni))	
		zm[,1,] <- 1
		zm[,2,1] <- Zjt[,1]

	for(i in 2: ni){ # The recursive part
			for(s in 2:sc){
				zm[,s,i] <- zm[,s,i-1] + Zjt[,i]*zm[,s-1,i-1]
			}
	  }	
	
		out<- zm[,,ni]*H
				rownames(out) <- paste("theta = ",round(ut,3))
				colnames(out) <-paste("x = " , 0:ni)
	return(out)
	}

