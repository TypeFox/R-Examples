setting.interval <-
function(y, delta, s, J){
		
	# Define the indicator matrices for risk sets and failure sets 
	
	n <- length(y)
	
	smax	<- max(s)
	
	case0 <- which(delta == 0)
	case1 <- which(delta == 1)	
	
	case0yleq <- which(delta == 0 & y <= smax)
	case0ygeq <- which(delta == 0 & y > smax)
	case1yleq <- which(delta == 1 & y <= smax)
	case1ygeq <- which(delta == 1 & y > smax)
		

	ind.d <- ind.r <- matrix(0, n, J)

	for(i in case1yleq){
		d.mat.ind	<- min(which(s - y[i] >=0))
		ind.d[i, d.mat.ind]	<- 1
		ind.r[i, 1:d.mat.ind] <- 1		
		}
		
	for(i in case0yleq){
		cen.j <- min(which(s - y[i] >=0))		
		ind.r[i, 1:cen.j]	<- 1			
		}	
		
	if(length(union(case1ygeq, case0ygeq)) > 0){
		ind.r[union(case1ygeq, case0ygeq),]	<- 1
		}		
		
#	ind.r[,1]	<- 1		
		
	ind.r_d	<- ind.r - ind.d;

	d	<- colSums(ind.d)
	
	list(ind.r = ind.r, ind.d = ind.d, d = d, ind.r_d = ind.r_d)
	}

