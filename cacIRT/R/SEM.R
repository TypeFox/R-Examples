SEM <-
function(ip, x, D = 1.7)
	{
		ti<-tif(ip, x, D)$f	
   		sem <- sqrt(1/ti)
   		return(sem)
   	}

