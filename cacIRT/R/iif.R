iif <-
function (ip, x, D = 1.7) 
	{
  
        p3 = irf(ip, x, D)$f
        q3 = 1 - p3
        
        f<-matrix(NA, nrow = length(x), ncol = dim(ip)[1])
        
        for(i in 1:length(x)){ for(j in 1:dim(ip)[1]){
        		
        f[i,j] = (D*ip[j, 1])^2*q3[i,j]/p3[i,j]*((p3[i,j]-ip[j,3])^2/(1-ip[j,3]^2))  }}
       
	    r = list(x = x, f = f)
	    return(r)
	}

