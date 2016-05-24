getsets = function(x){
	s = x[,1]
	s2 = x[,2]
	ind = seq_along(s)
	ind2 = seq_along(s2)
	while(length(ind) > 0){
		
		s[ind] = sapply(s[ind],function(z){
			c(z[z<0], x[ z[z>0], ])
		},simplify=FALSE)
		ind = which(sapply(s, function(v){
				any(v > 0)
			}))
		}
	while(length(ind2) > 0){	
		s2[ind2] = sapply(s2[ind2],function(z){
			c(z[z<0], x[ z[z>0], ])
		},simplify=FALSE)
		ind2 = which(sapply(s2, function(v){
				any(v > 0)
			}))
		}
		return(list(s,s2))
}


