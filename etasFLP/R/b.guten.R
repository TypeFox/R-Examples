
b.guten <-
function(magn,m0=min(magn)){
            x     =subset(magn,subset=(magn>m0))
	    x     =x-m0
            n	  =length(x)
            beta  =n/sum(x)
            samplevar =beta/sum(x)
	    b	  =beta/log(10)
	    se	  =sqrt(samplevar/log(10))
	    
            return(c(b,se))
            }

