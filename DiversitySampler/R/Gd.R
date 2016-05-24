Gd <-
function(obs,scrub=FALSE){
				if(scrub==TRUE)obs <- obs[obs>0]
				nsp <- length(obs)
                nsamp <- sum(obs)
                obs = obs/nsamp
                Gd<-(1-sum(obs^2))*(nsp/(nsp-1))
                return(Gd)
                }

