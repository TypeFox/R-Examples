Hs <-
function(obs,base=exp(1),corr=FALSE, scrub=TRUE){
				if (scrub==TRUE) obs <- obs[obs>0]
				nsp <- length(obs)
                nsamp <- sum(obs)
                obs = obs/nsamp
                lnobs = obs * 0
                for (k in 1:nsp) {
                  if (obs[k] == 0) {
                    lnobs[k] = 0
                  }
                  else {
                    lnobs[k] = log(obs[k], base)
                  }
                }
                if (corr==FALSE) {H = -sum(obs * lnobs)}
                if (corr==TRUE) {H = -sum(obs * lnobs) - ((nsp-1)/(2*nsamp))}
                return(H)
                }

