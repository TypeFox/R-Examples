ygupdate <-
function(mod,q,G,yg,p,sg,lg,betag,thetag,pig,om,tri){
	if(substring(mod,1,1)=="U"){
		if(substring(mod,2,2)=="U"){
			if(substring(mod,3,3)=="C"){
				### UUC model
				if(q==1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*%
							t(betag[,,g]) %*% sg[,,g]))*diag(p)
					}
				}
				if(q>1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*% 
							betag[,,g] %*% sg[,,g]))*diag(p)
					}
				}
			}
			if(substring(mod,3,3)=="U"){
				### UUU model
				if(q==1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% 
								t(betag[,,g]) %*% sg[,,g])
					}
				}
				if(q>1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*%
								betag[,,g] %*% sg[,,g])
					}
				}
			}
		}
		if(substring(mod,2,2)=="C"){
			if(substring(mod,3,3)=="C"){
				### UCC model
			  ygc <- matrix(0,p,p)
				if(q==1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*%
							t(betag[,,g]) %*% sg[,,g]))*diag(p)
						ygc <- ygc + pig[g]*yg[,,g]
					}
					for(g in 1:G){
						yg[,,g] <- ygc
					}
				}
				if(q>1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*% 
							betag[,,g] %*% sg[,,g]))*diag(p)
						ygc <- ygc + pig[g]*yg[,,g]
					}
					for(g in 1:G){
						yg[,,g] <- ygc
					}
				}
			}
			if(substring(mod,3,3)=="U"){
				### UCU model
			  ygc <- matrix(0,p,p)
				if(q==1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% 
								t(betag[,,g]) %*% sg[,,g])
						ygc <- ygc + pig[g]*yg[,,g]
					}
					for(g in 1:G){
						yg[,,g] <- ygc
					}
				}
				if(q>1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*%
								betag[,,g] %*% sg[,,g])
						ygc <- ygc + pig[g]*yg[,,g]
					}
					for(g in 1:G){
						yg[,,g] <- ygc
					}
				}
			}
		}	
	}
	if(substring(mod,1,1)=="C"){
		if(substring(mod,2,2)=="U"){
			if(substring(mod,3,3)=="C"){
				### CUC model
				if(q==1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - 2*lg[,,g] %*%
							t(betag[,,g]) %*% sg[,,g] + (lg[,,g] * thetag[,,g]) %*% 
							t(lg[,,g])))*diag(p)
					}
				}
				if(q>1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - 2*lg[,,g] %*%
							betag[,,g] %*% sg[,,g] + lg[,,g] %*% 
							thetag[,,g] %*% t(lg[,,g])))*diag(p)
					}
				}
			}
			if(substring(mod,3,3)=="U"){
				### CUU model
				if(q==1){
					for(g in 1:G){
						yg[,,g] <- diag(sg[,,g] - 2*lg[,,g] %*%
							t(betag[,,g]) %*% sg[,,g] + (lg[,,g] * thetag[,,g]) %*% 
							t(lg[,,g])) * diag(p)
					}
				}
				if(q>1){
					for(g in 1:G){
						yg[,,g] <- diag(sg[,,g] - 2*lg[,,g] %*%
							betag[,,g] %*% sg[,,g] + lg[,,g] %*% 
							thetag[,,g] %*% t(lg[,,g])) * diag(p)
					}
				}
			}
		}
		if(substring(mod,2,2)=="C"){
			if(substring(mod,3,3)=="U"){
				### CCU model
				if(q==1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% 
								t(betag[,,g]) %*% sg[,,g])
					}
				}
				if(q>1){
					for(g in 1:G){
						diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*%
								betag[,,g] %*% sg[,,g])
					}
				}
			}
			if(substring(mod,3,3)=="C"){
				### CCC model
				if(q==1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*%
							t(betag[,,g]) %*% sg[,,g]))*diag(p)
					}
				}
				if(q>1){
					for(g in 1:G){
						yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*% 
							betag[,,g] %*% sg[,,g]))*diag(p)
					}
				}
			}
		}
	}
	yg
}
