triupdate <-
function(mod,q,G,yg,p,sg,lg,betag,thetag,pig,om,tri,ng){
	if(substring(mod,3,3)=="1"){
		E <- matrix(0,p,p)
		if(q==1){
			for(g in 1:G){
				E <- E + (ng[g]/om[g])*(sg[,,g]-2*lg[,,g]%*%t(betag[,,g])%*%sg[,,g]+
							lg[,,g]%*%t(thetag[,,g])%*%t(lg[,,g]))
			}
		}
		if(q>1){
			for(g in 1:G){
				E <- E + (ng[g]/om[g])*(sg[,,g]-2*lg[,,g]%*%betag[,,g]%*%sg[,,g]+
										lg[,,g]%*%t(thetag[,,g])%*%t(lg[,,g]))
			}
		}
		k <- (prod(diag(E))^(1/p)-sum(ng))/2
		findum <- diag(E)/(sum(ng) + 2*k)
		tri[,,] <- diag(findum)
	}
	if(substring(mod,3,3)=="2"){
		E <- matrix(0,p,p)
		if(q==1){
				for(g in 1:G){
					E <- E + (ng[g]/om[g])*(sg[,,g]-lg[,,g]%*%t(betag[,,g])%*%sg[,,g])
				}
		}
		if(q>1){
				for(g in 1:G){
					E <- E + (ng[g]/om[g])*(sg[,,g]-lg[,,g]%*%betag[,,g]%*%sg[,,g])
				}
		}
		k <- (prod(diag(E))^(1/p)-sum(ng))/2
		findum <- diag(E)/(sum(ng) + 2*k)
		tri[,,] <- diag(findum)  
	}
	if(substring(mod,3,3)=="3"){
		E <- array(0, dim=c(p,p,G))
		k <- rep(0,G)
		if(q==1){
			for(g in 1:G){
				diag(E[,,g]) <- diag(sg[,,g] - 2*lg[,,g]%*%t(betag[,,g])%*%sg[,,g] + 
							lg[,,g]%*%t(thetag[,,g])%*%t(lg[,,g]))
				k[g] <- (ng[g]/2)*((prod(diag(E[,,g]))^(1/p))/om[g]-1)
				tri[,,g] <- diag(ng[g]/(om[g]*(ng[g]+2*k[g]))*diag(E[,,g]))
			}
		}
		if(q>1){
			for(g in 1:G){
				diag(E[,,g]) <- diag(sg[,,g] - 2*lg[,,g]%*%betag[,,g]%*%sg[,,g] + 
							lg[,,g]%*%thetag[,,g]%*%t(lg[,,g]))
				k[g] <- (ng[g]/2)*((prod(diag(E[,,g]))^(1/p))/om[g]-1)
				tri[,,g] <- diag(ng[g]/(om[g]*(ng[g]+2*k[g]))*diag(E[,,g]))
			}
		}
	}
	if(substring(mod,3,3)=="4"){
		E <- array(0, dim=c(p,p,G))
		k <- rep(0,G)
		if(q==1){
			for(g in 1:G){
				E[,,g] <- sg[,,g]-lg[,,g]%*%t(betag[,,g])%*%sg[,,g]
				k[g] <- (ng[g]/2)*((prod(diag(E[,,g]))^(1/p))/om[g]-1)
				tri[,,g] <- diag(diag(E[,,g])/(om[g]*(1+2*k[g]/ng[g])))
			}
		}
		if(q>1){
			for(g in 1:G){
				E[,,g] <- sg[,,g]-lg[,,g]%*%betag[,,g]%*%sg[,,g]
				k[g] <- (ng[g]/2)*((prod(diag(E[,,g]))^(1/p))/om[g]-1)
				tri[,,g] <- diag(diag(E[,,g])/(om[g]*(1+2*k[g]/ng[g])))
			}
		}
		
#		k <- (prod(diag(E))^(1/p)-sum(ng))/2
#		findum <- diag(E)/(sum(ng) + 2*k)
#		tri[,,] <- diag(findum)
	}
	tri
}
