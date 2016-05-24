omupdate <-
function(mod,q,G,yg,p,sg,lg,betag,thetag,pig,om,tri){
	if(substring(mod,3,3)=="1"){
		if(q==1){
			for(g in 1:G){
				om[g] <- (1/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
								2*lg[,,g]%*%t(betag[,,g])%*%sg[,,g] +
								lg[,,g]%*%t(thetag[,,g])%*%t(lg[,,g]))))
			}
		}
		if(q>1){
			for(g in 1:G){
				om[g] <- (1/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
								2*lg[,,g]%*%betag[,,g]%*%sg[,,g] +
								lg[,,g]%*%thetag[,,g]%*%t(lg[,,g]))))
			}
		}													 
	}
	if(substring(mod,3,3)=="2"){
		if(q==1){
			for(g in 1:G){
				om[g] <- (1/p)*sum(diag(solve(tri[,,g]) %*% sg[,,g] - 
										solve(tri[,,g]) %*% lg[,,g] %*% 
										t(betag[,,g]) %*% sg[,,g])) 
			}
		}
		if(q>1){
			for(g in 1:G){
				om[g] <- (1/p)*sum(diag(solve(tri[,,g]) %*% sg[,,g] - 
										solve(tri[,,g]) %*% lg[,,g] %*% 
										betag[,,g] %*% sg[,,g])) 
			}
		}
	}
	if(substring(mod,3,3)=="3"){
		du <- 0
		if(q==1){
			for(g in 1:G){
				du <- du + (pig[g]/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
							 2*lg[,,g]%*%t(betag[,,g])%*%sg[,,g] +
							 lg[,,g]%*%t(thetag[,,g])%*%t(lg[,,g]))))
			}
			om[] <- du
		}
		if(q>1){
			for(g in 1:G){
				du <- du + (pig[g]/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
							2*lg[,,g]%*%betag[,,g]%*%sg[,,g] +
							lg[,,g]%*%thetag[,,g]%*%t(lg[,,g]))))
			}
			om[] <- du
		}
	}
	if(substring(mod,3,3)=="4"){
		du <- 0
		if(q==1){
			for(g in 1:G){
				du <- du + (pig[g]/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
								lg[,,g]%*%t(betag[,,g])%*%sg[,,g])))
			}
			om[] <- du
		}
		if(q>1){
			for(g in 1:G){
				du <- du + (pig[g]/p)*sum(diag(solve(tri[,,g]) %*% (sg[,,g] - 
								lg[,,g]%*%betag[,,g]%*%sg[,,g])))
			}
			om[] <- du
		}
	}
	om
}
