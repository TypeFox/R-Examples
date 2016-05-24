yginitu <-
function(p,G,sg,lg,mod,pig,dumg){
	yg <- array(0, dim=c(p, p, G))
	ygc <- matrix(0,p,p)
	if(substring(mod,1,3)=="CUU"){
		for(g in 1:G){
			diag(yg[,,g]) <- diag(sg[,,g] - dumg[,,g] %*% t(dumg[,,g]))
		}
	}
	else{
		if(substring(mod,2,2)=="U"){
			for(g in 1:G){
				diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% t(lg[,,g]))
			}
		}
		if(substring(mod,2,2)=="C"){
			for(g in 1:G){
				diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% t(lg[,,g]))
				ygc <- ygc + pig[g]*yg[,,g]
			}
			for(g in 1:G){
				yg[,,g] <- ygc
			}
		}
	}
	if(substring(mod,1,1)=="M"){
		if(substring(mod,3,3)=="1"|substring(mod,3,3)=="3"){
			for(g in 1:G){
				diag(yg[,,g]) <- diag(sg[,,g] - dumg[,,g] %*% t(dumg[,,g]))
			}
		}
		if(substring(mod,3,3)=="2"|substring(mod,3,3)=="4"){
			for(g in 1:G){
				diag(yg[,,g]) <- diag(sg[,,g] - lg[,,g] %*% t(lg[,,g]))
			}
		}
	}
	yg
}
