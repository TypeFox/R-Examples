yginitc <-
function(p,G,sg,lg,sgc,mod,pig){
	yg <- array(0, dim=c(p, p, G))
	ygc <- matrix(0,p,p)
	if(substring(mod,1,1)=="U"){
		if(substring(mod,2,2)=="U"){
			for(g in 1:G){
				yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*% t(lg[,,g])))*diag(p)
			}
		}
		if(substring(mod,2,2)=="C"){
			for(g in 1:G){
				yg[,,g] <- (1/p)*sum(diag(sg[,,g] - lg[,,g] %*% t(lg[,,g])))*diag(p)
				ygc <- ygc + pig[g]*yg[,,g]
			}
			for(g in 1:G){
				yg[,,g] <- ygc
			}
		}
	}
	if(substring(mod,1,1)=="C"){
		if(substring(mod,2,2)=="U"){
			for(g in 1:G){
				yg[,,g] <- (1/p)*sum(diag(sgc - lg[,,g] %*% t(lg[,,g])))*diag(p)
			}
		}
		if(substring(mod,2,2)=="C"){
			for(g in 1:G){
				yg[,,g] <- (1/p)*sum(diag(sgc - lg[,,g] %*% t(lg[,,g])))*diag(p)
				ygc <- ygc + pig[g]*yg[,,g]
			}
			for(g in 1:G){
				yg[,,g] <- ygc
			}
		}
	}
	yg
}
