lgupdate <-
function(mod,p,q,G,ng,yginv,sg,betag,thetag,om,tri,lg){

	if(substring(mod,1,1)=="C"){
		if(substring(mod,2,2)=="U"){
			if(substring(mod,3,3)=="C"){
				#### CUC MODEL
				lgdum <- matrix(0,p,q)
				twodum <- matrix(0,q,q)
				if(q==1){
					for(g in 1:G){
						lgdum <- lgdum + (ng[g]*yginv[1,1,g] * sg[,,g] %*% betag[,,g]) 
						twodum <- twodum + ng[g]*yginv[1,1,g] * thetag[,,g]
					}
				}
				if(q>1){
					for(g in 1:G){
						lgdum <- lgdum + (ng[g]*yginv[1,1,g] * sg[,,g] %*% t(betag[,,g])) 
						twodum <- twodum + (ng[g]*yginv[1,1,g] * thetag[,,g])
					}
			
				}
				for(g in 1:G){
					lg[,,g] <- lgdum %*% solve(twodum)
				}
			}
			if(substring(mod,3,3)=="U"){
				### CUU MODEL
				lgsum <- matrix(0,p,q)
				lgdum <- matrix(0,p,q)
				if(q>1){
					for(g in 1:G){
						lgsum <- lgsum + ng[g]*yginv[,,g] %*% sg[,,g] %*% t(betag[,,g]) 
					}
				}
				if(q==1){
					for(g in 1:G){
						lgsum <- lgsum + ng[g]*yginv[,,g] %*% sg[,,g] %*% betag[,,g] 
					}
				}
				for(i in 1:p){
					thetsum <- matrix(0,q,q)
					for(g in 1:G){
						thetsum <- thetsum + (ng[g]*yginv[i,i,g]) * thetag[,,g]
					}
					lgdum[i,] <- lgsum[i,] %*% solve(thetsum)
				}
				for(g in 1:G){
					lg[,,g] <- lgdum
				}
			}
		}
		if(substring(mod,2,2)=="C"){
			### CCC, CCU, Mt2, Mt4 models
			for(g in 1:G){
				if(q==1){
					lg[,,g] <- sg[,,g] %*% betag[,,g] %*% solve(thetag[,,g])
				}
				if(q>1){
					lg[,,g] <- sg[,,g] %*% t(betag[,,g]) %*% solve(thetag[,,g])
				}
			}
		}	
	}
	if(substring(mod,1,1)=="U"){
		for(g in 1:G){
			if(q==1){
				lg[,,g] <- sg[,,g] %*% betag[,,g] %*% solve(thetag[,,g])
			}
			if(q>1){
				lg[,,g] <- sg[,,g] %*% t(betag[,,g]) %*% solve(thetag[,,g])
			}
		}
	}
	if(substring(mod,3,3)=="2"|substring(mod,3,3)=="4"){
		### Mt2, Mt4 models
		for(g in 1:G){
			if(q==1){
				lg[,,g] <- sg[,,g] %*% betag[,,g] %*% solve(thetag[,,g])
			}
			if(q>1){
				lg[,,g] <- sg[,,g] %*% t(betag[,,g]) %*% solve(thetag[,,g])
			}
		}
	}
	if(substring(mod,3,3)=="1"){
		### Mt1 model
		lgdum <- lgsum <- matrix(0,p,q)
		if(q>1){
			for(g in 1:G){
				lgsum <- lgsum + (ng[g]/om[g]) * sg[,,g] %*% t(betag[,,g]) 
			}
		}
		if(q==1){
			for(g in 1:G){
				lgsum <- lgsum + (ng[g]/om[g]) * sg[,,g] %*% betag[,,g] 
			}
		}
		thetsum <- matrix(0,q,q)
		for(g in 1:G){
			thetsum <- thetsum + (ng[g]/om[g]) * thetag[,,g]
		}
		lgdum <- lgsum %*% solve(thetsum)
		lg[,,] <- lgdum
	}
	if(substring(mod,3,3)=="3"){
		### Mt3 model
		lgsum <- matrix(0,p,q)
		lgdum <- matrix(0,p,q)
		if(q>1){
			for(g in 1:G){
				lgsum <- lgsum + ng[g]*solve(tri[,,g]) %*% sg[,,g] %*% t(betag[,,g]) 
			}
			for(i in 1:p){
				thetsum <- matrix(0,q,q)
				for(g in 1:G){
					thetsum <- thetsum + ng[g]*(1/tri[i,i,g]) * thetag[,,g]
				}
				lgdum[i,] <- lgsum[i,] %*% solve(thetsum)
			}
		}
		if(q==1){
			for(g in 1:G){
				lgsum <- lgsum + ng[g]*solve(tri[,,g]) %*% sg[,,g] %*% betag[,,g] 
			}
			for(i in 1:p){
				thetsum <- 0
				for(g in 1:G){
					#thetsum <- thetsum + ng[g]*solve(tri[i,i,g]) %*% t(thetag[,,g])
				  thetsum <- thetsum + ng[g]*(1/tri[i,i,g]) %*% t(thetag[,,g])
				}
				lgdum[i,] <- lgsum[i,] %*% solve(thetsum)
				#lgdum[i,] <- lgsum[i,] %*% (1/thetsum)
			}
		}
		
		lg[,,] <- lgdum
	}
	lg
}
