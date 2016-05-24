"idistBox" <-
function(X,R){
					N <- nrow(X)
					P <- ncol(X)
					tmp <-	.C("bidist",
										arg1 = as.double(X), 
										arg2 = as.double(R), 
										arg3 = as.double(array(0,dim=c(N,N))), 
										arg4 = as.double(array(0,dim=c(N,N))), 
										arg5 = as.integer(N), 
										arg6 = as.integer(P)
									)
					tmpidist <- array(0,dim=c(2,N,N))
					tmpidist[1,,] <- tmp$arg3
					tmpidist[2,,] <- tmp$arg4
					return(tmpidist)
				}