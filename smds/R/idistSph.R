"idistSph" <-
function(X,r){
					N <- nrow(X)
					P <- ncol(X)
					D <- as.matrix(dist(X))
					ldist <- D - r%*%t(rep(1,N)) - rep(1,N)%*%t(r)
					ldist[which(ldist<0)] <- 0
					udist <- D + r%*%t(rep(1,N)) + rep(1,N)%*%t(r)
					diag(udist) <- 0
					tmpidist <- array(0,dim=c(2,N,N))
					tmpidist[1,,] <- ldist
					tmpidist[2,,] <- udist
					return(tmpidist)
				}