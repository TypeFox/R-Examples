"obj.sph" <-
function(vecxr,ldiss,udiss, P){
					N <- nrow(ldiss)
					X <- matrix(vecxr[1:(N*P)],nrow=N,ncol=P)
					r <- vecxr[(N*P)+1:N]^2
					D <- dist(X)
					IDM <- idistSph(X,r)
					tmp <- array(0,dim=c(2,N,N))
					tmp[1,,] <- ldiss
					tmp[2,,] <- udiss
					return(sum( (tmp-IDM)^2 )/2)
				}