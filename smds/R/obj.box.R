"obj.box" <-
function(vecxr,ldiss,udiss, P){
					N <- nrow(ldiss)
					X <- matrix(vecxr[1:(N*P)],nrow=N,ncol=P)
					R <- matrix(vecxr[(N*P)+1:(N*P)]^2,nrow=N)
					IDM <- idistBox(X,R)
					tmp <- array(0,dim=c(2,N,N))
					tmp[1,,] <- ldiss
					tmp[2,,] <- udiss
					return(sum( (tmp-IDM)^2 )/2)
				}