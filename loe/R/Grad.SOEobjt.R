"Grad.SOEobjt"<-
function(vecx, cm,n, P, C=1){
							N <- n
							NC <- nrow(cm)
							X <- matrix(vecx,nrow=N,ncol=P)
							D<- make.distmat(X)
							igrad <- rep(0,N*P)
							.C("SOEgrad", 
								arg1=as.double(igrad),
								arg2=as.double(X),
								arg3=as.double(D),
								arg4=as.integer(cm),
								arg5=as.double( C ),
								arg6=as.integer( N ),
								arg7=as.integer( P ),
								arg8=as.integer(NC)
							)$arg1
						}