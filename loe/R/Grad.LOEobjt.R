"Grad.LOEobjt" <-
function(vecx,adm, P, C=1){
					N <- nrow(adm)
					X <- matrix(vecx,nrow=N,ncol=P)
					D<- make.distmat(X)
					igrad <- rep(0,N*P)
					.C("LOEgrad",
						arg1=as.double(igrad),
						arg2=as.double(X),
						arg3=as.double(D),
						arg4=as.integer(adm),
						arg5=as.integer(N),
						arg6=as.integer(P),
						arg7=as.double(C)
					)$arg1
				}

