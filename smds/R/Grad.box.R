"Grad.box" <-
function(vecxr,ldiss,udiss, P){
					N <- nrow(ldiss)
					X <- matrix(vecxr[1:(N*P)],nrow=N,ncol=P)
					R <- vecxr[(N*P)+1:(N*P)]
					igrad <- rep(0,2*N*P)
					.C("boxgrad",
						arg1=as.double(igrad),
						arg2=as.double(X),
						arg3=as.double(R),
						arg4=as.double(ldiss),
						arg5=as.double(udiss),
						arg6=as.integer(N),
						arg7=as.integer(P)
					)$arg1
				}

