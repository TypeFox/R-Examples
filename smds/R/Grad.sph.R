"Grad.sph" <-
function(vecxr,ldiss,udiss, P){
					N <- nrow(ldiss)
					X <- matrix(vecxr[1:(N*P)],nrow=N,ncol=P)
					r <- vecxr[(N*P)+1:N]
					D <- as.matrix(dist(X))
					IDM <- idistSph(X,r^2)
					ldm <- IDM[1,,]
					udm <- IDM[2,,]
					igrad <- rep(0,N*P+N)
					.C("sphgrad",
						arg1=as.double(igrad),
						arg2=as.double(X),
						arg3=as.double(r),
						arg4=as.double(D),
						arg5=as.double(ldm),
						arg6=as.double(udm),
						arg7=as.double(ldiss),
						arg8=as.double(udiss),
						arg9=as.integer(N),
						arg10=as.integer(P)
					)$arg1
				}

