"Objt.SOE"<-
function(vecx, cm,n,P,C=1){
								N <- n
								NC <- nrow(cm)
								X <- matrix(vecx,nrow=N,ncol=P)
								tmpD <- make.distmat(X)
								.C("SOEobjt",
									arg1 = as.double(tmpD),
									arg2 = as.integer(cm),
									arg3 = as.double( C ),
									arg4 = as.integer(N),
									arg5 = as.integer(NC),
									arg6 = as.double(0)
								)$arg6
							}