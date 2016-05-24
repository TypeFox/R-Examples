"Objt.LOE" <-
function(vecx, adm,P,C=1){
								N <- nrow(adm)
								X <- matrix(vecx,nrow=N,ncol=P)
								tmpD <- make.distmat(X)
								.Call("LOEobjt",
									as.double(tmpD),
									adm,
									as.double( C )
								)
							}
