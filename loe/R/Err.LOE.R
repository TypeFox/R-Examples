"Err.LOE" <-
function(X, ADM,c=1){
								p <- nrow(X)
								return(Objt.LOE(vecx=as.vector(X), adm=ADM,P=p,C=c))
							}
