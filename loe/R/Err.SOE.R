"Err.LOE" <-
function(X, CM,c=1){
								p <- nrow(X)
								return(Objt.SOE(vecx=as.vector(X), cm=CM,P=p,C=c))
							}
