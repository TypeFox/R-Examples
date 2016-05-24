"get.order"<-
function(DM){
						N <- nrow(DM)
						AM <- rep(-1,4*((N*(N-1))^2))
						avec <- 
									.C("getorder",#getorder(int *AM, double *D, int *N)
										arg1 = as.integer(AM),
										arg2 = as.double(DM),
										arg3 = as.integer(N)
									)$arg1
						A <- matrix(avec,
									byrow=TRUE,
									ncol =4
								)
						return(A[1:( length(which(avec>-1))/4 ),])
					}