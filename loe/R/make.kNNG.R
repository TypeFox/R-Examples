"make.kNNG" <-
function(DM, k=as.integer( 2*log(nrow(DM)) ), symm =FALSE, weight=FALSE ) {#begin function
					N <- nrow(DM)
					ADM <- matrix(0, N, N)
					#Search kNN point
					if(weight==TRUE){
						for (i in 1:N) {
							nid <- order(DM[i,])
							ADM[ i, nid[2:(k+1)] ] <- DM[ i, nid[2:(k+1)] ]
						}
					}else{
						for (i in 1:N) {
							nid <- order(DM[i,])
							ADM[i,nid[2:(k+1)] ] <- 1
						}
					}
					if(symm==TRUE){
						SADM <- ADM+t(ADM)
						SADM[SADM==2*ADM] <- ADM[SADM==2*ADM]
						ADM <- SADM
					}
					return(ADM)
			}
