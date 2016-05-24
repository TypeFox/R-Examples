restartlsm <- function(N, D, randomZ, Y, psi2){
	
	lsm <- NULL
	
	if(randomZ){
		
		lsm$lsmEZ <- matrix(rnorm(N * D), ncol = D)
		lsm$lsmVZ <- diag(D)
		
		} else {
			
			if(D %in% 2:3){ # Fruchterman-Reingold
				
				lsm$lsmEZ <- layout.fruchterman.reingold(graph.adjacency(Y), dim = D)
				lsm$lsmEZ <- lsm$lsmEZ / apply(lsm$lsmEZ, 2, sd)
				lsm$lsmVZ <- diag(D)

			} else { # Multidimensional Scaling
				lsm$lsmEZ<-cmdscale(as.dist(1-Y), D)
				lsm$lsmVZ<-diag(D)
		}
	
	}
	
	lsm$xiT <- glm(c(Y)~c(as.matrix(dist(lsm$lsmEZ)^2)))$coeff[1]
	names(lsm$xiT) <- NULL
	lsm$psi2T <- psi2
	
	lsm
	
	}