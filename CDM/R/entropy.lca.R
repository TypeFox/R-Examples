
####################################################
# entropy for din, gdina and mcdina objects
entropy.lca <- function( object ){

	# posterior <- object$posterior
	posterior <- object$like
	data <- object$data
	weights <- object$control$weights
	pjk <- object$pjk	
	if ( class(object) == "mcdina" ){
		weights <- object$weights 
		data <- object$dat
		data <- data - 1
		pjk <- object$pik[,,,1]
			}
	
	skillspace <- object$attribute.patt.splitted

	# skillprobs <- object$attribute.patt[,1]
	q.matrix <- object$q.matrix
	data.resp <- 1-is.na(data)
	data[is.na(data) ] <- 0
	
	N <- length(weights)
	weights <- N / sum(weights) * weights
	L <- nrow(skillspace )
	K <- ncol(skillspace)
	I <- ncol(data)
	eps <- 1E-10
	posterior <- posterior / rowSums( posterior )
	
	#*** calculate entropy for whole test
	entropy.total <- 1 +  sum( weights * posterior * log( posterior + eps) ) / N / log(L)
	
	#*** maximum number of skills
	maxskill <- apply( skillspace , 2 , max )
	
	#*** entropy for each skill
	entropy.skill <- rep(0,K)
	for (kk in 1:K){
		# kk <- 1
		Nkk <- maxskill[kk]
		posterior.kk <- matrix(NA , nrow=N , ncol= Nkk+1 )		
		for (vv in 0:Nkk){
			posterior.kk[,vv+1] <- rowSums( posterior[ , skillspace[,kk] == vv ]  )
						}		
		entropy.skill[kk] <- 1 +  sum( weights * posterior.kk * log( posterior.kk + eps) ) / N / log(Nkk+1)
				}

	#******** entropy for each item
	entropyM <- matrix( NA , nrow= I+1 , ncol= K +1 )
	entropyM[1,] <- c( entropy.total , entropy.skill )
	for (ii in 1:I){
		# ii <- 5    
		weights.ii <- weights * data.resp[,ii]
		N.ii <- sum(weights.ii)		
		pjk.ii <- pjk[ii,,]
		pjkM <- pjk.ii[ data[,ii] +1 , ]
				
		posterior.ii <- pjkM 
		posterior.ii <- posterior.ii / rowSums( posterior.ii )
		entropyM[ii+1,1] <- 1 +  sum( weights.ii * posterior.ii * log( posterior.ii + eps) ) / N.ii / log(L)
		for (kk in 1:K){    
			# skill kk and item ii
		    Nkk <- maxskill[kk]
		    posterior.kk <- matrix(NA , nrow=N , ncol= Nkk+1 )		
		    for (vv in 0:Nkk){
			      posterior.kk[,vv+1] <- rowSums( posterior.ii[ , skillspace[,kk] == vv ]  )
						     }			
			entropyM[ii+1,kk+1] <- 1 +  sum( weights * posterior.kk * log( posterior.kk + eps) ) / N / log(Nkk+1)
						}
					  }

	res <- data.frame( "item" = c("test" , colnames(data) ) , entropyM )
	colnames(res)[-1] <- c("entr_test" ,paste0("entr_skill" , 1:K ) )
	res2 <- list( "entropy" = res )
	class(res2) <- "entropy.lca"
	return(res2)
		}
#################################################################################		
# summary S3 method
summary.entropy.lca <- function( object , digits=2, ...){
    obji <- object$entropy
    V <- ncol(obji)
    for (vv in 2:V){
        obji[,vv] <- round( obji[,vv] , digits)
                    }
    rownames(obji) <- NULL
    print(obji)
        }
#####################################################################################