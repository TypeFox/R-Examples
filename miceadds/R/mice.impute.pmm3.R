mice.impute.pmm3 <- function (y, ry, x, donors=3 , noise = 10^5 , 
		ridge = 10^(-5) , ...){
    x <- cbind(1, as.matrix(x))	
	Nmis <- sum(!ry)
# a0 <- Sys.time()	
#    parm <- .norm.draw(y, ry, x, ...)
# cat("norm draw") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
    parm <- .norm.draw3(y, ry, x, ridge=ridge ,  ...)		
# cat("norm draw3") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
	yobs <- y[ ry ]
#cat("calc yhatobs and yhatmis") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
	# define distance matrix
	M1 <- matrix( yhatobs[,1] , nrow=sum(!ry) , ncol=sum(ry) , byrow=TRUE )	
	M2 <- matrix( yhatmis[,1] , nrow=sum(!ry) , ncol=sum(ry) )
	disty <- abs( M2 - M1 )
#cat("define distance matrix") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
	donor.ind <- rowKSmallest2.sirt(matr=disty , K=donors )$smallind
#cat("rowMins") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
	N1 <- nrow(disty)
	# sampled indices
	ind.sample <- sample( 1:donors , N1 , replace=TRUE )
	# select index
	res1 <- donor.ind[ cbind( 1:N1 , ind.sample) ]	
	# search for imputed values
	imp <- yobs[ res1 ]
# cat("rest") ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1			
# stop("here")
    return(imp)
	}

	
	