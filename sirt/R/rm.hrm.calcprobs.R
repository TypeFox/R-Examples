
################################################################
# calculate probabilities
.rm.hrm.calcprobs  <- function(  c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				theta.k ,RR , prob.item=NULL , prob.rater = NULL ){
	# calculate probabilities for true ratings
	a <- a.item
	b <- tau.item
	if ( is.null( prob.item ) ){ 
		res <- .rm.pcm.calcprobs( a , b , Qmatrix=Qmatrix , theta.k , I=VV , K , TP )
				} else { 
		res <- prob.item 
				}
	# calculate probabilities for raters
	calc.rater <- FALSE
	if (is.null(prob.rater)){
		calc.rater <- TRUE
						}
# calc.rater <- TRUE						
#	prob.categ <- array( 0 , dim=c(I , K+1 , TP ) )	
#	if (calc.rater){
#		for (ii in 1:I){
	#		ii <- 1	
#			h1[ii,,] <- matrix( c.rater[ii,] , nrow=K+1 , ncol=K,byrow=T) - 
#				         matrix( (0:K) * d.rater[ii] , nrow=K+1 , ncol=K,byrow=F )
#		h1[ii,,] <- - outer( (0:K)*d.rater[ii] , as.vector(c.rater[ii,]) , "-" )
#						}
#		h1 <- plogis( h1 )
#		prob.rater[,1,] <- h1[,,1]
#		for (kk in 1:(K-1)){ prob.rater[,kk+1,] <- h1[,,kk+1] - h1[,,kk] }
#		prob.rater[,K+1,] <- 1-h1[,,K]
#				}
				
	dimA <- c(I , K+1, K+1 )
	res2 <- res[ item.index ,,]	
	dimB <- dim(res2)					
	BM <- matrix( res2 , dimA[1]*dimB[2] , dimB[3] )
	#****	
	# if prob.rater is calculated	
	if (calc.rater){
		res2 <- probraterfct1( crater=c.rater , drater=d.rater , 
				          dimA=dimA,B=BM,dimB=dimB)
		prob.categ <- array( res2$probtotal , dim= c(dimA[c(1,2)],dimB[3]) )
		prob.rater <- array( res2$PRA , dim=dimA )		
					}
				
					
	#***
	# if prob.rater is not calculated
	AM <- matrix( prob.rater , dimA[1]*dimA[2] , dimA[3] )
    if ( ! calc.rater ){
		y <- arraymult1( AM , dimA , BM , dimB )
		prob.categ <- array( y , dim= c(dimA[c(1,2)],dimB[3]) )    
					}
	
	res <- list("prob.total"=prob.categ , "prob.rater"=prob.rater , 
		"prob.item"	= res )
	return(res)	
			}
#############################################################
