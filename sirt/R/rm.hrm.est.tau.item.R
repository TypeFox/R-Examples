

			
#####################################################################
.rm.hrm.est.tau.item <- function( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1  , theta.k ,
				msteps, mstepconv , tau.item.fixed , prob.rater ){
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(c.rater)	
	Q0 <- matrix(0,nrow=VV, ncol=K)
	se.tau.item <- Q0
	cat("  M steps tau.item parameter   |")
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		tau.item0 <- tau.item		
		for (kk in 1:K){
	#		kk <- 1
			Q1 <- Q0
			Q1[,kk] <- 1
			r1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR, 
					prob.item=NULL , prob.rater=prob.rater )
			pjk <- r1$prob.total	

			pjk1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item+h*Q1 ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL , prob.rater=prob.rater)$prob.total				
			pjk2 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item-h*Q1 ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL , prob.rater=prob.rater)$prob.total		
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.b.increment , numdiff.parm )					
			increment <- Q1 * matrix( res$increment , nrow=VV , ncol=K)	
			tau.item <- tau.item + increment
			se.tau.item[,kk] <- sqrt(abs(-1/res$d2)	)
					}
		conv1 <- max( abs( tau.item - tau.item0 ) )
		it <- it+1
		cat("-") 
		if (!is.null(tau.item.fixed)){
			tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
								}
			}
	cat(" " , it , "Step(s) \n")
	res <- list("tau.item" = tau.item , "se.tau.item" = se.tau.item , 
			"ll" = sum(res$ll0) , "prob.item"=r1$prob.item )
    return(res)
					}
