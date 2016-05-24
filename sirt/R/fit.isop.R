
############################################################
# Fit ISOP Model
fit.isop <- function( freq.correct , wgt , conv = .0001 , 
		      maxit=100 , progress=TRUE , calc.ll=TRUE){
    #****
    # initializations
    # monotone regression
    M1 <- as.matrix(freq.correct)
	I <- ncol(M1)
	PP <- nrow(M1)
    wgt <- as.matrix(wgt)	
    # isotonic row regression
    res <- fit.isotonic.rows(X=M1 , wgt=wgt )
    RX <- res$RX
    IR <- res$IR
    # isotonic column regression
    res <- fit.isotonic.cols(X=RX , wgt=wgt )
    CX <- res$CX
    IC <- res$IC    
    X <- M1
    deviation <- 1
    iter <- 0
	if (progress){ cat("\n*******ISOP Model***********\n") }
	#****
    # ISOP algorithm
    while( ( iter < maxit) & ( deviation > conv ) ){    
        Xold <- X
        # isotonic rows
        res <- fit.isotonic.rows(X=X-IC , wgt=wgt )
        RX <- res$RX
        IR <- res$IR    
        # isotonic columns
        res <- fit.isotonic.cols(X=X-IR , wgt=wgt )
        CX <- res$CX
        IC <- res$IC
        # updated X
        X <- Xold - IR - IC
        # calculate deviation
        deviation <- sum( ( X - Xold )^2*wgt )
		iter <- iter + 1
        if (progress){ 
			cat( "Iteration" , iter , "- Deviation =" ,  round( deviation , 6 ) , "\n")
			utils::flush.console()			
					}
                    }  # end algorithm					
    ###############################################
	RR <- nrow(freq.correct)
	CC <- ncol(freq.correct)
	if ( calc.ll ){ 
		CX[ CX > 1 ] <- 1
		CX[ CX < 0 ] <- 0
				}
	dfr0 <- data.frame( "stud.index" = rep(1:RR , CC) , 
				"item.index" = rep(1:CC ,each=RR) , 
				"wgt" = matrix( as.matrix( wgt ) , RR*CC , 1 ),
				"freq" = matrix( as.matrix(freq.correct ) , RR*CC , 1 ) ,
				"freq.fitted" = matrix( as.matrix( CX ) , RR*CC , 1 )
						)
    dfr0 <- dfr0[ order( dfr0$stud.index * 1000 + dfr0$item.index ) , ]						
    # deviation criterion
    wgt1 <- ( wgt / colSums( wgt ) ) / ncol(wgt)
    fit <- sqrt( sum( ( M1-CX  )^2 * wgt1  ) )
	#****
	# calculate likelihood
	ll <- NULL
	if (calc.ll){
		ll <- list( 
			"ll.ind" = .calc.ll.isop( freq.correct , wgt , irtfitted =freq.correct )	
						)
		ll$ll.isop <- .calc.ll.isop( freq.correct , wgt , irtfitted =CX )
		NW <- mean( colSums(wgt) )
		ll$llcase.ind <- ll$ll.ind /NW	
		ll$llcase.isop <- ll$ll.isop /NW	
		}
	# collect item and person scores
	item.sc <- ( seq( 0 , I-1 ) + .5 ) / I
	person.sc <- ( seq( 0 , PP-1 ) + .5 ) / PP
    # output
    res <- list( "fX" = CX , "ResX" = M1 - CX , "fit" = fit ,
		"item.sc"=item.sc , "person.sc" = person.sc , "ll"=ll ,
		"freq.fitted"=dfr0)
    return(res)
        }
############################################################
# calculate likelihood for saturated and fitted models
.calc.ll.isop <- function( freq.correct , wgt , irtfitted , eps=10^(-20) ){
    ll <- sum( wgt * freq.correct * log( irtfitted + eps ) +
            wgt * (1-freq.correct) * log( 1 - irtfitted + eps ) )
    return(ll)
        }
########################################################


#****************************************
# fit isotonic rows
fit.isotonic.rows <- function( X , wgt ){
    M2 <- X
    RR <- nrow(M2)
#    for (rr in 1:RR){
#        M2[rr,] <- monoreg(x=X[rr,], w=wgt[rr,] )$yf
#                }
#    RX <- M2
	RX <- monoreg.rowwise(X,wgt)
    IX <- X - RX
    res <- list( "X" = X , "RX"  = RX , "IR" = IX )
    return(res)
        }
		
#****************************************
# fit isotonic columns
fit.isotonic.cols <- function( X , wgt ){
    M2 <- X
    RR <- ncol(M2)
#    for (rr in 1:RR){
#        M2[,rr] <- monoreg(x=X[,rr], w=wgt[,rr] )$yf
#                }
#    RX <- M2
    RX <- monoreg.colwise(X,wgt)
    IX <- X - RX
    res <- list( "X" = X , "CX"  = RX , "IC" = IX )
    return(res)
        }		
#******************************************