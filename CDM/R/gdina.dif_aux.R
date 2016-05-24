###############################################
# auxiliary calculations in gdina model for
# detection of differential item functioning
gdina.dif.aux <- function( ocontrol , gg , data ){
	aggr.patt.designmatrix <- ocontrol$aggr.patt.designmatrix
	Aj <- ocontrol$Aj
	Mj <- ocontrol$Mj
	Mj.index <- ocontrol$Mj.index
	linkfct <- ocontrol$linkfct
	rule <- ocontrol$rule
	method <- ocontrol$method
	aggr.attr.patt <- ocontrol$aggr.attr.patt
	IP <- ocontrol$IP
	p.aj.xi <- ocontrol$p.aj.xi[,,gg]
	item.patt.split <- ocontrol$item.patt.split
	resp.patt <- ocontrol$resp.patt
	item.patt.freq <- ocontrol$item.patt.freq[,gg]
	freq.pattern <- ocontrol$freq.pattern
	invM.list <- ocontrol$invM.list
	eps <- eps2 <- 10^(-10)
	J <- length(Mj)
	prob_exp <- varmat.delta <- varmat.palj <- delta <- ndj <- as.list(1:J)
	R.lj <- ocontrol$R.lj.gg[,,gg]
	I.lj <- ocontrol$I.lj.gg[,,gg]
    # calculation of expected counts
    R.ljM <- R.lj %*% aggr.patt.designmatrix
    I.ljM <- I.lj %*% aggr.patt.designmatrix
    for (jj in 1:J){    # begin item
        Ajjj <- Aj[[jj]]
        Mjjj <- Mj[[jj]][[1]]
        Rlj.ast <- R.ljM[ jj, Mj.index[jj,5]:Mj.index[jj,6] ]
        Ilj.ast <- I.ljM[ jj, Mj.index[jj,5]:Mj.index[jj,6] ]
        pjjj <- Rlj.ast / ( Ilj.ast + eps2 )
        if (linkfct == "logit" ){ 
                pjjj[ pjjj > 1-eps ] <- 1 - eps
                pjjj[ pjjj < eps ] <- eps
                pjjj <- stats::qlogis( pjjj ) 
                                }
        #*****
    if (linkfct == "log" ){ 
                pjjj[ pjjj < eps ] <- eps
                pjjj <- log( pjjj )                         
                                }
                                
        Wj <- diag( Ilj.ast )
        
        if ( ( rule[jj] == "GDINA" )| ( method == "ULS" ) ){ 
                invM <- invM.list[[jj]] 
                delta.jj <- invM %*% crossprod(Mjjj ,pjjj)					
                            } else { 
                invM <- solve( crossprod(Mjjj , Wj ) %*% Mjjj + diag( rep( eps2 , ncol(Mjjj) )) )               
                delta.jj <- tcrossprod( invM , Mjjj ) %*% Wj %*% pjjj				
                                }
								
				pjj_exp <- ( Mjjj %*% delta.jj )[,1]
				if ( linkfct == "logit" ){ pjj_exp <- stats::plogis( pjj_exp) }
				if ( linkfct == "log" ){ pjj_exp <- exp( pjj_exp) }	
		# data frame with counts and expected probabilities
		prob_exp[[jj]] <- data.frame( "freq" = Ilj.ast / sum( Ilj.ast) , "prob" = pjj_exp )
		rownames(prob_exp[[jj]]) <- apply( Ajjj , 1 , 
				FUN = function(ll){ paste0("S" , paste0( ll , collapse="") ) } )
        delta[[jj]] <- delta.jj[,1]

		#***********************************
        #*****
        # variance matrix
        PAJXI <-  p.aj.xi
                Ajjj <- Aj[[jj]]
                Mjjj <- Mj[[jj]][[1]]
                Mjj2 <- Mj[[jj]][[2]]
                apjj <- aggr.attr.patt[[jj]] 
                M1 <- max( apjj )
                p.ajast.xi <- matrix( 0 , nrow=IP , ncol = M1 )
                for (kk in 1:M1){
                    pg1 <-  PAJXI[ , apjj == kk  ]                  
                    if ( is.vector(pg1)){ 
                                p.ajast.xi[,kk] <- pg1 
                                    } else {
                                p.ajast.xi[,kk] <- rowSums( pg1 ) 
                                        }
                                }   
                Rlj.ast <- stats::aggregate( R.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
                Ilj.ast <- stats::aggregate( I.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
                pjjj <- Rlj.ast[,2] / Ilj.ast[,2]       
                pjjjM <- outer( rep(1,IP) , pjjj ) + 10^(-20)
                
                nM <- ncol(pjjjM) 
                x1 <- outer( item.patt.split[,jj] , rep(1,nM) )
                r1 <- outer( resp.patt[,jj] * item.patt.freq , rep(1,ncol(pjjjM) ) )
                # Formula (17) for calculating the standard error   
                mat.jj <- p.ajast.xi * ( x1 - pjjjM) / ( pjjjM * ( 1 - pjjjM ) )    
                infomat.jj <- matrix( 0 , nM , nM )
                for (kk1 in 1:nM){
                    for (kk2 in kk1:nM){ 

                        # frequency weights must be taken into account
                        hh1 <- sum( mat.jj[,kk1] * mat.jj[,kk2] * item.patt.freq * 
                                            resp.patt[,jj] * item.patt.split[,jj] )
                        infomat.jj[kk2,kk1] <- infomat.jj[kk1,kk2] <-  hh1
                                        }
                                    }

                try( a1 <- solve( infomat.jj + diag( eps2 , ncol(infomat.jj) ) ) )
                if ( is.null(a1)){ 
                        cat( "Item" , colnames(data)[jj] , "Singular item parameter covariance matrix\n")
                        a1 <- NA*infomat.jj 
                            }
                varmat.palj[[jj]] <- Ijj <- a1
                Wj <- diag( Ilj.ast[,2] )   

                if ( ( method == "ULS" ) ){             
                    x1 <- t(Mjjj) %*% Mjjj  
                    diag(x1) <- diag(x1) + 10^(-8)                      
                    Wjjj <- solve( x1 ) %*% t(Mjjj)
                                    } else {
                    x1 <- t(Mjjj) %*% Wj %*% Mjjj                                   
                    diag(x1) <- diag(x1) + 10^(-8)                  
                    Wjjj <- solve( x1 ) %*% t(Mjjj) %*% Wj
                                            }
                if ( linkfct == "logit" ){
                    pjjj.link <- 1 / ( ( pjjj * ( 1 - pjjj ) ) + eps2 )
                    pjjj.link <- diag( pjjj.link )
                    Wjjj <- Wjjj %*% pjjj.link
                        }
                if ( linkfct == "log" ){
                    pjjj.link <- 1 /  ( pjjj  + eps2 )
                    pjjj.link <- diag( pjjj.link )
                    Wjjj <- Wjjj %*% pjjj.link
                        }
                varmat.delta[[jj]] <- Wjjj %*% Ijj %*% t(Wjjj) 
                
        ndj[[jj]] <- length( delta[[jj]] )
                }   # end jj
	res <- list( "delta" = delta , "varmat.delta" = varmat.delta ,
		"ndj" = ndj , "prob_exp" = prob_exp )
    return(res)
		}
#####################################################################
