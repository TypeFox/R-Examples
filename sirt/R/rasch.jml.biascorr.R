

#########################################################################################
# analytical bias correction of JML estimated Rasch model
# rasch.jml.biascorr <- function( data , data.resp , theta , b , adjmin=.5 ){
rasch.jml.biascorr <- function( jmlobj , itemfac = NULL # , adjmin=.5 
				){
	mod <- jmlobj
	m1 <- mod$data.proc
	data <- m1$data
	data.resp <- m1$data.resp
	theta <- m1$theta
	b <- mod$item$itemdiff

    I <- length(b)
	I1 <- mean( rowSums( data.resp ) )
	if (! is.null( itemfac) ){ I1 <- itemfac }
    N <- nrow(data)
    p_it <- stats::plogis( outer( theta , b , "-" ) )
    q_it <- 1 - p_it
    u_it <- - data + p_it   # d/d(theta) f
    v_it <- - u_it      # d/d(alpha) f
    # auxiliary calculations
    u_it_theta <- - p_it * q_it
    u_it_alpha <- p_it * q_it
    v_it_alpha <- - p_it * q_it
    u_it_alpha2 <- p_it * q_it^2 - p_it^2 * q_it
    v_it_alpha2 <- - p_it * q_it^2 + p_it^2 * q_it
    
    # terms
    s2 <- rowWgtMeans( v_it_alpha , data.resp )
    t00 <- rowWgtMeans( v_it^2 , data.resp )
    t1z <- - 1 / t00 
    t0 <- - t00 / s2
    t1b <- rowWgtMeans( v_it * v_it_alpha , data.resp )
    t2b <- rowWgtMeans( v_it_alpha2 , data.resp )
    t2c2 <- t1c2 <- t2c2 <- s2
    t2z <- 1/2 / t1c2
    
    b.bias2 <- b.bias1 <- rep(NA , I )
    M1 <- matrix( 1 , nrow=N,ncol=I)
    
    for (ii in 1:I){
        # ii <- 1
        M1.ii <- M1
        M1.ii[,ii] <- 0
        M1.ii <- M1 - M1.ii
        data.resp.ii <- data.resp[,ii]
        
        #******
        # bias correction method 1
        s1 <- rowWgtMeans( u_it_alpha * M1.ii , data.resp )
        s0 <- rowWgtMeans( u_it_theta * M1.ii , data.resp )
        I_i <- sum( - ( s0 - s1^2 / s2 ) * data.resp.ii / sum(data.resp.ii) )    
        t1a <- rowWgtMeans( v_it * u_it_alpha * M1.ii , data.resp )
        t2c1 <- t1c1 <- rowWgtMeans( u_it_alpha* M1.ii , data.resp )
        t2a <- rowWgtMeans( u_it_alpha2* M1.ii , data.resp )    
        # term
        btheta <- t0 * ( - t1z * ( t1a - t1b * t1c1 / t1c2 )
                        - t2z * ( t2a - t2b  * t2c1 / t2c2 ) )            
        btheta <- sum( btheta * data.resp.ii / sum(data.resp.ii ) )
        b.bias1[ii] <- btheta / I_i
        
        #******
        # bias correction method 2
        a1 <- rowWgtMeans( u_it^2 * M1.ii , data.resp )
        a4 <- a2 <- rowWgtMeans( u_it * v_it * M1.ii , data.resp )
        a3 <- rowWgtMeans( v_it^2  , data.resp )
        I_i <- sum( - ( a1 - a2 * a3 * a4 ) )
        U_it <- u_it * M1.ii - a2 / a3 * v_it
        V_2it <- v_it^2 + v_it_alpha
        t1 <- rowSums( U_it * V_2it  * data.resp )
        t2 <- 2*rowSums( v_it_alpha  * data.resp )
        b.bias2[ii] <- - sum( t1/t2  ) / I_i
        
            }
    
    b.analytcorr1 <- b - b.bias1 / I1
    b.analytcorr2 <- b - b.bias2 / I1
    
#    ind <- which( abs(b) > adjmin )
#    b.analytcorr1b <- mean( ( b.analytcorr1a / b )[ind] )*b
#    b.analytcorr2b <- mean( ( b.analytcorr2a / b )[ind] )*b  
    
    dfr <- data.frame( "b.JML"=b , "b.JMLcorr" = (I1-1)/I1*b , 
		b.analytcorr1 , # b.analytcorr1b ,
        b.analytcorr2  #, b.analytcorr2b 
		)
    
    res <- list( "b.biascorr" = dfr , "b.bias1" = b.bias1 , "b.bias2" = b.bias2 ,
			"itemfac" = I1 )
    return(res)
        }
#########################################################################################
# auxiliary function for bias correction
rowWgtMeans <- function( M , W ){
    rowSums( M*W ) / rowSums(W)
            }
#........................................................................................
