widals.snow <-
function(j, rm.ndx, Z, Hs, Ht, Hst.ls, locs, lags, b.lag, cv=0, geodesic=FALSE, wrap.around=NULL, GP.mx, stnd.d=FALSE, ltco=-16) {
    
    xenvr <- as.environment(1)
    
	tau <- nrow(Z)
	n <- ncol(Z)
	k <- length(lags)
	
	rho <- GP.mx[ j, 1 ]
	reg <- GP.mx[ j, 2 ]
	alpha <- GP.mx[ j, 3 ]
	beta <- GP.mx[ j, 4 ]
	flatten <- GP.mx[ j, 5 ]
	
    locs.3D <- cbind( locs, rep(0, n) )
	locs.long.3D <- cbind( rep(locs[ ,1], k), rep(locs[ ,2], k)  ,   beta*rep( lags, each=n ) )
	
	z.lags.vec <- rep( lags, each=n )
	
#Z.obsAtLags <- als.prepare(NULL, Z, lags, Z.na=NULL, tt.rng=1:nrow(Z))$X.sub #### tt.rng not used
	
    use.Hst.ls <- !is.null(Hst.ls)
    
	
	if( cv <= 0 ) {
		
        Y.als <- H.als.b(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, rho=rho, reg=reg, b.lag=b.lag, Hs0=NULL, Ht0=NULL, Hst0.ls=NULL)$Z.hat
		assign( "Y.als", Y.als, pos=xenvr )
        
#yalsp <- als.prepare(NULL, Y.als, lags, Z.na=NULL, tt.rng=1:nrow(Y.als)) #### tt.rng not used
#rm( Y.als )
		
        Z.delta <- Z - Y.als
		Z.delta <- Z.clean.up( Z.delta )
        
#rm(Z.obsAtLags)
#Y0.als <- yalsp$Y.sub
#rm(yalsp)
		
        if( cv == -1 ) { ##################### zero only self ref sites at lag 0
            self.refs <- (which( lags == 0 ) - 1)*n
        }
        if( cv == -2 ) { ##################### zero all self ref sites across all lags
            self.refs <- I(0:(k-1)) * n
        }
        

        
 		Z.adj <- crispify( locs1=locs.3D, locs2=locs.long.3D, Z.delta=Z.delta, z.lags.vec=z.lags.vec, geodesic=geodesic, 
						  alpha=alpha, flatten=flatten, self.refs=self.refs, lags=lags, stnd.d=stnd.d, log10cutoff=ltco )
		
		
##		mean func pred  +  ( obs @ lags  -  mean func pred @ lags )  %*%  Weights
        
#Y.idw <- Y0.als  +   Z.adj
		
		Z.wid <- Y.als  +   Z.adj

	}
	
	if( cv == 1 | cv == 2 ) { #### cv=1 drops only site at lag 0, cv=2 drops site for all time lags #### only uses cv = 2
		
        Z.wid <- matrix(NA, tau, n )
        
		loc.0 <- which( lags == 0 )
		
		for( kk in 1:length(rm.ndx) ) {
			ii <- rm.ndx[[kk]]
			if( cv == 1 ) { drop.ndx <- (loc.0-1)*n + ii } else { drop.ndx <- ( rep( 1:k, each=length(ii) ) - 1 )*n + ii }
			
			
            if( use.Hst.ls ) { 
                red.Hst.ls <- list()
                for(i in 1:tau) { 
                    red.Hst.ls[[i]] <- Hst.ls[[i]][ -ii, , drop=FALSE]  
                } 
            } else { 
                red.Hst.ls <- NULL
            }
            
            Y.als <- H.als.b(Z=Z[ , -ii, drop=FALSE], Hs=Hs[ -ii, , drop=FALSE], Ht=Ht, Hst.ls=red.Hst.ls, 
							 rho=rho, reg=reg, b.lag=b.lag, Hs0=Hs, Ht0=Ht, Hst0.ls=Hst.ls)$Z0.hat
			
			
#yalsp <- als.prepare(NULL, Y.als, lags, is.na(Y.als), tt.rng=1:nrow(Y.als)) #### tt.rng not used
#           rm(Y.als)
            
            Z.delta.drop <- Z[ , -ii, drop=FALSE] - Y.als[ , -ii, drop=FALSE]
			Z.delta.drop <- Z.clean.up( Z.delta.drop )
			
			z.lags.vec.drop <- z.lags.vec[ -drop.ndx ]
            
#Y0.als <- yalsp$Y.sub[ , ii]
#           rm(yalsp)
            

            Z.adj <- crispify( locs1=locs.3D[ ii, , drop=FALSE ], locs2=locs.long.3D[ -drop.ndx , , drop=FALSE ],
							  Z.delta=Z.delta.drop, z.lags.vec=z.lags.vec.drop, geodesic=geodesic, 
							  alpha=alpha, flatten=flatten, self.refs=c(-1), lags=lags, stnd.d=stnd.d, log10cutoff=ltco )
            

			
			Z.wid[ , ii] <- Y.als[ , ii ]  +   Z.adj
			
			
		}
	}
	
	return(Z.wid)
}
