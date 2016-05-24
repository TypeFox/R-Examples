mi.anova <- function( mi.res , formula , type = 2 ){
    # INPUT:
    # mi.res  ... mids object (from mice imputation function)
    # formula ... formula for ANOVA model (variable names must be in colnames(mi.list[[1]]), ...
    # converting MICE object to mi.list
    mi.list <- mi.res
    if( class(mi.list) == "mids.1chain" ){	
		mi.list <- mi.list$midsobj
			}
    if( class(mi.list) == "mids" ){
            # number of imputations 
            m <- mi.list$m
            # list of completed datasets
            h1 <- list( rep("", m ))
            for (ii in 1:m){ 
				h1[[ii]] <- as.data.frame( mice::complete( mi.list , ii ) ) 
							}										
            mi.list <- h1
            }
    # converting mi.norm objects
    if (class(mi.res) == "mi.norm" ){ mi.list <- mi.list$imp.data }
	#**** type II sum of squares
	if ( type==2){
		anova.imp0 <- lapply( mi.list , FUN = function(dat){ 
						stats::lm( formula, data = dat ) } )                 
		anova.imp <- lapply( anova.imp0 , FUN = function( obj){ 
					summary( stats::aov(obj)) 
							} )
				}
	#**** type III sum of squares
	if (type==3){
	     Nimp <- length(mi.list)
#		 Nimp <- 1
	     vars <- base::all.vars( as.formula( formula ))[-1]
		 VV <- length(vars)
		 # define contrasts
         ma_contrasts <- as.list(1:VV)		 
		 names(ma_contrasts) <- vars
		 dat <- mi.list[[1]]
         for (vv in 1:VV){	
#			names(ma_contrasts)[[vv]] <- vars[vv] 
			ma_contrasts[[ vars[vv] ]] <- "contr.sum"
			if ( ! is.factor( dat[, vars[vv] ]  ) ){
					ma_contrasts[[ vars[vv] ]] <- NULL
								}
						}					
#		 ma_contrasts <- list( age = contr.sum , hyp = contr.sum )	
		 # estimate linear model				
		 anova.imp0 <- lapply( as.list( 1:Nimp) , FUN = function(ii){				
			dat <- mi.list[[ii]]
			mod1 <- stats::lm( formula , data=dat , contrasts=ma_contrasts)
			return(mod1)
				} )
		 # compute summary
		 anova.imp <- lapply( as.list( 1:Nimp) , FUN = function( ii ){ 
					obj <- anova.imp0[[ii]]
					car::Anova(obj , type=3) 
						} 
							)

				
				}

						

    # number of F statistics to be evaluated
	if (type==2){ 
			FF <- nrow( anova.imp[[1]][[1]] ) - 1
				}
	if (type ==3){
			FF <- nrow(anova.imp[[1]]["F value"])-2							
				}
    anova.imp.inf <- t( sapply( 1:FF , FUN = function(ff){
            micombine.F( sapply( 1:( length(anova.imp) ) , FUN = function(ii){ 
						
						if ( type ==2 ){
						  r1 <- anova.imp[[ii]][[1]]$'F value'[ff] 						
										}
						if ( type ==3 ){
						  r1 <- anova.imp[[ii]]$'F value'[ff+1] 						
										}
						return(r1)			
									} ) ,
                            df1 = ifelse (type==2 , anova.imp[[1]][[1]]$Df[ff] ,
										 anova.imp[[1]]["Df"][ff+1,1]     )  , 							
							display = FALSE )
				} ) )

				
    # ANOVA results
    res <- anova.imp.inf[ , c(3,4,1,2) ]
    res <- matrix( res , ncol = 4 )
    res[,3] <- round( res[,3] , 4 )
    res[,4] <- round( res[,4] , 6 )
    g1 <- rownames( anova.imp[[1]][[1]] )[1:FF] 
	if (type==3){ g1 <-	rownames( anova.imp[[1]] )[1 + 1:FF] }
    rownames(res) <- g1
    res <- data.frame(res)	
		
    # compute eta squared and partial eta squared coefficients
	if (type==2){
		SS <- rowMeans( matrix( unlist( lapply( anova.imp , FUN = function(ll){ 
				ll[[1]][,2] } ) ) , ncol = length(mi.list) )  )
						}
     if (type==3){
		SS <- rowMeans( matrix( unlist( lapply( anova.imp , FUN = function(ll){ 
				l2 <- ll["Sum Sq"][-1,1] 			
				return(l2)
						} ) ) , ncol = length(mi.list) )  )
						}
						
    # calculate (average) R squared
    r.squared <-  sum(SS[ - (FF+1) ]) / sum(SS) 
    res$eta2 <- round( SS[ - ( FF + 1 ) ] / sum( SS ) , 6 )
    res$partial.eta2 <- round( SS[ - (FF+1) ] / ( SS[ - (FF+1) ] + SS[ FF + 1 ] ) , 6 )
	g1	 <- c("F value" , "Pr(>F)" )
#    colnames(res)[3:4] <- colnames( anova.imp[[1]][[1]] )[ c(4,5) ]
    colnames(res)[3:4] <- g1
    colnames(res)[1:2] <- c("df1" , "df2")
    c1 <- colnames(res)
    res <- rbind( res , res[1,] )
    rownames( res)[ nrow(res) ]  <- "Residual"
    res[ nrow(res) , ] <- NA
    res <- data.frame( "SSQ" = SS , res )
    colnames(res)[-1] <- c1
    cat("Univariate ANOVA for Multiply Imputed Data ", paste0("Type " , type , ")" ) , " \n\n")
    cat("lm Formula: ", formula  )
    cat( paste( "\nR^2=" , round(r.squared , 6 ) , sep="") , "\n" )
    cat("..........................................................................\n")
    cat("ANOVA Table \n " )
    print( round( res ,5) )	
    invisible( list( "r.squared" = r.squared , "anova.table" = res , type=type ) )
    }
