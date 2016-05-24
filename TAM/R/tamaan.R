
##########################################################
# tamaan function
tamaan <- function( tammodel , resp , tam.method=NULL, 
             control=list() , doparse=TRUE , ... ){
    #******************************
	# process syntax with tamaanify
	cl <- match.call()
	cl1 <- paste(cl)
	# cl1con <- cl1["con"]
    incr.fac <- FALSE	
	if ( length( grep("increment.factor" , cl1 ) ) > 0 ){
		incr.fac <- TRUE 
			}	

	s0 <- Sys.time()	
	res0 <- tamaanify( tammodel=tammodel , resp=resp , tam.method=tam.method ,
			doparse=doparse )
	
			
	anal.list <- res0$ANALYSIS.list
	resp <- res0$resp
	#*** attach control elements (see tam.mml)
    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 , QMC=TRUE , 
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 4 ,            
                 maxiter = 1000 , max.increment = 1 , 
                 min.variance = .001 , progress = TRUE , ridge=0 ,
                 seed = NULL , xsi.start0=FALSE , increment.factor=1 , fac.oldxsi=0)
		if ( anal.list$type %in% c("LCA","OLCA") ){
			#con$increment.factor <- 1.05 
			con$increment.factor <- 1.01
						}
		if ( anal.list$type %in% c("TRAIT") ){
			# con$increment.factor <- 1.02 
			con$increment.factor <- 1.01
						}
		if ( anal.list$type %in% c("MIXTURE") ){
#			con$increment.factor <- 1.035
			con$increment.factor <- 1.01
						}						
		if ( anal.list$type %in% c("LOCLCA") ){
			con$increment.factor <- 1.01
						}
	if ( incr.fac ){
		    con$increment.factor <- control$increment.factor 	
					}
					
	#a0 <- Sys.time()			   
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
        }
		
	#******************************
	# tam.mml
    if ( res0$method == "tam.mml" ){
		res <- tam.mml( resp=res0$resp , A=res0$A , xsi.fixed=res0$xsi.fixed ,
					Q=res0$Q , variance.fixed=res0$variance.fixed , 
					control=con , ... )
		res$tamaan.method <- "tam.mml"			
						}
						
	#******************************
	# tam.mml.2pl
    if ( res0$method == "tam.mml.2pl" ){	

		res <- tam.mml.2pl( resp=res0$resp , A=res0$A , xsi.fixed=res0$xsi.fixed ,
					Q=res0$Q , variance.fixed=res0$variance.fixed ,
					B.fixed=res0$B.fixed , est.variance=res0$est.variance,
					control=con, ... )
		res$tamaan.method <- "tam.mml.2pl"			
						}	
	#******************************
	# 3PL: latent class analysis
	if ( ( res0$method == "tam.mml.3pl" ) & ( anal.list$type == "LCA" ) ){
		res <- tamaan.3pl.lca( res0=res0 , anal.list=anal.list , con=con , ... )
								}
	#***********************************

	#******************************
	# 3PL: ordered latent class analysis
	if ( ( res0$method == "tam.mml.3pl" ) & ( anal.list$type == "OLCA" ) ){
		res <- tamaan.3pl.olca( res0=res0 , anal.list=anal.list , con=con , ... )
								}
	#***********************************

	#******************************
	# 3PL: trait model
	if ( ( res0$method == "tam.mml.3pl" ) & ( anal.list$type == "TRAIT" ) ){	
		res <- tamaan.3pl.trait( res0=res0 , anal.list=anal.list , con=con , ... )
								}
	#***********************************
	
	#******************************
	# 3PL: located latent class analysis
	if ( ( res0$method == "tam.mml.3pl" ) & ( anal.list$type == "LOCLCA" ) ){	
		res <- tamaan.3pl.loclca( res0=res0 , anal.list=anal.list , con=con , ... )
								}
	#***********************************

	#******************************
	# 3PL: mixture distribution models
	if ( ( res0$method == "tam.mml.3pl" ) & ( anal.list$type == "MIXTURE" ) ){	
		res <- tamaan.3pl.mixture( res0=res0 , anal.list=anal.list , con=con , ... )
								}
	#***********************************	
	
	
	s1 <- Sys.time()
	time1 <- c( s0 , s1 , s1-s0)
	res$time <- time1
		
	# add tamaanify object to output
	res0$resp <- NULL
	res$tamaanify <- res0
	
    class(res) <- "tamaan"
	res$CALL <- cl
	
	return(res)
				}
###########################################################


