
###################################################
# create parameter table for GDINA model
gdina.partable <- function(res){

	item <- res$item
	
	
	#**************************
	# item parameters
	dfr <- data.frame( "partype" = "item" , "parindex" = NA , "item" = item$itemno ,
				"item.name" = item$item  )
	dfr$skillclass <- 0
	dfr$varyindex <- dfr$item	
	pa <- paste(item$partype.attr)
	pa <- ifelse(  pa == "" , "Int" , pa )
	dfr$parnames <- paste0( dfr$item.name , "_" , gsub( "-" , "_" , pa ) )
	dfr$value <- item$est
	# fixed vs. free parameters	
	dfr$fixed <- FALSE
	dfr$free <- TRUE
	if ( ! is.null( res$control$delta.fixed ) ){
			z1 <- res$control$delta.fixed
			m1 <- unlist( lapply( z1 , FUN = function(ll){ sum( is.na( ll ) ) > 0 } ) )
			m1 <- which( ! m1 )			
			dfr[ dfr$item %in% m1 , "free"] <- FALSE
			dfr[ dfr$item %in% m1 , "fixed"] <- TRUE
					}
	dfr$rule <- item$rule
	dfr$group <- 0
	dfr$totindex <- NA
	dfr0 <- dfr	
	
	#*************************************
	# skill class distribution
	G <- res$G
	ap0 <- ap <- res$attribute.patt
	if (G==1){
	   ap <- matrix( ap[,1] , ncol=1 )
			 }
	dfr1 <- NULL
	for (gg in 1:G){	
		# gg <- 1		 			 
		ap.names <- paste0("prob_class" , rownames(ap0) , "_group" , gg)	
		L <- nrow(ap)
		dfr <- data.frame( "partype" = rep("probs",L)  , "parindex" = NA , "item" = 0 ,
							  "item.name" = "")
		dfr$skillclass <- 1:L
		dfr$varyindex <- NA
		dfr$parnames <- ap.names
		dfr$value <- ap[,gg]			
		dfr$fixed <- FALSE
		dfr$free <- TRUE
		dfr$rule <- ""
		dfr$group <- gg
		dfr$totindex <- NA
		dfr1 <- rbind( dfr1 , dfr )
					}
	dfr0 <- rbind( dfr0 , dfr1 )

	#*************************************
	# marginal skill distribution
	G <- res$G
	ap0 <- ap <- res$skill.patt
	K <- nrow(ap)
	V <- ncol(ap)
		# gg <- 1		 			 
	ap.names <- paste0("prob_skill" , rownames(ap0) ) #  , "_group" , gg)	
	apnames <- NULL
	l1 <- strsplit( colnames(ap0) , split="prob" , fixed=TRUE )
	l1 <- unlist( lapply( l1 , FUN = function(ll){ substring(ll[2],1,1) } ) )
	l2 <- strsplit( colnames(ap0) , split="group" , fixed=TRUE )
	l2 <- unlist( lapply( l2 , FUN = function(ll){ substring(ll[2],1,1) } ) )
	if (G ==1 ){
		l2 <- rep("1",V)
				}
	for (vv in 1:V){
		apnames <- c( apnames , 
					paste0( ap.names , "_lev" , l1[vv] , "_group" , l2[vv])
						)
					}
	    L <- K*V
		dfr <- data.frame( "partype" = rep("margprobs",L)  , "parindex" = NA , "item" = 0 ,
							  "item.name" = "")
		dfr$skillclass <- 0
		dfr$varyindex <- NA

		dfr$parnames <- apnames
		dfr$value <- as.vector(ap)			
		dfr$fixed <- FALSE
		dfr$free <- TRUE
		dfr$rule <- ""
		dfr$group <- rep( l2 , each=K)
		dfr$totindex <- NA
	dfr0 <- rbind( dfr0 , dfr )

	dfr0$totindex <- seq( 1 , nrow(dfr0) )
		
	return(dfr0)
		}
######################################################		
