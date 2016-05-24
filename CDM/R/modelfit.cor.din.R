
#############################################################
# Model fit for din object
modelfit.cor.din <- function( dinobj , jkunits = 0 ){
    mod <- dinobj
	
	# classes din and gdm
	if ( class(dinobj) %in% c("din","gdina") ){
		data <- as.matrix( mod$data )
		posterior <- mod$posterior
		probs <- mod$pjk
				}
	# class gdm
    if (class(mod) == "gdm"){
		jkunits <- 0
		probs <- aperm( mod$pjk , c(2,3,1) )
		posterior <- mod$posterior
		data <- as.matrix(mod$data)
			}
##			
##   R>  print( str (probs)) 
##    num [1:12, 1:2, 1:6] 0.816 0.722 0.749 0.788 0.847 ...
##   [ items , categ , skills ]
##
##   R>  print( str (posterior)) 
##    num [1:1000, 1:6] 0.986 0.9 0.986 0.986 0.88 ...
##  [ persons , skills ]
#    res <- modelfit.cor( data , posterior , probs )

 
    res <- modelfit.cor2( data , posterior , probs ) 

	#******************************************
	# Jackknife
	HJJ <- sum( abs( jkunits) )
	if ( HJJ > 0 ){
		data <- dinobj$data 
		weights <- dinobj$control$weights
		q.matrix <- dinobj$q.matrix
		guess.init <- dinobj$guess$est
		slip.init <- dinobj$slip$est
		c1 <- dinobj$control  
		N <- nrow(data)		
		if ( length(jkunits) == 1){
		    jkunits <- min( N , jkunits )
			jkunits <- floor( jkunits  * ( 1:N ) / (N+1) ) + 1
							}
		JJ <- length( unique( jkunits ) )
		jkunits <- match( jkunits , unique(jkunits ) )
		ms <- res$modelfit.stat
		ms.jack <- matrix( NA , nrow= nrow(ms) , ncol=JJ )
		rownames(ms.jack) <- rownames(ms)
		cat( paste0("|" , paste( rep("*",20) , collapse="" ) , "|\n|") )
		progressbar_init <- 1:JJ
		progressbar_init <- diff( floor( 20 * ( 1:JJ ) / (JJ+1) ) + 1  )
		progressbar_init <- c(20-sum(progressbar_init), progressbar_init )
		for (jj in 1:JJ){
	#		jj <- 1
			data.jj <- data[ jkunits != jj , ]
			weights.jj <- weights[ jkunits != jj ]  
		#print( paste0("Iteration" , jj ) )			
			#****
			# DINA/DINO model			
			if (class(dinobj)=="din"){
				mod.jj <- din( data=data.jj , q.matrix =q.matrix , 
							skillclasses = c1$skillclasses , conv.crit = c1$conv.crit, 
							dev.crit = c1$dev.crit , maxit = c1$maxit ,
							constraint.guess = c1$constraint.guess , constraint.slip = c1$constraint.slip ,
							guess.init = guess.init , slip.init = slip.init,
							guess.equal = c1$guess.equal , slip.equal = c1$slip.equal , 
							zeroprob.skillclasses = c1$zeroprob.skillclasses , 
							weights = weights.jj ,  rule = c1$rule , 
							wgt.overrelax = c1$wgt.overrelax , 
							wgtest.overrelax = c1$wgtest.overrelax , 
							param.history = FALSE , 
							progress = FALSE )
													}
			#****
			# GDINA model
			if (class(dinobj)=="gdina"){
				mod.jj<- gdina( data=data.jj, q.matrix, skillclasses=c1$skillclasses , 
						conv.crit = c1$conv.crit , 
						dev.crit = c1$dev.crit , maxit = c1$maxit ,
						linkfct = c1$linkfct , Mj = c1$Mj , 
						group = c1$group[ jkunits != jj  ] , 
						method = c1$method , 
						delta.designmatrix = c1$delta.designmatrix , 
						delta.basispar.lower = c1$delta.basispar.lower , 
						delta.basispar.upper = c1$delta.basispar.upper , 					
#						delta.basispar.init = unlist(dinobj$delta)  , 
						zeroprob.skillclasses = c1$zeroprob.skillclasses , 
						reduced.skillspace= c1$reduced.skillspace , 
						HOGDINA = c1$HOGDINA , 
						Z.skillspace = c1$Z.skillspace , 
						weights = weights.jj ,  rule = c1$rule , 
						progress = FALSE ,	progress.item = FALSE )			
					}						
			#*** evaluate model fit
#			f1jj <- modelfit.cor( data=mod.jj$data , posterior=mod.jj$posterior , 
#							probs=mod.jj$pjk )		
			f1jj <- modelfit.cor2( data=mod.jj$data , posterior=mod.jj$posterior , 
							probs=mod.jj$pjk )									
			ms.jack[,jj] <- f1jj$modelfit.stat[,1]
			if ( progressbar_init[jj] == 1 ){ cat("-") ; utils::flush.console() }
					}
		cat("|\n")
		res$modelfit.stat.jack <- ms.jack
		# pseudo values
		ms1 <- ms[,1]
		psx <- ms1 + ( JJ-1 )* ( ms1 - ms.jack )
		# jackknife estimate
		ms$jkunits <- JJ	
		ms$jk_est <- rowMeans( psx )
		ms$jk_se <- sqrt( rowSums( ( psx - ms$jk_est )^2 ) / (JJ-1 ) / JJ  )
		ms$jk_se <- sqrt( apply( ms.jack , 1 , FUN = function(ll){ 
				sum( ( ll - mean(ll) )^2 ) } ) * (JJ-1) / JJ )
		ms$est_low <- ms$jk_est - 1.96 * ms$jk_se
		ms$est_upp <- ms$jk_est + 1.96 * ms$jk_se
		res$modelfit.stat <- ms	
			}
	###**** output
	class(res) <- "modelfit.cor.din"
	return(res)	
                        }
##############################################################
# summary
summary.modelfit.cor.din <- function( object , ... ){	
	cat("Test of Global Model Fit\n")
	obji <- object$modelfit.test
	for (vv in seq(2,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)
	cat("\nFit Statistics\n")
	obji <- object$modelfit.stat
	for (vv in seq(1,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)		
		}
#################################################################	

