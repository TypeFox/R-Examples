
#############################################################
# plot p values permutation test for LSEM
plot.lsem.permutationTest <- function( x , type="global" , 
		stattype = "SD" , parindex = NULL , 
		sig_add = TRUE , sig_level=.05 , sig_pch=17 ,
		nonsig_pch = 2 , 	sig_cex = 1 , 
		sig_lab = "p value" , 
		stat_lab = "Test statistic" , 
		moderator_lab = NULL , digits=3 , 
		title=NULL , parlabels = NULL , 
		ask=TRUE , ... ){

	if ( is.null(parindex) ){
		NP <- max( x$parameters$parindex )
		parindex <- 1:NP
							}		
		
	#######################################			
	# global test statistic
    if ( type == "global"){
			
		teststat <- x$teststat[ parindex , , drop=FALSE]
		test_p <- teststat[ , paste0(stattype,"_p") ]
		
		NP <- nrow(teststat)
		labs <- paste(teststat$par)
		if ( ! is.null(parlabels) ){
			labs <- parlabels					
				}
		if (sig_add){
			nc <- nchar(labs)
			NC <- max(nc)
			for (pp in 1:NP){
#			    labs[pp] <- paste0( labs[pp] , rep( " " , NC-nc[pp] ) , collapse="" )
					       }
			labs <- paste0( labs , " (p=" , formatC( test_p , digits=digits , format="f") , ")")
				}
				
	
		main <- title 		
		if ( is.null(title)){
     		main <- paste0( "p values " , stattype)
							}
							
		xlab1 <- sig_lab
		
		h1 <- seq(NP , 1 )
		labs <- labs[ h1  ]
		test_p <- test_p[ h1 ]
		
		graphics::dotchart( test_p , labels=labs , xlim=c(0,1) , pch= sig_pch ,
						xlab=xlab1  , main= main , lwd=1.2)
		graphics::abline(v= sig_level , col=2 , lty=4 , lwd=2)
						
						}	
	#**********************************************************
	if (type == "pointwise"){
	
		   ppt <- x$parameters_pointwise_test

	for (pp in parindex ){	
		graphics::par(mfrow=c(2,1))		   
			# pp <- parindex[2]
			ind.pp <- which( parindex == pp)
			
			x.pp <- ppt[ ppt$parindex == pp , ]
			ylab1 <- sig_lab
			if ( is.null(moderator_lab) ){
				moderator_lab <- x$moderator 
									}
			if ( is.null(parlabels)){						
				t1 <- paste( x.pp$par[1] )
						} else {
				t1 <- parlabels[ind.pp]
							}

			graphics::plot( x.pp$moderator , x.pp$est , xlab= moderator_lab , 
					ylab= stat_lab , main = t1 , type="o" , pch=16)							
			graphics::abline( h=0 , lty=3 , col="gray")
			
			# p value
			pch1 <- ifelse( x.pp$p < sig_level , sig_pch , nonsig_pch )			
			graphics::plot( x.pp$moderator , x.pp$p , xlab= moderator_lab , 
					ylab= ylab1 , main = t1 , type="o" , pch=pch1 , 
					ylim=c(0,1) , cex=sig_cex)
					
		    # lines( spline( modgrid[,1] , y = x.pp$p , n=100 )  )					
			graphics::abline(h= sig_level , col=2 , lty=4 , lwd=2)

			graphics::par(mfrow=c(1,1))				
			graphics::par(ask=ask)			
					}
				
						}
			}
###################################################################			