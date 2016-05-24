###############################################################################
plotDevianceTAM   <- function ( tam.obj , omitUntil = 1, reverse = TRUE ,
			change=TRUE) {
        stopifnot(class(tam.obj) %in% c("tam.mml","tam.mml.2pl","tam.mml.mfr","tam.mml.3pl","tamaan") )
		
		devhistory <- tam.obj$deviance.history
        if(omitUntil>0)  {
				devChange <- devhistory[-c(1:omitUntil),2]
						} else { 
				devChange <- devhistory[,2]
							}
		if ( change ){ 	
			devChange <- diff(devChange) 	
				ylab1 <- "Deviance Change"
					} else {
				ylab1 <- "Deviance"
						}
		
        if(reverse)      {devChange <- -1 *  devChange }
        devChange <- data.frame ( nr = omitUntil + 1:length(devChange), devChange)
		xm        <- ceiling( max(devChange[,1])/10 )*10      
		xt        <- NULL;	for ( i in c( 1:30 ) ) xt <- c ( xt , (xm/10) %% i == 0 )
		xt        <- max ( which ( xt ) )                         
		cex       <- 0.85 - ( length(devChange[,1]) / 1000 )       
		if ( cex < 0.40 ) cex <- 0.40                             
		graphics::plot ( devChange[,1] , devChange[,2] , type = "o" , 
				main = "Deviance Change Plot", xlab = "Iteration" , 
				xlim = c(min(devChange[,1]) ,max(devChange[,1])) ,  xaxp = c(0,xm,xt) , 
				ylab = ylab1 , pch = 20 , cex = cex , lwd = 0.75 )
		graphics::abline( a=0 , b=0 )                                   
		dcr       <- devChange[devChange[,2]<0,]               
        graphics::points( dcr[,1] , dcr[,2] , pch=20, cex = cex , col="red") 
			}
###############################################################################