
#---------------------------------------------------
# summary of jackknife statistic
summary.R2noharm.jackknife <- function(object,logfile=NULL,...){
	# INPUT:
	# object ... object of class jackknife.R2noharm
	#.........................
	if ( ! is.null(logfile) ){ sink( paste0(logfile,".Rout") , split=TRUE ) }	
	NJ <- length(object$u.jackunits)
	cat("Jackknife NOHARM Model with" , NJ , "Jackknife Units\n\n")
	dfr <- object$partable
	dfr[,1] <- gsub( "\\.stat" , "" , paste( dfr[,1] ))
	dfr1 <- dfr
	dfr1[,3] <- round( dfr1[,3] , 4 )
	dfr1[,4] <- round( dfr1[,4] , 4 )
	dfr1[,5] <- round( dfr1[,5] , 4 )	
	print( dfr1 )
	if ( ! is.null(logfile) ){ sink() }		
	invisible(dfr)
	}
#----------------------------------------------------

