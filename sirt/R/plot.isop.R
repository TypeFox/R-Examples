
#############################################################
# plot ISOP item response functions
plot.isop <- function( x , ask=TRUE , ... ){
	if ( x$model == "isop.dich"){
		plotisopdich( x , ask=ask , ... )	
						}
	if ( x$model == "isop.poly"){
		plotisoppoly( x , ask=ask , ... )	
						}						
				}
#############################################################	

##################################################
# plot isop class for polytomous data	
plotisoppoly <- function(x , ask=TRUE , ... ){
		mod3 <- x
		mod <- x
		pers <- mod3$person
		pers <- pers[ order( pers$mpsc ) , ]
		I <- ncol( mod3$dat )
		K <- ncol(x$p.itemcat)
		for (kk in 2:K){
		#kk <- 2
		#********************
		# plot all item response functions
		graphics::par( mfrow=c(2,2))
			# Labels for plots
			kv <- substitute( expression( a ) , list(a=kk) )
			xlabplot <- expression( paste( "Modified P Score " , rho[p] ))
#			ylabplot <- expression(paste(P,"(",X[i], ">=" , kv , "|", rho[p] , ")" ))    
#			ylabplot <- substitute( paste(P,"(",X[i], ">=" , kv , "|", rho[p] , ")" )   ,
#							list(kv=kk-1) )
			ylabplot <- substitute( paste(P,"(",X[i] >= kv , "|", rho[p] , ")" )   ,
							list(kv=kk-1) )
			# saturated IRT model
			p1 <- mod$surv.saturated
			mainplot1 <- paste0( "Survivor Functions Category " , kk - 1 )
			mainplot <- paste0( mainplot1 ,  "\nSaturated Model")			
			plot.all.icc.isop.poly( pers , p1 , xlabplot , ylabplot , mainplot,I,kk)
			# isop
			p1 <- mod$surv.isop
			mainplot <- paste0( mainplot1 ,  "\nISOP Model")
			plot.all.icc.isop.poly( pers , p1 , xlabplot , ylabplot , mainplot,I,kk)
			# adisop
			p1 <- mod$surv.adisop
			mainplot <- paste0( mainplot1 ,  "\nADISOP Model")
			plot.all.icc.isop.poly( pers , p1 , xlabplot , ylabplot , mainplot,I,kk)
			# logistic function
			p1 <- mod$surv.grm
			mainplot <- paste0( mainplot1 ,  "\nGraded Response Model")
			plot.all.icc.isop.poly( pers , p1 , xlabplot , ylabplot , mainplot,I,kk)
		graphics::par( mfrow=c(1,1))
		graphics::par(ask=ask)    
				}
	}
#########################################################################	

##################################################
# plot isop class for dichotomous data	
plotisopdich <- function(x , ask=TRUE , ... ){
		mod3 <- x
		mod <- x
		pers <- mod3$person
		pers <- pers[ order( pers$mpsc ) , ]
		I <- ncol( mod3$dat )
		#********************
		# plot all item response functions
		graphics::par( mfrow=c(2,2))
			# Labels for plots
			xlabplot <- expression( paste( "Modified P Score " , rho[p] ))
			ylabplot <- expression(paste(P,"(",X[i]==1, "|", rho[p] , ")" ))    
			# saturated IRT model
			p1 <- mod$prob.saturated
			mainplot <- paste0( "Item Response Functions \nSaturated Model")
			plot.all.icc.isop.dich( pers , p1 , xlabplot , ylabplot , mainplot,I)
			# isop
			p1 <- mod$prob.isop
			mainplot <- paste0( "Item Response Functions \nISOP Model")
			plot.all.icc.isop.dich( pers , p1 , xlabplot , ylabplot , mainplot,I)
			# adisop
			p1 <- mod$prob.adisop
			mainplot <- paste0( "Item Response Functions \nADISOP Model")
			plot.all.icc.isop.dich( pers , p1 , xlabplot , ylabplot , mainplot,I)
			# logistic function
			p1 <- mod$prob.logistic
			mainplot <- paste0( "Item Response Functions \nLogistic Model (Rasch Model)")
			plot.all.icc.isop.dich( pers , p1 , xlabplot , ylabplot , mainplot,I)
		graphics::par( mfrow=c(1,1))
		graphics::par(ask=ask)    
		#*********************
		# separate item response functions
		for (ii in 1:I){
		# ii <- 1
		graphics::par( mfrow=c(2,2))
			xlabplot <- expression( paste( "Modified P Score " , rho[p] ))
			ylabplot <- expression(paste(P,"(",X[i]==1, "|", rho[p] , ")" ))  
			mainplot <- paste0("Item " , colnames(mod3$dat)[ii] , " | Saturated Model" )
			graphics::plot( pers$mpsc , mod3$prob.saturated[ii,2, pers$scoregroup ] , type="l" ,
						 xlab=xlabplot , ylab=ylabplot , main=mainplot , ylim=c(0,1)  )
			mainplot <- paste0("Item " , colnames(mod3$dat)[ii] , " | ISOP Model" )               
			graphics::plot( pers$mpsc , mod3$prob.isop[ii,2, pers$scoregroup ] , type="l" ,
						 xlab=xlabplot , ylab=ylabplot , main=mainplot , ylim=c(0,1)  )
			mainplot <- paste0("Item " , colnames(mod3$dat)[ii] , " | ADISOP Model" )                              
			graphics::plot( pers$mpsc , mod3$prob.adisop[ii,2, pers$scoregroup ] , type="l" ,
						 xlab=xlabplot , ylab=ylabplot , main=mainplot , ylim=c(0,1)  )
			mainplot <- paste0("Item " , colnames(mod3$dat)[ii] , " | Logistic Model" )
			graphics::plot( pers$mpsc , mod3$prob.logistic[ii,2, pers$scoregroup ] , type="l",
						 xlab=xlabplot , ylab=ylabplot , main=mainplot , ylim=c(0,1)  )
		   graphics::par( mfrow=c(1,1))    
		   graphics::par(ask=ask)    
			}		
	}
#########################################################################	
	
#########################################################################
# aux plot function all item response functions
plot.all.icc.isop.dich <- function( pers , p1 , 
	xlabplot , ylabplot , mainplot , I ){
	# begin plot
	graphics::plot( pers$mpsc , p1[1,2, pers$scoregroup ] , type="l" , ylim=c(0,1),
		 xlab= xlabplot , ylab=ylabplot , main = mainplot)    
	for (ii in 2:I){
		graphics::lines( pers$mpsc , p1[ii,2, pers$scoregroup ] , lty=ii , col=ii )
					}
				}
#########################################################################
# aux plot function all item response functions
plot.all.icc.isop.poly <- function( pers , p1 , 
	xlabplot , ylabplot , mainplot , I , kk){
	# begin plot
	graphics::plot( pers$mpsc , p1[1,kk, pers$scoregroup ] , type="l" , ylim=c(0,1),
		 xlab= xlabplot , ylab=ylabplot , main = mainplot)    
	for (ii in 2:I){
		graphics::lines( pers$mpsc , p1[ii,kk, pers$scoregroup ] , lty=ii , col=ii )
					}
				}