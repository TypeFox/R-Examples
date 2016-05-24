


###############################################################################
###############################################################################
###############################################################################
mice.1chain <- function(data, burnin = 10 , iter = 20 , Nimp = 10 , 
		method = vector("character", length = ncol(data)),
		predictorMatrix = (1 - diag(1, ncol(data))),
		visitSequence = (1:ncol(data))[apply(is.na(data), 2, any)],
		form = vector("character", length = ncol(data)),
		post = vector("character", length = ncol(data)),
		defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
		diagnostics = TRUE, printFlag = TRUE,
		seed = NA, imputationMethod = NULL,
		defaultImputationMethod = NULL, data.init = NULL, ...){
    ###################################################################------
	CALL <- match.call()
    # burnin phase
    iterstep0 <- round( seq( burnin , iter , length = Nimp+1 ) )
    iterstep <- diff( iterstep0  )
	implist <- list(  1:(Nimp+1) )
	datlist <- list( 1:Nimp )
#     data <- as.matrix(data )
	# indicator for factors
	VV <- ncol(data)
	ind.factor <- sapply( 1:VV , FUN = function(vv){ is.factor( data[,vv]) } )
	ind.num <- setdiff( 1:VV , which( ind.factor ) )	
	meansObs <- colMeans( data[,ind.num] , na.rm=TRUE )
	NObs <- colSums( 1 - is.na(data[,ind.num]))	
	VarObs <- ( colSums( data[,ind.num]^2 , na.rm=TRUE ) - NObs * meansObs^2 ) / ( NObs  ) 
	NMiss <- colSums( is.na(data[,ind.num]))		

	#******************************
    cat("************ BURNIN PHASE |", "Iterations 1 -" , burnin , 
			"******************\n") 
			utils::flush.console()
    implist[[1]] <- imp0 <- mice::mice( data , maxit=burnin , m=1 , 
                imputationMethod = imputationMethod ,
                predictorMatrix = predictorMatrix , method=method , visitSequence=visitSequence ,
                form=form , post=post , defaultMethod=defaultMethod ,  diagnostics=diagnostics ,
                printFlag=printFlag , seed=seed , defaultImputationMethod=defaultImputationMethod ,
                data.init = data.init , ... )
    dat0 <- mice::complete( imp0 , 1 )
	chainMean <- t( imp0$chainMean[,,1] )
	chainVar <- t( imp0$chainVar[,,1] )
	rownames(chainVar)[ nrow(chainVar) ] <- rownames(chainMean)[ nrow(chainMean) ] <- "end_burnin"
    # imputation phase
	
    cat("\n\n************ IMPUTATION PHASE |", "Iterations", burnin+1 , "-" , iter , 
			"******************\n") ; 
	utils::flush.console()
    for ( mm in 1:Nimp){
	    cat("\n --- Imputation" , mm , "| Iterations" ,
			iterstep0[mm] , "-" , iterstep0[mm+1] , "\n" ) ; 
		utils::flush.console()
		implist[[mm+1]] <- imp1 <- mice::mice( data , maxit=iterstep[mm] , m=1 , 
					imputationMethod = imputationMethod ,
					predictorMatrix = predictorMatrix , method=method , visitSequence=visitSequence ,
					form=form , post=post , defaultMethod=defaultMethod ,  diagnostics=diagnostics ,
					printFlag=printFlag , seed=seed , defaultImputationMethod=defaultImputationMethod ,
					data.init = dat0 , ... ) 		
        datlist[[mm]] <- dat0 <- mice::complete( imp1 , 1 )
		chainMean <- rbind( chainMean , t( imp1$chainMean[,,1] ) )
		chainVar <- rbind( chainVar , t( imp1$chainVar[,,1] ) )
		rownames(chainVar)[ nrow(chainVar) ] <- 
			rownames(chainMean)[ nrow(chainMean) ] <- paste0( "imp_" , mm )	
						}
        
    ###################################################--------------------------
    # post-processing

	# conversion into a mids object
	ids <- seq( 1 , nrow(data) )
	datalong <- data.frame( ".imp"=0 , "id"= ids , data )
    for ( mm in 1:Nimp ){ 
		datalong <- rbind( datalong , 
				data.frame( ".imp"=mm , "id" = ids , datlist[[mm]] ) )
						}
	# midsobj <- datalist2mids(dat.list=datlist, progress = FALSE)
	# as.mids in mice	
	midsobj <- mice::as.mids(data = datalong , .imp=1, .id=2)
	
	# adjust M and Var for chains
	vars <- colnames(chainMean)
	cM <- nrow(chainMean)
	eps <- 10^(-20)
	meansObs[ is.na(meansObs) ] <- 0
	chainMPar <- ( chainMean * matrix( NMiss[vars] , nrow=cM , ncol=length(vars) , byrow=TRUE ) +
				matrix( NObs[vars] * meansObs[vars] , nrow=cM , ncol=length(vars) , byrow=TRUE ) ) /
						( eps+  matrix( NObs[vars] + NMiss[vars] , nrow=cM , ncol=length(vars) , byrow=TRUE ) )
	vars <- colnames(chainVar)	
	VarObs[ is.na(VarObs) ] <- 0	
	chainVarPar <- ( chainVar * matrix( (NMiss[vars]) , nrow=cM , ncol=length(vars) , byrow=TRUE ) +
				matrix( ( NObs[vars])* VarObs[vars] , nrow=cM , ncol=length(vars) , byrow=TRUE ) ) /
						( eps+ matrix( NObs[vars] + NMiss[vars] , nrow=cM , ncol=length(vars) , byrow=TRUE ) )						
	midsobj <- mice::as.mids(datalong , .imp=1)
#	midsobj$m <- Nimp
	
	imm <- implist[[mm+1]]
	midsobj$predictorMatrix <- imm$predictorMatrix
	midsobj$method <- imm$method
#	midsobj$m <- Nimp
	
	midsobj$call <- CALL
	
	res <- list( "midsobj"=midsobj , "datlist"=datlist , "datalong" = datalong , 
					"implist" = implist ,
					 "chainMPar"=chainMPar , 
					"chainVarPar" =  chainVarPar   ,
					"chainMean" = chainMean , "chainVar" = chainVar ,
					"burnin" = burnin , "iter"=iter , "Nimp"=Nimp ,
					"time" = Sys.time()
					)
	res$nobs <- nrow(data)
	res$nvars <- ncol(data)					
	class(res) <- "mids.1chain"
    return(res)
            }
###############################################################################
###############################################################################
###############################################################################
# S3 method summary
summary.mids.1chain <- function( object , ... ){
   cat("Multiply imputed dataset using one chain\n\n")
#   cat("***" , paste0("Date : " , Sys.time() ) , "\n\n")
   cat( paste0( "Number of iterations = " , object$iter) , "\n")
   cat( paste0( "Number of burnin-iterations = " , object$burnin) , "\n")   
   cat( paste0( "Number of imputations = " , object$Nimp) , "\n")   
   cat( paste0( object$nobs , " cases and " , object$nvars , " variables \n\n" ) )   
   cat("---------\n")   
   print( summary( object$midsobj ) )
	}
# print method	
print.mids.1chain <- function( x , ...){
	summary.mids.1chain( object = x , ...)
	}	
	
################################################################################
# S3 plot method
plot.mids.1chain <- function( x , plot.burnin=FALSE , ask=TRUE , ... ){
    VV <- ncol( x$chainMean )
    graphics::par(mfrow=c(2,2))
	if (!plot.burnin){
		iter_vec <- seq( x$burnin + 1 , x$iter )
				} else { iter_vec <- 1:x$iter }
    
    for (vv in 1:VV){
        plot( iter_vec , x$chainMean[ iter_vec , vv ] , type="l" , 
                    xlab="Iterations" , ylab="M" , 
                    main=paste0("Mean " , colnames(x$chainMean)[vv] )
                                )
        plot( iter_vec , sqrt( x$chainVar[ iter_vec , vv ] ) , type="l" , 
                    xlab="Iterations" , ylab="SD" , 
                    main=paste0("SD " , colnames(x$chainMean)[vv] )
                                )
        if (vv %% 2 == 0){ par(ask=ask) }
                }
    graphics::par(mfrow=c(1,1))
        }
