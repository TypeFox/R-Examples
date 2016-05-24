
###########################################################################
IRT.jackknife.gdina <- function( object , repDesign , ... ){
	
	    rdes <- repDesign
		# read control arguments
		ctl <- object$control
		# replicate weights
		wgtrep <- repDesign$wgtrep
		RR <- ncol(wgtrep)
		# convert polychoric correlations
		pcvec <- CDM.polychorList2vec(polychorList=object$polychor)
		# assess modelfit
		fmod <- IRT.modelfit( object )$modelfit.stat
		fitvars <- c("MADcor","SRMSR" , "100*MADRESIDCOV")	
		# read parameter table
		partable <- object$partable
		
		#**** create parameter table	
		jpartable <- partable
		jpartable$skillclass <- jpartable$varyindex <- jpartable$fixed <-
				jpartable$free <- jpartable$totindex <- jpartable$rule <- NULL
		jpartable$group <- NULL
		#---- polychoric correlation
		dfr <- data.frame("partype" = "polychor" , "parindex" =NA , "item" = 0 ,
				"item.name" = "" , "parnames" = names(pcvec)  ,
				"value" = pcvec )
		rownames(dfr) <- NULL
		jpartable <- rbind( jpartable , dfr )
		#----- model fit
		dfr <- data.frame("partype" = "modelfit" , "parindex" =NA , "item" = 0 ,
				"item.name" = "" , "parnames" = fitvars  ,
				"value" = fmod[ fitvars , "est" ] )
		rownames(dfr) <- NULL
		jpartable <- rbind( jpartable , dfr )		
		
		
		#---------------------------------------
		# read control parameters for GDINA function
		args1 <- c( "skillclasses" ,  "q.matrix" , "conv.crit" ,
					"dev.crit" ,  "maxit"  ,   "linkfct" , "Mj"   ,
  				    "group","method","delta.designmatrix","delta.basispar.lower",
					"delta.basispar.upper","delta.basispar.init","zeroprob.skillclasses",
					"reduced.skillspace","HOGDINA","Z.skillspace",
					"weights","rule", "increment.factor",
					"fac.oldxsi","avoid.zeroprobs" , "delta.fixed" 
							)

		inputlist <- list( "data" = object$data ,  
						   "delta.init" = object$delta , "attr.prob.init" = ctl$attr.prob )
		A1 <- length(args1)
		for (vv in 1:A1){ 
				a.vv <- args1[vv]
				inputlist[[ a.vv ]] <- ctl[[ a.vv ]] 
					}
		inputlist$progress <- FALSE
		inputlist$calc.se <- FALSE
		
		#-----------------------------------
		# create jackknife parameter table
		NP <- nrow(jpartable)
		parsM <- matrix( NA , nrow=NP , ncol=RR)
		rownames(parsM) <- jpartable$parnames
		NP1 <- nrow( object$partable )
				
		#*****************************************
		# loop over datasets
		# cat("0  |" )
		for (rr in 1:RR){
		#	rr <- 1
			if ( rr %% 10 == 1 ){
				cat( paste0( "\n",rr-1 , " |" )) 
				utils::flush.console()
						}
			inputlist$weights <- rdes$wgtrep[,rr]
			mod1a <- do.call( gdina , inputlist )
			# convert polychoric correlations
			pcvec <- CDM.polychorList2vec(polychorList=mod1a$polychor)
			# assess modelfit
			fmod <- IRT.modelfit( mod1a)$modelfit.stat					
			parsM[ 1:NP1 , rr ] <- mod1a$partable$value			
			parsM[ jpartable$partype == "polychor"  , rr ] <- pcvec
			parsM[ jpartable$partype == "modelfit"  , rr ] <- fmod[ fitvars , "est" ] 		
			cat("-") ; utils::flush.console()
			if ( rr %% 10 == 0 ){
				cat( paste0( "|" )) 
				utils::flush.console()
						}			
						}
			cat("\n")
		fayfac <- rdes$fayfac
		res0 <- jkestimates( est=jpartable$value , parsM , fayfac  )			
		jpartable$jkest <- res0$dfr$jkest		
		jpartable$jkse <- res0$dfr$jkse		
		jpartable$t <- jpartable$jkest / jpartable$jkse
		jpartable$p <- 2*stats::pnorm( - abs( jpartable$t ) )
		
		res <- list( "parsM" = parsM , "vcov" = res0$vcov_pars , 
				      "jpartable" = jpartable , "fayfac" = fayfac )
		class(res) <- "IRT.jackknife"
		return(res)

				}
###########################################################################		
summary.IRT.jackknife <- function( object , digits = 3 ,  ... ){	
		obji <- object$jpartable
		V1 <- which( colnames(obji) == "value" )
		V2 <- ncol(obji)
		for (vv in V1:V2){
			obji[,vv] <- round( obji[,vv]  , digits )
						}
		print( obji )		
	
		}
################################################		
coef.IRT.jackknife <- function( object ,  bias.corr=FALSE , ... ){
		obji <- object$jpartable
		# bias.corr <- FALSE
		if ( bias.corr ){ 
				vari <- "jkest"
					} else {
				vari <- "value"
						}
		vec <- obji[ , vari ]
		names(vec) <- obji$parnames
		return(vec)
			}
##################################################			
vcov.IRT.jackknife <- function( object , ... ){
		return(object$vcov)
			}
##################################################			