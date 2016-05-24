###################################################################################################
#' Algorithm Interface to ES + GLG
#'
#' Interface for mixed model tuning of an Evolution Strategy ES optimizing Gaussian Landscapes.
#' For each repeat of an ES parameter setting, a different seed for the ES is used. At the same time,
#' this seed is used to create different Gaussian Landscapes with the Gaussian Landscape Generator GLG, see \code{\link{spotGlgCreate}}.
#'
#' @param spotConfig Contains the list of spot configurations, results of the algorithm can be passed to this list instead of the .res file.
#'		  spotConfig defaults to "NA", and will only be passed to the Algorithm if spotConfig$spot.fileMode=FALSE. See also: \code{\link{spotGetOptions}}
#'			Items used are: \cr \cr
#'			alg.currentDesign: data frame holding the design points that will be evaluated \cr
#'			io.apdFileName: name of the apd file \cr
#'			io.desFileName: name of the des file \cr
#'			io.resFileName: name of the res file, for logging results (if spotConfig$spot.fileMode==TRUE)\cr
#'			spot.fileMode: boolean, if selected with true the results will also be written to the res file, otherwise it will only be saved in the spotConfig returned by this function\cr
#' @return this function returns the \code{spotConfig} list with the results in spotConfig$alg.currentResult
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}} \code{\link{optim}}
#' \code{\link{spotFuncStartBranin}} \code{\link{spotAlgStartSannVar}}
#' @export
###################################################################################################
spotAlgStartEsGlg <- function(spotConfig){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	io.apdFileName=spotConfig$io.apdFileName;
	io.desFileName=spotConfig$io.desFileName;
	io.resFileName=spotConfig$io.resFileName;	
	#default Values that can be changed with apd file
	# ES defaults
	control=list()
	control$maxit <- 100
	control$sigmaInit=1.0;
	control$nSigma=1;
	control$tau0=0.0; #
	control$tau=1.0;
	control$stratReco=3; #inter reco     or 1
	control$objReco=2; #inter reco       or 2
	control$kappa=-1; #1
	control$mue = 5;
	control$nu = 2;
	control$sigmaRestart=0;
	control$prescanmult = 1; 	
	control$mutation<-2;
	control$rho<-"bi";
	control$maxGen<-Inf;
	#control$glgSeed<-0 #This is not an ES setting, moved to GLG defaults below
	#GLG defaults
	dim=2
	lb <- rep(-1,dim)
	ub <- rep(1,dim)
	ngauss= 10
	maxval = 1
	ratio = 0.8 
	npinst = 9 #number of random instances
	glgSeed = 0 #starting seed for random problem instances
	## read problem design file
	if(file.exists(io.apdFileName)){
		source(io.apdFileName,local=TRUE)
	}
	else{
		spotWriteLines(spotConfig$io.verbosity,1,"apd File not found, defaults used");
	}		
	if (spotConfig$spot.fileMode){ 
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  io.desFileName), con=stderr());
		## read doe/dace etc settings:
		des <- read.table( io.desFileName, sep=" ", header = TRUE);	
	}else{
		des <- spotConfig$alg.currentDesign; 
	}
	pNames <- names(des);	
	config<-nrow(des);	
	for (j in 1:npinst){
		pinst <- glgSeed + j	#SEED passed to GLG
		fn <- spotGlgCreate(dimension=dim,nGaussian=ngauss,lower=lb, upper=ub, globalvalue=maxval,ratio=ratio,seed=pinst)
		fun <- function(x) maxval -fn(x) #use difference to opimum as target value
		for (k in 1:config){
			for (i in 1:des$REPEATS[k]){
				## PARAMETER VALUES OF THE ES, WHICH ARE TUNED
				if (is.element("NPARENTS", pNames)){
					control$mue <- des$NPARENTS[k]
				}
				if (is.element("NU", pNames)){
					control$nu <- des$NU[k]
				}
				if (is.element("NSIGMA", pNames)){
					control$nSigma <- des$NSIGMA[k]
				}
				if (is.element("TAU0", pNames)){
					control$tau0 <- des$TAU0[k]
				}
				if (is.element("TAU", pNames)){
					control$tau <- des$TAU[k]
				}
				if (is.element("KAPPA", pNames)){
					control$kappa <- des$KAPPA[k]
				}
				if (is.element("SIGMARESTART", pNames)){
					control$sigmaRestart <- des$SIGMARESTART[k]
				}
				if (is.element("SIGMAINIT", pNames)){
					control$sigmaInit <- des$SIGMAINIT[k]
				}
				if (is.element("PRESCANMULT", pNames)){
					control$prescanmult <- des$PRESCANMULT[k]
				}				
				if (is.element("OBJRECO", pNames)){
					control$objReco <- des$OBJRECO[k]
				}
				if (is.element("STRATRECO", pNames)){
					control$stratReco <- des$STRATRECO[k]
				}	  
				##	OUTPUT PARAMETERS
				conf <- k
				if (is.element("CONFIG", pNames)){
					conf <- des$CONFIG[k]
				}
				spotStep<-NA
				if (is.element("STEP", pNames)){
					spotStep <- des$STEP[k]
				}			
				seed <- des$SEED[k]+i-1	#SEED passed to ES, as well as GLG
				control$seed=seed
				optimres <- spotOptimEs(par= rep(NA,dim), fn = fun, lower= lb, upper= ub, control=control)	
				#### WRITE RESULTS
				res <- list(Y=optimres$value, # last value     #TODO: should be alltime best?
					NPARENTS=control$mue,
					NU=control$nu,
					KAPPA=control$kappa,
					InitSgm=control$sigmaInit,
					#TauMult=tauMult,###########################TODO?
					Rho=control$rho,
					NSIGMA=control$nSigma,
					SIGMARESTART=control$sigmaRestart,
					Mutation=control$mutation,
					TAU0=control$tau0,
					TAU=control$tau,
					SIGMAINIT=control$sigmaInit,
					PRESCANMULT=control$prescanmult,
					OBJRECO = control$objReco,
					STRATRECO = control$stratReco,
					#Function=fName,   #can lead to read/write errors, since it is a function and not a string now
					MaxIter=control$maxit,
					Dim=dim,
					PINST=pinst,          
					SEED=seed,
					CONFIG=conf
				)
				if (is.element("STEP", pNames)){
					res=c(res,STEP=spotStep)
				} 
				res <-data.frame(res)
				if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
					colNames = TRUE
					if (file.exists(io.resFileName)){
						colNames = FALSE
					}				
					write.table(res, file = io.resFileName, row.names = FALSE, 
						col.names = colNames, sep = " ", append = !colNames, quote = FALSE)
					colNames = FALSE					
				}
				spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res)#always log the results in spotConfig			
			}
		}	
	}	
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	spotConfig
}

