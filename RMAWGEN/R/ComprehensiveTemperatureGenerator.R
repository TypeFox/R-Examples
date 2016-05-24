NULL
#' 
#' The Comprehensive Temperature Generator
#' 
#' @param station see respective input parameter on \code{\link{setComprehensiveTemperatureGeneratorParameters}}
#' @param Tx_all,Tn_all,mean_climate_Tn,mean_climate_Tx,Tx_spline,Tn_spline see respective input parameter on \code{\link{setComprehensiveTemperatureGeneratorParameters}}
#' @param year_max,year_min,leap,nmonth,verbose see respective input parameter on \code{\link{setComprehensiveTemperatureGeneratorParameters}}
#' @param p,type,lag.max,ic,activateVARselect see respective input parameter on \code{\link{getVARmodel}}
#' @param year_max_sim last year of the simulation period. Default is equal to \code{year_max} 
#' @param year_min_sim first year of the simulation period. Default is equal to \code{year_min}
#' @param mean_climate_Tn_sim monthly averaged daily minimum temperatures for the simulated scenario and used by the random generator .  Default is \code{mean_climate_Tn}
#' @param mean_climate_Tx_sim monthly averaged daily maximum temperatures for the simulated scenario and used by the random generator .  Default is \code{mean_climate_Tx}
#' @param Tx_spline_sim daily timeseries (from the first day of \code{year_min_sim} to the last day of \code{year_max_sim}) of averaged maximum temperature which can be obtained by a spline interpolation of monthly mean values (for the generation period). Default is \code{Tx_spline}. See for spline interpolation utilized \code{\link{splineInterpolateMonthlytoDailyforSeveralYears}}.
#' @param Tn_spline_sim daily timeseries (from the first day of \code{year_min_sim} to the last day of \code{year_max_sim}) of averaged minimum temperature which can be obtained by a spline interpolation of monthly mean values (for the generation period). Default is \code{Tn_spline}. See for spline interpolation utilized \code{\link{splineInterpolateMonthlytoDailyforSeveralYears}}.

#' @param onlygeneration logical variable. If \code{TRUE} the VAR model \code{varmodel} is given as input and only random generation is done, otherwise (default) is calculated from measured data 

#' @param varmodel the comprehensinve VAR model as a \code{\link{varest2}} or \code{\link{GPCAvarest2}} S4 object or a \code{NULL} object. If \code{NULL} (default), the comprehensinve VAR is estimated from measured data within the function, otherwise it is given as input and only random generation is done.
# #### @param varmodel the VAR model as a \code{\link{varest2}}  or a \code{\link{GPCAvarest2}} object. If \code{NULL}, it is  given as input and only random generation is done, otherwise (default) is calculated from measured data 
#' @param normalize,sample,extremes see \code{\link{normalizeGaussian_severalstations}} or \code{\link{setComprehensiveTemperatureGeneratorParameters}}
#' @param type_quantile see \code{type} on \code{\link{quantile}}
#' @param option integer value. If 1, the generator works with minimun and maximum temperature, if 2 (default) it works with the average value between maximum and minimum temparature and the respective daily thermal range.
#' @param n_GPCA_iteration number of iterations of Gaussianization process for data. Default is 0 (no Gaussianization) 
#' @param n_GPCA_iteration_residuals number of iterations of Gaussianization process for VAR residuals. Default is 0 (no Gaussianization)
#' @param exogen data frame or matrix containing the (normalized or not) exogenous variables (predictors) for the recorded (calibration) period. Default is \code{NULL}.
#' @param exogen_sim  data frame or matrix containing the (normalized or not) exogenous variables (predictors) for the simulation period. Default is \code{NULL}. If it is \code{NULL}, \code{exogen_sim} is set equal to \code{exogen} within the function.
#' @param is_exogen_gaussian logical value, If \code{TRUE}, \code{exogen_sim} and \code{exogen} are given as already normalized variables, otherwhise they are not normalized. Default is \code{FALSE}
#' @param exogen_all data frame containing exogenous variable formatted like \code{Tx_all} and {Tn_all}. Default is \code{NULL}. 
#' It is alternative to \code{exogen} and if it not \code{NULL},\code{is_exogen_gaussian} is automatically set to \code{FALSE}	
#' @param exogen_all_col vector of considered  columns of \code{exogen_all}. Default is \code{station}.
#' @param nscenario number of generated scenarios for daily maximum and minimum temperature
#' @param yearly  logical value. If \code{TRUE} the monthly mean values are calculated for each year from \code{year_min} to \code{year_max} separately. Default is \code{FALSE}.  
#' @param yearly_sim logical value. If \code{TRUE} the monthly mean values are calculated for each year from \code{year_min_sim} to \code{year_max_sim} separately. Default is \code{yearly}. 
#' @param seed seed for stochastic random generation see \code{\link{set.seed}}
#' @param noise stochastic noise to add for variabile generation. Default is \code{NULL}. See \code{\link{newVARmultieventRealization}}. Not used in case that \code{nscenario>1}.
#' 
#'   
#' @export 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'    
#' @seealso \code{\link{setComprehensiveTemperatureGeneratorParameters}}, \code{\link{generateTemperatureTimeseries}} ,\code{\link{generateTemperatureTimeseries}},\code{\link{splineInterpolateMonthlytoDailyforSeveralYears}}. 
#' 
#'        
#'
#' @note It pre-processes series and generates multi-site temperature fields by using \code{\link{setComprehensiveTemperatureGeneratorParameters}},\code{\link{getVARmodel}} and \code{\link{generateTemperatureTimeseries}}. Detailed examples can be viewed of this function in \href{https://docs.google.com/file/d/0B8xDtMCnW3dJU2JIemVqMnpKTHc/edit}{this presentation}.
#' 
#'  
#' @return  A list of the following variables: 
#' 
#' \code{input}   list of variables returned by  \code{\link{setComprehensiveTemperatureGeneratorParameters}}
#' 
#' \code{var}     varest object containing the used VAR model (if useVAR is true), \code{NULL} (otherwise)
#' 
#' \code{output}  list variables returned by  \code{\link{generateTemperatureTimeseries}} (i.e. generated timeseries)
#' 
#' 
#' 
#' @examples 
#' 
#' data(trentino)
#' 
#' set.seed(1222) # set the seed for random generations!
#' year_min <- 1961
#' year_max <- 1990
#' 
#' year_min_sim <- 1982
#' year_max_sim <- 1983
#' 
#' n_GPCA_iter <- 5 
#' n_GPCA_iteration_residuals <- 5
#' p <- 1
#' vstation <- c("B2440","B6130","B8570","B9100","LAVIO","POLSA","SMICH","T0001",
#'  "T0010","T0014","T0018","T0032","T0064","T0083","T0090","T0092",
#' "T0094","T0099","T0102","T0110","T0129","T0139","T0147","T0149",
#' "T0152","T0157","T0168","T0179","T0189","T0193","T0204","T0210",
#' "T0211","T0327","T0367","T0373")		
#' ## Not Run: the call to ComprehensiveTemperatureGenerator may elapse 
#' ## too long time (more than 5 eseconds) and is not executed  by CRAN check.  
#' ## Please uncomment the following line to run the example on your own PC.
#' # generation00 <-ComprehensiveTemperatureGenerator(station=vstation[16],
#' # Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,
#' # p=p,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,
#' # sample="monthly",year_min_sim=year_min_sim,year_max_sim=year_max_sim)
#' 
#' 
#' 
#' 
#' 	

ComprehensiveTemperatureGenerator <-
function(
		station=c("T0001","T0010","T0099"),
		Tx_all,
		Tn_all,
		mean_climate_Tn=NULL,
		mean_climate_Tx=NULL,
		Tx_spline=NULL,
		Tn_spline=NULL,
		year_max=1990,
		year_min=1961,
		leap=TRUE,
		nmonth=12,
		verbose=TRUE,
		p=1,
		type="none",
		lag.max=NULL,
		ic="AIC",
		activateVARselect=FALSE,
		year_max_sim=year_max,
		year_min_sim=year_min,
		mean_climate_Tn_sim=NULL,
		mean_climate_Tx_sim=NULL,
		Tn_spline_sim=NULL,
		Tx_spline_sim=NULL,
		onlygeneration=FALSE,
		varmodel=NULL,normalize=TRUE,
		type_quantile=3,
		sample=NULL,
		extremes=TRUE,
		option=2,
		yearly=FALSE,
		yearly_sim=yearly,
		n_GPCA_iteration=0,
		n_GPCA_iteration_residuals=n_GPCA_iteration,
		exogen=NULL,
		exogen_sim=exogen,
		is_exogen_gaussian=FALSE,
		exogen_all=NULL,
		exogen_all_col=station,
		nscenario=1,
		seed=NULL,
		noise=NULL
	
) {
	

	
	useVAR=TRUE	
	if (option==2) normalize=TRUE # set NORMALIZE TRUE always for option=2
	
	origin <- paste(year_min,"1","1",sep="/") # Must start from Jan 1 
	origin_sim <- paste(year_min_sim,"1","1",sep="/") # Must start from Jan 1 
	
	
#	if (!onlygeneration){
	param <- setComprehensiveTemperatureGeneratorParameters(station=station,
				Tx_all=Tx_all,
				Tn_all=Tn_all,
				mean_climate_Tn=mean_climate_Tn,
				mean_climate_Tx=mean_climate_Tx,
				Tx_spline=Tx_spline,
				Tn_spline=Tn_spline,
				year_max=year_max,
				year_min=year_min,
				leap=leap,
				nmonth=nmonth,
				verbose=verbose,
				cpf=NULL,
				normalize=normalize,sample=sample,option=option,yearly=yearly)
		
		# PUT HERE CHECK FOR exogen 
	if (!onlygeneration){
		if (!is.null(exogen_all)) {
			
			 exogen <- as.data.frame(extractyears(exogen_all,year_min=year_min,year_max=year_max,station=exogen_all_col))
			 is_exogen_gaussian=FALSE
			 if (is.null(exogen_sim)) exogen_sim <- exogen
			
		}
		if (!is.null(exogen) & (!is_exogen_gaussian))  {
					
			exogen0 <- exogen
			exogen <- normalizeGaussian_severalstations(x=exogen0,data=exogen0,sample=sample,cpf=NULL,origin_x=origin,origin_data=origin,extremes=extremes)					
			
			
		}	
		var <- getVARmodel(data=param[['data_for_var']],suffix=c("_T1","_T2"),sep="",p=p,type=type,lag.max=lag.max,ic=ic,activateVARselect=activateVARselect,exogen=exogen,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,n_GPCA_iteration=n_GPCA_iteration,extremes=extremes) 
		
		if (activateVARselect) return(list(input=param,varselect=var))

	} else {

		var <- varmodel
	}
	
	if(is.null(Tx_spline)) Tx_spline <- param[['Tx_spline']]
	if(is.null(Tn_spline)) Tn_spline <- param[['Tn_spline']]
	

	
	
	nyear_sim <- year_max_sim-year_min_sim+1
	
	
	
# TO BE MODIFIED 	
	if (is.null(Tx_spline_sim)) {
		if (is.null(mean_climate_Tx_sim)) mean_climate_Tx_sim <- param[['monthly_mean_Tx']]
		Tx_spline_sim <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val=mean_climate_Tx_sim,start_year=year_min_sim,nyear=nyear_sim,leap=leap,yearly=yearly_sim))
#		str(Tx_spline_sim)
#		str(mean_climate_Tx_sim)
#		str(names(Tx_spline_sim))
#		str(names(mean_climate_Tx_sim))
#		names(Tx_spline_sim) <- names(mean_climate_Tx_sim)
		
		if (yearly_sim) {
			
			names(Tx_spline_sim) <- colnames(mean_climate_Tx_sim[[1]])
		} else {
			names(Tx_spline_sim) <- colnames(mean_climate_Tx_sim)
		}
		
	}
	
	if (is.null(Tn_spline_sim)) {
		
		if (is.null(mean_climate_Tn_sim)) mean_climate_Tn_sim <- param[['monthly_mean_Tn']]
		Tn_spline_sim <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val=mean_climate_Tn_sim,start_year=year_min_sim,nyear=nyear_sim,leap=leap,yearly=yearly_sim))
	#	names(Tn_spline_sim) <- names(mean_climate_Tn_sim)
		if (yearly_sim) {
		
			names(Tn_spline_sim) <- colnames(mean_climate_Tn_sim[[1]])
		} else {
			names(Tn_spline_sim) <- colnames(mean_climate_Tn_sim)
		}
	
	
	}

#	if (is.null(Tx_spline_sim)) Tx_spline_sim <- Tx_spline
#	if (is.null(Tn_spline_sim)) Tn_spline_sim <- Tn_spline	
	
	# THIS LINE IS TEMPORARY AND ONLY FOR TESTING
	if (!is.null(noise)) {if (noise=="residuals") noise <- residuals(var) }
	
	 
	
	SplineAdvTm_sim <- (Tx_spline_sim+Tn_spline_sim)/2.0
	SplineAdvDeltaT_sim <- (Tx_spline_sim-Tn_spline_sim)
	if (!is.null(exogen_sim) & (!is_exogen_gaussian))  {
		
		exogen0_sim <- exogen_sim
		exogen_sim <- normalizeGaussian_severalstations(x=exogen0_sim,data=exogen0_sim,sample=sample,cpf=NULL,origin_x=origin_sim,origin_data=origin_sim,extremes=extremes)					
	}
	if (!is.null(seed))	set.seed(seed)
	original_data <- param[['data_original']]
	
	if (option==2) { 
	#	print(class(original_data))
		
		
		
	
		ntall <- as.integer(ncol(original_data))
		ntn <- as.integer(ncol(original_data)/2)+1
#		str(Tx_spline)
#		str(Tx_spline_sim)
#		str(Tn_spline)
#		str(Tn_spline_sim)
	
		if (nrow(Tx_spline_sim)<nrow(original_data)) {
			
			nfrac <- as.integer(nrow(original_data)/nrow(Tx_spline_sim))+1
			Tx_spline_sim2 <- mapply(rep,Tx_spline_sim,nfrac)
			Tn_spline_sim2 <- mapply(rep,Tn_spline_sim,nfrac)
		} else {
			
			Tx_spline_sim2 <- Tx_spline_sim
			Tn_spline_sim2 <- Tn_spline_sim
		}
	
		ntn_rows <- 1:nrow(original_data)
#		str(ntn_rows)
	
	
	
	
	
		original_data[,ntn:ntall] <- original_data[,ntn:ntall]/(Tx_spline[ntn_rows,station]-Tn_spline[ntn_rows,station])*(Tx_spline_sim2[ntn_rows,station]-Tn_spline_sim2[ntn_rows,station])
	
	}
	

	
	results <- generateTemperatureTimeseries(std_tn=param[['stdTn']],std_tx=param[['stdTx']],SplineTx=Tx_spline_sim,SplineTn=Tn_spline_sim,SplineTm=SplineAdvTm_sim,SplineDeltaT=SplineAdvDeltaT_sim,std_tm=param[['stdTm']],var=var,normalize=normalize,type=type_quantile,sample=sample,option=option,original_data=original_data,origin_x=origin_sim,origin_data=origin,exogen=exogen_sim,extremes=extremes,noise=noise)	
	
	if (nscenario>1) {
		
		# TO BE PARALLELIZED SOMEHOW
		
		for (kk in 2:nscenario) {
			results_temp <- generateTemperatureTimeseries(std_tn=param[['stdTn']],std_tx=param[['stdTx']],SplineTx=Tx_spline_sim,SplineTn=Tn_spline_sim,SplineTm=SplineAdvTm_sim,SplineDeltaT=SplineAdvDeltaT_sim,std_tm=param[['stdTm']],var=var,normalize=normalize,type=type_quantile,sample=sample,option=option,original_data=param[['data_original']],origin_x=origin_sim,origin_data=origin,exogen=exogen_sim,extremes=extremes)	
			Tx_index <- sprintf("Tx_gen%05d",kk)
			Tn_index <- sprintf("Tn_gen%05d",kk)
			results[[Tx_index]] <- results_temp$Tx_gen
			results[[Tn_index]] <- results_temp$Tn_gen
			
		}
	} 

	if (onlygeneration) {
		
		return(list(output=results))
		
	} else {
		
		return(list(input=param,var=var,output=results,temporary=original_data))
		
	}
	return(list(input=param,var=var,output=results))	
	
	
}

