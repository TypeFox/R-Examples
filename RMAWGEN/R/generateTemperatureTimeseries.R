
NULL


#'  
#'  Returns time series of Daily Maximum and Minimum with a random multi-realization obtained by using \code{\link{newVARmultieventRealization}}. This function is called by  \code{\link{ComprehensiveTemperatureGenerator}}.
#'  
#'  @param std_tn vector containing standard deviation of daily minimum temperature anomalies. \code{stdTn} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#' 	@param std_tx vector containing standard deviation of daily maximum temperature anomalies. \code{stdTx} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#' 	@param SplineTx matrix containing the averaged daily maximum temperature  obtained by a spline interpolation of monthly means . \code{SplineAdvTx} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#'  @param SplineTn matrix containing the averaged daily minimum temperature  obtained by a spline interpolation of monthly means . \code{SplineAdvTn} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#'  @param SplineTm matrix containing the averaged daily "mean" temperature   obtained by a spline interpolation of monthly means . \code{SplineAdvTm} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#'  @param SplineDeltaT matrix containing the rescaled averaged daily temperature range obtained by a spline interpolation of monthly means . 
#'         \code{SplineAdvDelta_T_sim/SplineAdvDelta_T} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#'  @param std_tm vector containing standard deviation of daily "mean" temperature anomalies. \code{stdTn} is default, see \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#'  @param var A VAR model represented by a \code{varest} object as returned by \code{\link{getVARmodel}} or \code{\link{VAR}}
#'  @param exogen see \code{\link{VAR}}
#'  @param normalize logical variable If \code{TRUE} \code{\link{normalizeGaussian_severalstations}} is used, otherwise not. If \code{option} is 2, it is always \code{TRUE}.
#'  @param sample,origin_x,origin_data,extremes see \code{\link{normalizeGaussian_severalstations}}
#'  @param type see \code{\link{quantile}}
#'  @param option integer value. If 1, the generator works with minimum and maximum temperature, if 2 (Default) it works with th average value between maximum and minimum temparature and the respective daily Thermal Range.
#'  @param original_data matrix containing the measured standardized temperature anomalies
#'  @param noise stochastic noise to add for variabile generation. Default is \code{NULL}. See \code{\link{newVARmultieventRealization}}. 
#' 
#'  @export 
#' @return  This function returns a list of the following variables: 
#' 
#' \code{res_multigen} matrix containing standardized values of daily maximum and minimum temperature anomalies
#' 
#' \code{Tx_spline} matrix containing climatic "spline-interpolated" daily maximum temperature 
#' 
#' \code{Tn_spine} matrix containing climatic "spline-interpolated" daily minimum temperature 
#' 
#' \code{Tx_gen} matrix containing generated daily maximum daily temperature (\eqn{Tx_{gen}})
#' 
#' \code{Tn_gen} matrix containing generated daily minimum daily temperature (\eqn{Tn_{gen}})
#' 
#' \code{Tm_gen} matrix containing generated "mean" daily temperature defined as  \eqn{\frac{Tx_{gen}+Tn_{gen}}{2}}
#' 
#' \code{DeltaT_gen} matrix containing generated daily thermal range defined as  \eqn{Tx_{gen}-Tn_{gen}}
#' 
#' See the R code for further details
#' 
#' @seealso \code{\link{newVARmultieventRealization}},\code{\link{normalizeGaussian_severalstations}}
#' 
#' @author Emanuele Cordano, Emanuele Eccel
#'        
#' 



generateTemperatureTimeseries <-
function (std_tn,std_tx,SplineTx,SplineTn,
		SplineTm,SplineDeltaT,std_tm,
		var=NULL,exogen=NULL,normalize=TRUE,type=3,extremes=TRUE,sample=NULL,option=1,original_data,origin_x=NULL,origin_data=NULL,noise=NULL) {
	
	#  @author  Emanuele Cordano
	#    
	#   
	#
	#        @note Calculated complete generated time series of Daily Maximum Temperature and Daily Minum Temperature with a random multi-realization obtained by using NewMultiRealizations function and saves the results as global variables 
	# @return 0 in case of success, -1 otherwise
	
	if (is.null(var)) { 
		
		print("Error in generateTemperatureTimeseries: VAR model (varest object) is missing!")
		return(-1)
		
	} else {
		
		res_multigen0 <- newVARmultieventRealization(var=var,exogen=exogen,nrealization=nrow(SplineTx),type=type,extremes=extremes,noise=noise)
		
		
		
		if (normalize) {
			
			
			res_multigen <- normalizeGaussian_severalstations(x=res_multigen0,data=original_data,inverse=TRUE,type=type,sample=sample,origin_x=origin_x,origin_data=origin_data,extremes=extremes)
			
			
		} else {
			
			res_multigen <- res_multigen0
			
		}
		
		
	}

	if (option==1) {
		
		Tx_gen <- extractTxFromAnomalies(res_multigen,std=std_tx,SplineAdv=SplineTx) # first-helf column
		Tn_gen <- extractTnFromAnomalies(res_multigen,std=std_tn,SplineAdv=SplineTn) # second-half column
		
		
	} else if (option==2) {
		
		
		Tm_gen <- extractTxFromAnomalies(res_multigen,std=std_tm,SplineAdv=SplineTm)
		
		ntall <- as.integer(ncol(res_multigen))
		ntn <- as.integer(ncol(res_multigen)/2)
		
		
		DeltaT_gen <- res_multigen[,(ntn+1):ntall] #*SplineDeltaT removed by EC 20100524
		
		
		Tx_gen <- Tm_gen+DeltaT_gen/2.0
		Tn_gen <- Tm_gen-DeltaT_gen/2.0
		
		
	} else if (option==3) {	
		# extract_average_temperature 
		Tm_gen <- res_multigen%*% diag(std_tm) + SplineTm
		
		Tx_gen <- Tm_gen-SplineTm+SplineTx
		Tn_gen <- Tm_gen-SplineTm+SplineTn
	}
	
	out <- list(res_multigen,SplineTx,SplineTn,Tx_gen,Tn_gen,original_data)
	
	names(out) <- c("res_multigen","Tx_spline","Tn_spline","Tx_gen","Tn_gen","data_original")
	return(out)
	
}

