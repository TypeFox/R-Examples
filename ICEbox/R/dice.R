dice = function(ice_obj, DerivEstimator){

	#error checking:
	if(class(ice_obj) != "ice"){
		stop("ice_obj is not a valid ice object.")
	}

	gridpts = ice_obj$gridpts

	if (missing(DerivEstimator)){
		EstimatorWrapper = function(y){
			D1tr( x = gridpts, y = supsmu(x=gridpts,y=y)$y)  #numerical derivative of supersmooth.
		}
	} else{
		#argument checking???
		EstimatorWrapper = function(y){
			DerivEstimator(y=y,x=gridpts)		
		}
	}

	#compute derivatives
	dice_obj = ice_obj
	dice_obj$d_ice_curves = t(apply(ice_obj$ice_curves, 1, FUN = EstimatorWrapper))
	dice_obj$ice_curves = NULL

	#do it for the pdp as well.
	dice_obj$dpdp = EstimatorWrapper(ice_obj$pdp)

	#figure out point on each curve that corresponds to observed X.
    col_idx_of_actual = c(1, 1 + cumsum(diff(ice_obj$xj)>0))
	row_idx = 1:nrow(dice_obj$d_ice_curves)
	actual_deriv_idx = cbind(row_idx, col_idx_of_actual)
	dice_obj$actual_deriv = dice_obj$d_ice_curves[actual_deriv_idx]

	#compute the sd of the derivatives at each gridpt.
	dice_obj$sd_deriv = apply(dice_obj$d_ice_curves, 2, sd)
	
	#clean up, make it of class 'dice'	
	dice_obj$actual_prediction = NULL
	dice_obj$predictfcn = NULL
	dice_obj$frac_to_build = NULL
	dice_obj$indices_to_build = NULL
	class(dice_obj) = "dice"
	
	invisible(dice_obj)
}
