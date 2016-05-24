transform_params <-
function(pars){
	## Transforms parameter values (from log scale to physical scale)
	
	## Define the order in which the parameters should be passed to the solver
	params <- c(mrna6_degradation_rate =1, mrna7_degradation_rate =1, mrna8_degradation_rate =1, p_degradation_rate = 1,  r6_Kd = 1, r6_h = 1, r11_Kd=1, r11_h = 1, r12_Kd = 1, r12_h = 1, pro6_strength = 1,  pro7_strength = 1, pro9_strength = 1, rbs6_strength = 1, rbs7_strength = 1, rbs8_strength = 1)
	
	## Add mrna degradation rates which are all ones
	temp <- rep(1, 3)
	names(temp) <- c("mrna6_degradation_rate",
	                 "mrna7_degradation_rate",
	                 "mrna8_degradation_rate"
	                 )
	pars <- c(temp,pars)
	
	## Rescale dissociation constants (concentrations) on log scale
	temp_select <- names(pars) %in% c("r6_Kd", "r11_Kd")
	pars[temp_select] <- 10^( pars[temp_select]  / 100 * 6 - 3)
	
	
	## Rescale strength on log scale
	temp_select <- names(pars) %in% c("pro6_strength", "pro7_strength", "pro9_strength", "rbs6_strength", "rbs7_strength", "rbs8_strength")
	pars[temp_select] <- 10^( pars[temp_select]  / 100 * 6 - 3)
	
	## Add known coefficients
	temp2 <- c(4,2,2,0.2)
	names(temp2) <- c("r6_h", "r11_h", "r12_h","r12_Kd")
	pars <- c(temp2,pars)
	
	## Rescale protein degradation rate on log scale
	pars[names(pars) == "p_degradation_rate"] <- 10^(pars[names(pars) == "p_degradation_rate"] / 100 * 4 - 3) 
	
	## Return parameters in correct order
	return(pars[names(params)])
}
