impacts.splm<-function(obj, listw = NULL, time = NULL, ..., tr=NULL, R=200, type="mult", empirical=FALSE, Q=NULL){

if(is.null(listw) && is.null(tr)) stop("either listw or tr should be provided")

	
if(!is.null(listw) ){	
	if(listw$style != "W") stop("Only row-standardised weights supported")
	if(is.null(time) && is.null(tr)) stop("time periods should be provided")
}


if(is.null(tr)){
	
sparse.W <- listw2dgCMatrix(listw)
s.lws <- kronecker(Diagonal(time) , sparse.W)
tr <- trW(s.lws, type= type)
	
	}
	
if(is.na(match(obj$type, c("fixed effects lag","fixed effects sarar","random effects ML", "fixed effects GM","lag GM","fixed effects GM")))) stop("object type not recognized")
	
	if(obj$type == "fixed effects lag"){
		
class(obj)<- "gmsar"	
obj$type <- "SARAR"
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2
obj$secstep_var <- obj$vcov
imp <- impacts(obj, tr=tr, R=R, ...)

	}
	
	if(obj$type == "fixed effects sarar"){

class(obj)<- "gmsar"	
obj$type <- "SARAR"
rho <- obj$coefficients[2]
obj$coefficients <- obj$coefficients[-2]
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2
obj$secstep_var <- obj$vcov[-2,-2]
imp <- impacts(obj, tr=tr, R=R,...)		
		
	}

	if(obj$type == "fixed effects error") stop("Impacts Estimates are not available for Error Model")

	if(obj$type == "random effects ML")	{

if(!is.null(obj$arcoef)) {
class(obj)<- "gmsar"	
obj$type <- "SARAR"

obj$coefficients <- c(obj$arcoef, obj$coefficients)
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2
obj$secstep_var <- matrix(0,nrow(obj$vcov)+1,nrow(obj$vcov)+1)
obj$secstep_var[1,1] <- obj$vcov.arcoef
obj$secstep_var[(2:(nrow(obj$vcov)+1)),(2:(nrow(obj$vcov)+1))] <- obj$vcov
imp <- impacts(obj, tr=tr, R=R, ...)		
		}
		else stop("Impacts Estimates are not available for Error Model")		
		
	}
	

	if(obj$type == "fixed effects GM"){
		
		if(is.null(obj$endog)) {
obj$secstep_var <- vcov(obj)			
class(obj)<- "gmsar"	
obj$type <- "SARAR"
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2

imp <- impacts(obj, tr=tr, R=R, ...)		
			
			
		}
				
		else stop("No impacts estimates when endogenous variables are present in the system")
					
	}

if(obj$type == "lag GM")			{
	
		if(is.null(obj$endog)) {

class(obj)<- "gmsar"	
obj$type <- "SARAR"
obj$secstep_var <- obj$var			
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2

imp <- impacts(obj, tr=tr, R=R, ...)		
			
			
		}
				
		else stop("No impacts estimates when endogenous variables are present in the system")
					

	
}


if(obj$type == "random effects GM")			{
	
		if(is.null(obj$endog)) {

class(obj)<- "gmsar"	
obj$type <- "SARAR"
obj$secstep_var <- obj$vcov			
obj$data <- as.vector(obj$model)
obj$s2 <- obj$sigma2

imp <- impacts(obj, tr=tr, R=R, ...)		
			
			
		}
				
		else stop("No impacts estimates when endogenous variables are present in the system")
					

	
}
		



	
return(imp)	
	
}