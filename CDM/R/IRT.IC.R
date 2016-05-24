
##############################################################################
# information criteria
IRT.IC <- function( object ){
    ll <- logLik(object)
	res <- c( ll , -2*ll , attr(ll, "df") , attr(ll,"nobs" ) )
	names(res) <- c("loglike" , "Deviance" , "Npars" , "Nobs" )
	p <- Npars <- res["Npars"]
	n <- res["Nobs"]
	res["AIC"]  <- -2*ll + 2*Npars
	res["BIC"]  <- -2*ll + log(n)*p 
	res["AIC3"] <- -2*ll + 3*Npars
    res["AICc"] <- -2*ll + 2*p + 2*p*(p+1)/(n-p-1)
	res["CAIC"] <- -2*ll + (log(n)+1)*p
    return(res)
                }
##############################################################################				
#			"   | AIC = -2*LL + 2*p  \n" )    
#		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
#			"   | BIC = -2*LL + log(n)*p  \n" )  
#		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   
##############################################################################
				