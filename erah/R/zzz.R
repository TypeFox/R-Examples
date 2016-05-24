
.onAttach <- function(libname, pkgname) {
	        
	erahV <- utils::packageVersion("erah")     
   
	msg <- paste("Welcome to eRah. This is a early release of eRah (V",erahV,"). For bugs, problems and issues, please use the eRah forum at: http://erah.lefora.com/. Describe your problem and include the output of sessionInfo().", sep="")
    
    packageStartupMessage(msg)            
	
}    
