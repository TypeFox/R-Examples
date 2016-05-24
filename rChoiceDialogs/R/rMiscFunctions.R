

#' Check if Tcl/Tk graphics can be used
#'
#' @name canUseTclTk
#' @return TRUE if Tcl/Tk graphics can be used and FALSE otherwise.
#' @export canUseTclTk
#' @author  Alex Lisovich, Roger Day

canUseTclTk<-function(){
	if(!interactive())
		return(FALSE);
	if(getOption("menu.graphics") && 
		.Platform$OS.type == "windows" &&
		 capabilities("tcltk")) 
			return(TRUE);
			
	if(getOption("menu.graphics") &&
		capabilities("tcltk") &&
		capabilities("X11")) 
			return(TRUE);
			
	return(FALSE);
}


#' Check if Java graphics can be used
#'
#' @return TRUE if Java graphics can be used and FALSE otherwise.
#' @export canUseJava
#' @author  Alex Lisovich, Roger Day

canUseJava<-function(){
	if(!interactive())
		return(FALSE);

	javaInstalled<-suppressWarnings(system("java -version",ignore.stderr = TRUE));
	if(javaInstalled!=0)
		return(FALSE);

	if(substring(version$os, 1, 6)=="darwin")
		return(FALSE);

	return(TRUE);
}


#' Check if modal Java dialogs can be used
#'
#' @name canUseJavaModal
#' @return TRUE if modal Java dialogs can be used and FALSE otherwise.
#' @export canUseJavaModal
#' @author  Alex Lisovich, Roger Day

canUseJavaModal<-function(){
	if(!interactive())
		return(FALSE);
	if(substring(version$os, 1, 6)=="darwin" && .Platform$GUI=="AQUA")
		return(FALSE);
	if(.Platform$OS.type == "windows" && .Platform$GUI[1]=="Rgui")
		return(FALSE);

	return(TRUE);
}


