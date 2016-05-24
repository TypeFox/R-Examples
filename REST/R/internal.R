.onAttach <- function(libname, pkgname){
	if (!interactive()){
		return()
	}
	else{
		Commander()
	}
	
}

if(getRversion() >= "2.15.1") globalVariables(c("top","button.command"))