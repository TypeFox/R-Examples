# Initializes Java

.onLoad <- function(libname, pkgname) {
	if ((!is.null(Sys.info())) && (Sys.info()[1]=="Darwin")) {
#		options(java.parameters=c("-Xms1G","-Xmx1G","-Djava.awt.headless=true"))
		options(java.parameters=c("-Djava.awt.headless=true"))
	} else {
#		options(java.parameters=c("-Xms1G","-Xmx1G"))
	}
	.jpackage(pkgname);
	.jcall("freak/rinterface/control/LogRegInterface", "V", "setRMode");   ## Freak R init		
	.jcall("freak/rinterface/control/RFlags", "V", "setUseCase",.jfield("freak/rinterface/control/RFlags","I","R"));   	
}

#Example:
#none needed
