#'@useDynLib YplantQMC
.onAttach <- function(lib, pkg) {
    
	  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- sprintf("\n-- Plant modelling with Yplant - QuasiMC for R (version %s)\n\tRefer to the manual for detailed instructions, \nor see the example in ?YplantDay.\n", 
		as.character(ver))
    packageStartupMessage(msg)
	
    # MC 4/12/2012 - updated to include Mac OS X	
  	if(.Platform$OS.type != "windows" && (Sys.info()[['sysname']] != "Darwin"))
  		packageStartupMessage("-- You are not using a Windows or Mac machine - YplantDay() and runYplant() will not work.")
	
	if(!checkInstallation()){
		packageStartupMessage("!- To Install QuasiMC, type installQuasiMC().")
                if (.Platform$OS.type == "windows"){
		        packageStartupMessage("   You need an internet connection - two files will be copied to c:/QuasiMC")
                }
                else if (Sys.info()[['sysname']] != "Darwin") {
                        qmc_dir = paste(path.expand("~"),"/QuasiMC/",sep="")
                        packageStartupMessage(paste("   You need an internet connection - two files will be copied to", qmc_dir,"\n"))
                }
	}
  
}

# .onUnload <- function(libpath){
#     library.dynam.unload("YplantQMC", libpath)
# 	#unloadNamespace("YplantQMC")
# }
