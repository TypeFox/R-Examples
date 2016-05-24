	
	runquasimc <- function(cfgfile, inputfile=NA, outputfile=NA, debug=FALSE, intern=TRUE){
	     if(.Platform$OS.type == "windows"){
         
			if(debug){
				cmd <- paste("c:\\QuasiMC\\QuasiMC.exe -e c:\\QuasiMC\\enviro.e",cfgfile, 
                     "-debug <", inputfile, ">", outputfile)
			} else {
				cmd <- paste("c:\\QuasiMC\\QuasiMC.exe -e c:\\QuasiMC\\enviro.e", cfgfile, 
                     "-no_xserver <", inputfile, ">", outputfile)
			}
      
			if(intern)
				shell(cmd, intern=TRUE)
			else
				shell(cmd, intern=FALSE)
		
		} else if (Sys.info()[['sysname']] == "Darwin"){
                        qmc_app = paste(path.expand("~"),"/QuasiMC/QuasiMC.app/Contents/MacOS/QuasiMC",sep="")
                        qmc_env = paste(path.expand("~"),"/QuasiMC/enviro.e",sep="")
      
      if(debug)
				cmd <- paste(qmc_app, "-e", qmc_env, cfgfile, "-debug <", inputfile, ">", outputfile)
			else
				cmd <- paste(qmc_app, "-e", qmc_env, cfgfile, "-no_xserver <", inputfile, ">", outputfile)
             
			system(cmd)

		   } else {
			stop("QuasiMC is not currently available for your OS.")
		}
	}
