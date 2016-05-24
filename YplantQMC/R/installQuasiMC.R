#'Download and install the QuasiMC executable
#'
#'@description Downloads the Windows/Mac executable necessary for running
#'\code{\link{runYplant}} and \code{\link{YplantDay}}.
#'
#'QuasiMC is developed by Cieslak et al. (2008), which is part of the Virtual
#'Laboratory plant modelling software created at the University of Calgary
#'(\url{http://www.algorithmicbotany.org}).
#'
#'Both files are placed in the directory (for windows, \code{c:/QuasiMC}),
#'which is created if it does not exist already.
#'
#'This function can also be used to update an existing installation.
#'
#'
#'@param proxy Use TRUE if connecting through a proxy. This invokes
#'\code{setInternet2} - only works on Windows.
#'@author QuasiMC by Mik Cieslak, R installation code by Remko Duursma.
#'@references Cieslak, M., C. Lemieux, J. Hanan and P. Prusinkiewicz. 2008.
#'Quasi-Monte Carlo simulation of the light environment of plants. Functional
#'Plant Biology. 35:837-849.
#'
#'See also \url{http://www.remkoduursma/yplantqmc}
#'@keywords misc
#'@export
installQuasiMC <- function(proxy=FALSE){

	if(.Platform$OS.type == "windows"){
		
		if(!file.exists("c:/QuasiMC")){
			dir.create("c:/QuasiMC")
			message("Created the directory : c:/QuasiMC")
		}

		if(proxy)utils::setInternet2(TRUE)
		
		if(file.exists("c:/quasimc/QuasiMC.zip")){
			message("Updating existing QuasiMC installation.")
			unlink("c:/quasimc/QuasiMC.zip")
		}
		download.file("http://www.remkoduursma.com/quasimc/QuasiMC.zip",
			"c:/quasimc/QuasiMC.zip")
		unzip("c:/quasimc/QuasiMC.zip", exdir="C:/quasimc")
		message("'QuasiMC.exe' and 'enviro.e' have been placed in the folder 'c:/QuasiMC'.\n  !- Please do not rename or move!\n")

  } else if (Sys.info()[['sysname']] == "Darwin") {

               qmc_dir <- paste(path.expand("~"),"/QuasiMC/",sep="")
               qmc_tgz <- "QuasiMC_Mac.tgz"
               qmc_file <- paste(qmc_dir,qmc_tgz,sep="")

               if(!file.exists(qmc_dir)){
                  dir.create(qmc_dir)
                  message(paste("Created the directory :",qmc_dir))
               }

               if(file.exists(qmc_file)){
                        message("Updating existing QuasiMC installation.")
                        unlink(qmc_file)
               }
               download.file(paste("http://dl.dropbox.com/u/28862459/",qmc_tgz,sep=""), qmc_file)
               untar(qmc_file, exdir=qmc_dir)
               message("'QuasiMC.app' and 'enviro.e' have been placed in the folder '~/QuasiMC'.\n !- Please do not rename or move!\n")
               message("Current QuasiMC.app is for Mac OS >= 10.5")
        } else {
	       stop("QuasiMC is (currently) available for Windows and Mac OS X only.")
	}
}

