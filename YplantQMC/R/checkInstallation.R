checkInstallation <- function(){

  if(.Platform$OS.type == "windows"){
	chk <- file.exists("c:/QuasiMC/enviro.e") & file.exists("c:/QuasiMC/QuasiMC.exe")
  }
  else if (Sys.info()[['sysname']] == "Darwin") {
        qmc_dir = paste(path.expand("~"),"/QuasiMC/",sep="")
        chk <- file.exists(paste(qmc_dir,"enviro.e",sep="")) & file.exists(paste(qmc_dir,"QuasiMC.app",sep=""))
  }
  else {
        return(FALSE)
  }
return(chk)
}
