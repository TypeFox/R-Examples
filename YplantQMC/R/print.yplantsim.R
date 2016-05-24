#'@method print yplantsim
#'@S3method print yplantsim
print.yplantsim <- function(x, ...){
	

	cat("Yplant - simulation result ('yplantsim' object).\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	cat("Plant file: ", x$plant$pfile, "\n")
	cat("Leaf file: ", x$plant$lfile, "\n")
	cat("Number of time steps: ", x$nsteps, "\n\n")
	
	if(x$met$daylength != "not calculated"){
		if("A" %in% names(x$psrdata))
      cat("Total photosynthesis (mol CO2 d-1) =", 
			round(sum(x$psrdata$A * x$psrdata$timestep),3)*10^-6, "\n")
		if("E" %in% names(x$psrdata)){
			cat("Total transpiration (mol H2O d-1) =", 
			round(sum(x$psrdata$E * x$psrdata$timestep),3)*10^-3, "\n\n")
		}
	}
	cat("To view plant totals by timestep, use: psrdata() on your simulation result.\n")
	
}
