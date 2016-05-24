# Function to convert export TRiP compatible DDD files
# 
# Author: greilich
# based on original script by botas, greilich, Jul11
###############################################################################
#
# TODO: Transfer target / projectile information to DDD file, currently only energy,
# TODO: the rest is hard-coded. sgre, 2011-09-18
#
AT.SPC.convert.to.DDD <- function(file.name.spc, file.name.ddd = NULL, endian = 'little', plot = TRUE, write = TRUE, ...){

	# Read data
	data           <- AT.SPC.read( file.name = file.name.spc, 
								   endian    = endian,
								   ...)
								  
	if(!is.null(file.name.ddd)){
		file.name.ddd <- gsub(".spc", ".ddd", file.name.spc, fixed = TRUE)
	}
	
	# Convert to dose
	DDD           <- AT.SPC.tapply( spc                  = data$spc, 
									INDEX                = c("depth.g.cm2"), 
									FUN                  = AT.total.D.Gy, 
									additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
																c("stopping.power.source.no", "2", FALSE)),
									names.results        = "D.Gy")
									  
	# Translate to ddd
	conversion.factor <- 6.2415e9
	DDD$D.Gy          <- DDD$D.Gy * conversion.factor

	if(plot == TRUE){
		lattice::xyplot( D.Gy ~ depth.g.cm2,
				DDD,
				type = 'o',
				main = file.name.ddd)
		}

	output <- "!filetype    ddd
	!fileversion    19980520
	!filedate    "
	output <- paste(output, date(), "\n", sep = "")
	output <- paste(output, "!projectile    12C
	!material      H2O
	!composition   H2O
	!density 1
	!energy ",
	round(data$energy.MeV.u, 0),"
	#   DDD file from libamtrack free spc
	#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[g/cm**2] factor FWHM2[g/cm**2]
	!ddd
	", sep = "")

	for(i in 1:nrow(DDD)){
		output <- paste(output,
						" ",
						sprintf("%5.4f", DDD$depth.g.cm2[i]),
						" ",
						sprintf("%5.4f", DDD$D.Gy[i]),
						"\n",
						sep = "")
	}

	if(write == TRUE){
		write(output, file = file.name.ddd)
	}
	
	return(DDD)
}
