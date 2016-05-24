#############################################################
# Export of stopping powers from libamtrack for use in TRiP98
#
# S. Greilich, Jul 2011 / Sep 2011
#############################################################
AT.SPC.export.DEDX <- function(stopping.power.source.no, file.name.DEDX = NULL, element.names = NULL, energy.MeV.u = NULL, plot = TRUE, write = TRUE)
{
	if(is.null(file.name.DEDX)){
		file.name.DEDX    <- "libamtrack.dedx"
	}

	# If no elements given, use ICRU49/73 - and A = 2*Z for simplicity
	if(is.null(element.names)){
		element.names    <- c(	"H",  "He", "Li", "Be", "B", 
								"C",  "N",  "O",  "F",  "Ne", 
								"Na", "Mg", "Al", "Si", "P", 
								"S",  "Cl", "Ar")
		particle.names    <- c(	"1H",   "4He",  "6Li",  "8Be",  "10B", 
								"12C",  "14N",  "16O",  "18F",  "20Ne", 
								"22Na", "24Mg", "26Al", "28Si", "30P", 
								"32S",  "34Cl", "36Ar")
	}else{
	# TODO: Add 2*Z as A to particle.names
	}
	
	particle.nos <- AT.particle.no.from.particle.name(particle.name = particle.names)
	
	# if no energy grid given, use ICRU49/73
	if(is.null(energy.MeV.u)){
		energy.MeV.u  	  <- c( 2.5e-02, 3.0e-02, 4.0e-02, 5.0e-02, 6.0e-02, 7.0e-02, 8.0e-02, 9.0e-02, 1.0e-01, 1.5e-01, 2.0e-01, 2.5e-01, 3.0e-01, 4.0e-01,
								5.0e-01, 6.0e-01, 7.0e-01, 8.0e-01, 9.0e-01, 1.0e+00, 1.5e+00, 2.0e+00, 2.5e+00, 3.0e+00, 4.0e+00, 5.0e+00, 6.0e+00, 7.0e+00,
								8.0e+00, 9.0e+00, 1.0e+01, 1.5e+01, 2.0e+01, 2.5e+01, 3.0e+01, 4.0e+01, 5.0e+01, 6.0e+01, 7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02,
								1.5e+02, 2.0e+02, 2.5e+02, 3.0e+02, 4.0e+02, 5.0e+02, 6.0e+02, 7.0e+02, 8.0e+02, 9.0e+02, 1.0e+03)
	}
	
	# Get stopping power data
	df <- expand.grid( energy.MeV.u      = energy.MeV.u,
	                   particle.name     = particle.names,
					   particle.no       = 0,
					   stopping.power.MeV.cm2.g = 0)

	df$particle.no   <- particle.nos[match(df$particle.name, particle.names)]
  df$stopping.power.MeV.cm2.g <- AT.Mass.Stopping.Power.with.no( stopping.power.source.no = stopping.power.source.no,
                                                              E.MeV.u     = df$energy.MeV.u,
                                                              particle.no = df$particle.no,
                                                              material.no = AT.material.no.from.material.name(material.name = "Water, Liquid"))$Stopping.Power.MeV.cm2.g
	
	output <- "!filetype    dEdx
	!fileversion    19980515
	!filedate    "

	output <- paste(output, date(), "\n", sep = "")
	output <- paste(output, "!material    H2O
	!density    1
	#############################################################
	#Export of SHIELD-HIT Bethe stopping powers for use in TRiP98
	#
	# S. Greilich, Jul 2011
	#############################################################
	\n", sep = "")


	for(cur.particle in unique(df$particle.name)){
		# cur.particle <- unique(df$particle.name)[1]
		ii <- df$particle.name == cur.particle
		output <- paste(output, 
						"!projectile    ", cur.particle, "\n", 
						"# E/(MeV/u) dE/dx(MeVcm**2/g)\n",
						"!dedx\n", 
						sep = "")

		for(i in 1:nrow(df[ii,])){
			output <- paste(output,
							df[ii,]$energy.MeV.u[i], " ",
							df[ii,]$stopping.power.MeV.cm2.g[i], "\n",
							sep = "")
		}
	}

	if(write == TRUE){
		write(output, file = file.name.DEDX)
	}
 
  if(plot == TRUE){
    lattice::xyplot( log10(stopping.power.MeV.cm2.g) ~ log10(energy.MeV.u)|sprintf("data source no. %d", stopping.power.source.no),
                     df,
                     groups   = df$particle.name,
                     type     = 'o',
                     auto.key = list(space = 'right'))
  }

  if(sum(is.na(df$stopping.power.MeV.cm2.g))>0){
    print("Warning: NAs in stopping power data! Check your energy grid.")
  }
  
	return(df)
}
