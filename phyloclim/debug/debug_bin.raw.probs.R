library(Roxalis)
clim <- list.files("/Users/Stoffi/R_files/maxent/layers_Desert",
	pattern = ".asc", full.names = TRUE)
path_bioclim = clim[1]

path_model <- modeldir <- "/Users/Stoffi/R_files/maxent/output_Carnosae"