WriteMspFile <- function(spectra,
                         metadata,
                         filename = "library.msp",
                         comment = "") {

  zz <- file(description = filename, open = "w")

  spectra.names <- unique(spectra$filename)

  for(i in 1:length(spectra.names)) {

    spectrum.name <- as.character(spectra.names[i])

    spectrum <- droplevels(spectra[spectra$filename == spectrum.name, ])

    metadat <- metadata[metadata$filename == spectrum.name, ]

    cat(paste("NAME: ", na.omit(metadat$compound), "\r\n", sep = ""), file = zz)

    cat(paste("COMMENT: ", comment, "\r\n", sep = ""), file = zz)

    cat(paste("Num Peaks: ", nrow(spectrum), "\r\n", sep = ""), file = zz)

    cat(paste(paste(spectrum$mz, spectrum$intensity, sep = " "), ";", sep = ""),
        file = zz, sep = "\r\n")    
    
    cat("\r\n", file = zz)

  }
  
  close(zz)

}
