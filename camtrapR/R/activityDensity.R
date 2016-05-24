activityDensity <- function(recordTable,
                            species,
                            allSpecies = FALSE,
                            speciesCol = "Species",
                            recordDateTimeCol = "DateTimeOriginal",
                            recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                            plotR = TRUE,
                            writePNG = FALSE,
                            plotDirectory,
                            createDir = FALSE,
                            pngMaxPix = 1000,
                            add.rug = TRUE,
                            ...
)
{

  wd0 <- getwd()
  mar0 <- par()$mar
  on.exit(setwd(wd0))
  on.exit(par(mar = mar0), add = TRUE)

  stopifnot(is.logical(c(allSpecies, writePNG, plotR, createDir)))
  if(allSpecies == FALSE) {
    stopifnot(species %in% recordTable[,speciesCol])
    stopifnot(hasArg(species))
  }

  recordTable$DateTime2 <- strptime(recordTable[,recordDateTimeCol], format = recordDateTimeFormat, tz = "UTC")
  if("POSIXlt" %in% class(recordTable$DateTime2) == FALSE) stop("couldn't interpret recordDateTimeCol of recordTable using specified recordDateTimeFormat")
  if(any(is.na(recordTable$DateTime2))) stop("at least 1 entry in recordDateTimeCol of recordTable could not be interpreted using recordDateTimeFormat")
  recordTable$Date2 <- as.Date(recordTable$DateTime2)
  recordTable$Time2 <- format(recordTable$DateTime2, format = "%H:%M:%S", usetz = FALSE)


  # radians time
  recordTable$Time.rad <- (as.numeric(as.POSIXct(strptime(recordTable$Time2, format = "%H:%M:%S"))) -
                             as.numeric(as.POSIXct(strptime("0", format = "%S")))) / 3600 * (pi/12)

  if(isTRUE(writePNG)){
    if(hasArg(plotDirectory)){
      if(isTRUE(createDir)){
        dir.create(plotDirectory, recursive = TRUE, showWarnings = FALSE)
        setwd(plotDirectory)
      } else {
        stopifnot(file.exists(plotDirectory))
        setwd(plotDirectory)
      }
    } else {
      stop("writePNG is TRUE. Please set plotDirectory")}
  }

  pngWidth <- pngMaxPix
  pngHeight <- round(pngMaxPix * 0.8)

  if(allSpecies == FALSE){

    subset_species <- subset(recordTable, recordTable[,speciesCol] == species)

    if(isTRUE(writePNG)){
      png(filename = paste("activity_density_", species, "_", Sys.Date(), ".png", sep = ""),
          width = pngWidth, height = pngHeight, units = "px", res = 96, type = "cairo")
      densityPlot(subset_species$Time.rad,
                  main = paste("Activity of", species),
                  rug = add.rug,
                  ...)
      mtext(paste("number of records:", nrow(subset_species)), side = 3, line = 0)
      dev.off()
    }
    if(isTRUE(plotR)){
      densityPlot(subset_species$Time.rad,
                  main = paste("Activity of", species),
                  rug = add.rug,
                  ...)
    mtext(paste("number of records:", nrow(subset_species)), side = 3, line = 0)
    }

  } else {

    subset_species_list <- list()

    for(i in 1:length(unique(recordTable[,speciesCol]))){

      spec.tmp <- unique(recordTable[,speciesCol])[i]
      subset_species <- subset(recordTable, recordTable[,speciesCol] == spec.tmp)

      subset_species <- subset(recordTable, recordTable[,speciesCol] == spec.tmp)
      if(nrow(subset_species) == 1){
        warning(paste(spec.tmp, "had only 1 record, cannot estimate density"))
      } else {
        if(isTRUE(writePNG)){
          png(filename = paste("activity_density_", spec.tmp, "_", Sys.Date(), ".png", sep = ""),
              width = pngWidth, height = pngHeight, units = "px", res = 96, type = "cairo")
          densityPlot(subset_species$Time.rad,
                      main = paste("Activity of", spec.tmp),
                      rug = add.rug,
                      ...)
          mtext(paste("number of records:", nrow(subset_species)), side = 3, line = 0)
          dev.off()
        }

        if(isTRUE(plotR)){
          densityPlot(subset_species$Time.rad,
                      main = paste("Activity of", spec.tmp),
                      rug = add.rug,
                      ...)
          mtext(paste("number of records:", nrow(subset_species)), side = 3, line = 0)
        }
      }
      subset_species_list[[i]] <- subset_species$Time.rad
      names(subset_species_list)[i] <- spec.tmp
    }
  }
  if(allSpecies == FALSE){
    return(invisible(subset_species$Time.rad))
  } else {
    return(invisible(subset_species_list))
  }
}