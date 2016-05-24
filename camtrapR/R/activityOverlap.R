activityOverlap <- function(recordTable,
                            speciesA,
                            speciesB,
                            speciesCol = "Species",
                            recordDateTimeCol = "DateTimeOriginal",
                            recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                            plotR = TRUE,
                            writePNG = FALSE,
                            addLegend = TRUE,
                            legendPosition = "topleft",
                            plotDirectory,
                            createDir = FALSE,
                            pngMaxPix = 1000,
                            add.rug = TRUE,
                            ...){

  wd0 <- getwd()
  mar0 <- par()$mar
  on.exit(setwd(wd0))
  on.exit(par(mar = mar0), add = TRUE)

  stopifnot(is.logical(c(writePNG, plotR, createDir)))
  stopifnot(hasArg(speciesA))
  stopifnot(hasArg(speciesB))
  stopifnot(all(c(speciesA, speciesB) %in% recordTable[,speciesCol]))

  recordTable$DateTime2 <- strptime(recordTable[,recordDateTimeCol], format = recordDateTimeFormat, tz = "UTC")
  if("POSIXlt" %in% class(recordTable$DateTime2) == FALSE) stop("couldn't interpret recordDateTimeCol of recordTable using specified recordDateTimeFormat")
  if(any(is.na(recordTable$DateTime2))) stop("at least 1 entry in recordDateTimeCol of recordTable could not be interpreted using recordDateTimeFormat")
  recordTable$Date2 <- as.Date(recordTable$DateTime2)
  recordTable$Time2 <-   format(recordTable$DateTime2, format = "%H:%M:%S", usetz = FALSE)

  # convert time to radians
  recordTable$Time.rad <- (as.numeric(as.POSIXct(strptime(recordTable$Time2, format = "%H:%M:%S"))) -
                             as.numeric(as.POSIXct(strptime("0", format = "%S")))) / 3600 * (pi/12)

  subset_speciesA <- subset(recordTable, recordTable[,speciesCol] == speciesA)
  subset_speciesB <- subset(recordTable, recordTable[,speciesCol] == speciesB)

  if(any(is.na( c(subset_speciesA$Time.rad, subset_speciesB$Time.rad))))stop("NAs produced in converting to radial time. Sorry, that's a bug. Please report it.")

  n_record_string <- paste("number of records:", paste(nrow(subset_speciesA), nrow(subset_speciesB), sep = " / "))
  
  # set graphics  parameters and out directory
  overlapEst.tmp <- overlap::overlapEst(A = subset_speciesA$Time.rad, B = subset_speciesB$Time.rad)

  dhat.tmp <- paste(names(overlapEst.tmp)[1],
                    round(overlapEst.tmp, digits = 2)[1], sep = "=")
  cex.sub <- 0.75
  pngWidth <- pngMaxPix
  pngHeight <- round(pngMaxPix * 0.8)
  mar.tmp <-  c(7.1, 4.1, 4.1, 2.1)
  y_usr_factor <- 0.22

  dots <- list(...)
  if(!is.null(dots[['linetype']])){ lty.tmp <- dots[['linetype']]}  else {lty.tmp <- c(1, 2)}
  if(!is.null(dots[['linewidth']])){lwd.tmp <- dots[['linewidth']]} else {lwd.tmp <- c(1, 1)}
  if(!is.null(dots[['linecol']])){  col.tmp <- dots[['linecol']]}   else {col.tmp <- c("black", "blue")}


  if(isTRUE(writePNG)){
    if(is.null(plotDirectory) == FALSE){
      if(isTRUE(createDir)){
        dir.create(plotDirectory, recursive = TRUE, showWarnings = FALSE)
        setwd(plotDirectory)
      } else {
        stopifnot(file.exists(plotDirectory))
        setwd(plotDirectory)
      }
    } else {stop("writePNG is TRUE. Please set plotDirectory")}
  }

  if(isTRUE(writePNG)){

    png(filename = paste("activity_overlap_", speciesA, "-", speciesB, "_", Sys.Date(), ".png", sep = ""),
        width = pngWidth, height = pngHeight, units = "px", res = 96, type = "cairo")

    par(mar = mar.tmp)
    overlapPlot(A = subset_speciesA$Time.rad, B = subset_speciesB$Time.rad,
                main = paste("Activity overlap: ", speciesA, "-", speciesB),
                rug = add.rug,
                ...)
    legend( x = "top",
            legend = dhat.tmp,
            bty = "n",
            cex = cex.sub)

    if(isTRUE(addLegend)){
      legend(x = legendPosition,
             legend = c(speciesA, speciesB),
             lty = lty.tmp,
             lwd = lwd.tmp,
             col = col.tmp,
             bg = "white"
      )
    }
    mtext(n_record_string, side = 3, line = 0)            
    dev.off()
  }

  if(isTRUE(plotR)){

    par(mar = mar.tmp)
    plot.values <- overlapPlot(A = subset_speciesA$Time.rad, B = subset_speciesB$Time.rad,
                               main = paste("Activity overlap: ", speciesA, "-", speciesB),
                               rug = add.rug,
                               ... )

    legend(x = "top",
           legend = dhat.tmp,
           bty = "n",
           cex = cex.sub)

    if(isTRUE(addLegend)){
      legend(x = legendPosition,
             legend = c(speciesA, speciesB),
             lty = lty.tmp,
             lwd = lwd.tmp,
             col = col.tmp,
             bg = "white"
      )
    }
    mtext(n_record_string, side = 3, line = 0)            
    
    return(invisible(plot.values))
  }
}