# copyright (C) 2014-2016 A.Rebecq

# TODO : add xtable parameters to marginsToTeX parameters
marginsToTeX <- function(marginMatrix, names=NULL, pct=FALSE, popTotal=NULL,
                         scaleboxTeX=NULL, file=NULL,
                         label=NULL, caption=NULL) {
  
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Package xtable needed for export of margins in LateX to work. Please install it.",
         call. = FALSE)
  }
  
  if(!is.matrix(marginMatrix)) {
    stop("marginsToTeX input type has to be matrix.")
  }
  
  if(!is.null(names)) {
    if(length(names) != nrow(marginMatrix)) {
      stop("Name length must equal number of rows in marginMatrix")
    }
    
    marginMatrix[,1] <- names
  }
  
  if(pct) {
    numericPart <- marginMatrix[,3:(ncol(marginMatrix))]
    numericPart <- as.numeric(numericPart)
    numericPart <- 100*numericPart
    numericPart -> marginMatrix[,3:(ncol(marginMatrix))]
  }
  
  # Write zeros as NA
  marginMatrix[as.numeric(marginMatrix) == 0] <- NA

  marginDF <- as.data.frame(marginMatrix)
  
  # Heuristic rule for scalebox
  if(is.null(scaleboxTeX)) {
    if(ncol(marginDF) >= 10) {
      scaleboxTeX <- 1.4 - ncol(marginDF) / 20
    }
    
    if(ncol(marginDF) >= 28) {
      stop("Automatic scaleboxing not configured for more than 28 margins.")
    } 
  }
  
  captionTeX <- caption
  if(!is.null(popTotal)) {
    captionTeX <- paste(caption, " -- total population : ", round(popTotal,0),sep="")
  }
  
  print(xtable::xtable(marginDF, caption=captionTeX, label=label), include.rownames = FALSE, include.colnames = FALSE,
               floating = TRUE, scalebox=scaleboxTeX, file=file
        )

}