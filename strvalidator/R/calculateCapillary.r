################################################################################
# TODO LIST
# TODO: Hard-coded plate/injection information, should be put in a text file instead.

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 15.12.2014: Changed parameter names to format: lower.case
# 28.10.2013: First version.

#' @title Calculate Capillary Balance
#'
#' @description
#' Calculates the ILS inter capillary balance.
#'
#' @details
#' Calculates the inter capillary balance for the internal lane standard (ILS).
#' Require information from both the 'samples.table' and the 'plot.table'.
#' 
#' @param samples.table data frame containing at least the columns
#'  'Sample.File', 'Sample.Name', 'Size.Standard', 'Instrument.Type',
#'  'Instrument.ID', 'Cap', 'Well', and 'SQ'.
#' @param plot.table data frame containing at least the columns
#'  'Sample.File.Name', 'Size', and 'Height'.
#' @param sq numeric threshold for 'Sizing Quality' (SQ).
#' @param run character string for run name.
#' @param debug logical indicating printing debug information.
#'  
#' @return data.frame with with columns 'Instrument', 'Instrument.ID', 'Run',
#' 'Mean.Height', 'SQ', 'Injection', 'Capillary', 'Well', 'Comment'.
#' 
#' @export
#' 
#' @importFrom utils read.delim help head tail str
#' @importFrom graphics title
#' 

calculateCapillary <- function(samples.table, plot.table, sq=0, run="", debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Check samples.table ---------------------------------------------------------
  
  if(is.null(samples.table$Sample.File)){
    stop("'Sample.File' does not exist!")
  }
  
  if(is.null(samples.table$Sample.Name)){
    stop("'Sample.Name' does not exist!")
  }

  if(is.null(samples.table$Size.Standard)){
    stop("'Size.Standard' does not exist!")
  }
  
  if(is.null(samples.table$Instrument.Type)){
    stop("'Instrument.Type' does not exist!")
  }
  
  if(is.null(samples.table$Instrument.ID)){
    stop("'Instrument.ID' does not exist!")
  }
  
  if(is.null(samples.table$Cap)){
    stop("'Cap' does not exist!")
  }

  if(is.null(samples.table$Well)){
    stop("'Well' does not exist!")
  }

  if(is.null(samples.table$SQ)){
    stop("'SQ' does not exist!")
  }

  # Check plot.table ------------------------------------------------------------
  
  if(is.null(plot.table$Sample.File.Name)){
    stop("'Sample.File.Name' does not exist!")
  }
  
  if(is.null(plot.table$Size)){
    stop("'Size' does not exist!")
  }
  
  if(is.null(plot.table$Height)){
    stop("'Height' does not exist!")
  }
  
#  if(is.null(plot.table$Data.Point)){
#    stop("'Data.Point' does not exist!")
#  }
  
  # Read size standard --------------------------------------------------------

  .separator <- .Platform$file.sep # Platform dependent path separator.

  ilsName <- unique(samples.table$Size.Standard)
  
  if(!is.null(ilsName) && length(ilsName) == 1){
    # Get package path.
    packagePath <- path.package("strvalidator", quiet = FALSE)
    subFolder <- "extdata"
    fileName <- "ils.txt"
    
    filePath <- paste(packagePath, subFolder, fileName, sep=.separator)
    
    ils <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",
                           dec = ".", fill = TRUE, stringsAsFactors=FALSE)
    
    if(debug){
      print(paste("ILS definition loaded from",filePath))
      print(str(ils))
      print(head(ils))
      print(tail(ils))
    }
    
  } else {
    stop("Size standard name error: not found, or multiple.")
  }

  # Merge ---------------------------------------------------------------------
  
  # "Sample.File.Name" must be changed to "Sample.File" before mergin.
  # NB! Sample.File is the only common information.
  # NB! This mean that it is not possible to distinguish between re-injections (runs)!
  # NB! The only option if re-injections are used is to use positional information.
  # NB!  i.e. firs re-injection comes first followed by the second and so on...
  if(!"Sample.File" %in% names(plot.table)){
    
    if("Sample.File.Name" %in% names(plot.table)){
      
      if(debug){
        print("Change column names from:")
        print(names(plot.table))
      }

      # Change column name.
      colnames(plot.table)[which(colnames(plot.table) == "Sample.File.Name")] <- "Sample.File"
      
      if(debug){
        print("To:")
        print(names(plot.table))
      }
      
    }
    
  }
  
  # Merge data frames.  
  df <- merge(plot.table, samples.table, by="Sample.File")
  
  # Prepare -------------------------------------------------------------------
  
  # Check data type of Height.
  if(is.character(df$Height)){
    message("'Height' is character. Converting to numeric.")
    # Convert to numeric.
    df$Height <- as.numeric(df$Height)
  }

  # Check data type of Height.
  if(is.character(df$Size)){
    message("'Size' is character. Converting to numeric.")
    # Convert to numeric.
    df$Size <- as.numeric(df$Size)
  }

  # Add information -----------------------------------------------------------
  
  # Get instrument type.
  instrType <- unique(df$Instrument.Type)
  instrument <- ifelse(length(instrType) == 1, instrType, NA)

  if(debug){
    print("Instrument type:")
    print(instrument)
  }
  
  # Get specific information for the instrument type.
  # NB! Hard-coded information, should be put in a text file instead.
  
  # Lane number (plate position) and well is plate specific (her only 96-well).
  lane <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96)
  well <- c("A01","B01","C01","D01","E01","F01","G01","H01","A02","B02","C02","D02","E02","F02","G02","H02","A03","B03","C03","D03","E03","F03","G03","H03","A04","B04","C04","D04","E04","F04","G04","H04","A05","B05","C05","D05","E05","F05","G05","H05","A06","B06","C06","D06","E06","F06","G06","H06","A07","B07","C07","D07","E07","F07","G07","H07","A08","B08","C08","D08","E08","F08","G08","H08","A09","B09","C09","D09","E09","F09","G09","H09","A10","B10","C10","D10","E10","F10","G10","H10","A11","B11","C11","D11","E11","F11","G11","H11","A12","B12","C12","D12","E12","F12","G12","H12")
  
  # Capillary and injection is instrument dependent.
  if(instrument == "ABI3500"){
    
    capillary <- c(1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24,1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24,1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24,1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24)
    injection <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
    
  } else if(instrument == "ABI3130") {
    
    capillary <- c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16)
    injection <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)
  
  } else {
    
    stop(paste("Specific information not available for instrument of type:", instrument))
  }
  
  # Add extra information by looping over each well.
  for(w in seq(along=well)){
    
    df$Injection[df$Well == well[w]] <- injection[w]   # This is the important information.
    df$Lane[df$Well == well[w]] <- lane[w]             # Equals the well number when reading in column order (1-12).
    # df$Capillary[df$Well == well[w]] <- capillary[w] # Can be used for checking, should be identical with 'Cap'.
    
  }

  if(debug){
    print("Extra information added:")
    print(head(df))
  }
  
  # Nothing to separeate one run (i.e. 4 injections of the same plate run in sequence.)
  # Solution is to create separate projects in GMIDX and export on run per file.
  # Separate per capillary instead.
  
  # Declare temporary vector variables.
  vecHeight <- vector()
  vecSQ <- vector()
  vecInj <- vector()
  vecCap <- vector()
  vecWell <- vector()
  vecCom <- vector()
  
  # Initiate row/lane index.
  r <- 1
  
  # Get instrument type.
  id <- unique(df$Instrument.ID)
  if(length(id) == 1){
    instrumentID <- id
  } else {
    instrumentID <- NA
    message("Multiple instrument ID's present. Will use 'NA'.")
  }
  
  # Get all capillaries.
  cap <- unique(df$Cap)
  
  # Loop over all capillaries.
  for(c in seq(along=cap)){
    
    # Make selection for current capillary.
    selCap <- df$Cap == cap[c]
    
    # Make selection for current capillary.
    inj <- unique(df$Injection)
    
    # Loop over all injections.
    for(i in seq(along=inj)){
      
      # Make selection for current injection.
      selInj <- df$Injection == inj[i]
      
      # Make combined selection for current capillary and injection.
      selection <- selCap & selInj
      
      # Check if empty selection.
      if(sum(selection) != 0){
        
        # Check current sizing quality.
        csq <- unique(df[selection,]$SQ)
        
        # Expect only one value.
        if(length(csq) != 1){
          stop(paste("Error in dataset or selection",
                     "\nCurrent capillary:", cap[c],
                     "\nCurrent injection:", inj[i],
                     "\nSQ:", paste(csq, collapse="")))      
        }
        
        # Check if sizing quality pass threshold.
        # Default is 0, everything accepted.
        if(!(csq < sq)){
          
          # Select one sample.
          # Sum peaks, remember to filter first!
          
          # Get peak sizes.
          pz <- df[selection,]$Size
          
          # Get peak heights.
          ph <- df[selection,]$Height
          
          if(debug){
            print("Peak size:")
            print(pz)
            print("Peak heights:")
            print(ph)
          }
          
          # Filter peak heights not in ILS.
          ph <- ph[pz %in% ils$sizeStdDefinition]
          
          if(debug){
            print("Filtered peak heights:")
            print(ph)
          }
          
          # Calculate mean and save information in tmp vectors.
          vecHeight[r] <- mean(ph)
          vecSQ[r] <- unique(df[selection,]$SQ)
          vecInj[r] <- unique(df[selection,]$Injection)
          vecCap[r] <- unique(df[selection,]$Cap)
          vecWell[r] <- unique(df[selection,]$Well)
          vecCom[r] <- NA
          
        } else {
          # Size quality below threshold.
          
          vecHeight[r] <- NA
          vecSQ[r] <- unique(df[selection,]$SQ)
          vecInj[r] <- unique(df[selection,]$Injection)
          vecCap[r] <- unique(df[selection,]$Cap)
          vecWell[r] <- unique(df[selection,]$Well)
          vecCom[r] <- paste("SQ <", sq)
          
        }
        
      } else {
        # Data is missing.
        
        if(debug){
          print(paste("Data is missing for capillary", cap[c], "injection", inj[i]))
        }
        
        vecHeight[r] <- NA
        vecSQ[r] <- NA
        vecInj[r] <- inj[i]
        vecCap[r] <- cap[c]
        vecWell[r] <- well[r]
        vecCom[r] <- "Missing"
        
      }
      
      # Increase vector index.
      r <- r + 1
      
    }
    
  }

  # Create result data frame.
  res <- data.frame(Instrument=instrument,
                    Instrument.ID=instrumentID,
                    Run=run,
                    Mean.Height=vecHeight,
                    SQ=vecSQ,
                    Injection=vecInj,
                    Capillary=vecCap,
                    Well=vecWell,
                    Comment=vecCom,
                    stringsAsFactors=FALSE)
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(res)
  
}
