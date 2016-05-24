################################################################################
# TODO LIST
# TODO: Allow to use Data.Point

################################################################################
# CHANGE LOG (last 20 changes)
# 12.10.2015: First version.

#' @title Detect Spike
#'
#' @description
#' Detect samples with possible spikes in the DNA profile.
#'
#' @details Creates a list of possible spikes by searching for peaks aligned
#' vertically (i.e. nearly identical size).
#' 
#' @param data data.frame with including colums 'Sample.Name', 'Marker', 'Size'.
#' @param threshold numeric number of peaks of similar size in different dye
#' channels to pass as a possible spike (NULL = number of dye channels
#' minus one to allow for on unlabelled peak).
#' @param round.to numeric toleranse for Size.
#' @param kit string or numeric for the STR-kit used (NULL = auto detect).
#' @param debug logical indicating printing debug information.
#' 
#' @export
#' 
#' @importFrom data.table data.table := .N
#' 
#' @return data.frame
#' 
#' @seealso \code{\link{data.table}}


calculateSpike <- function(data, threshold=NULL, round.to=1, kit=NULL, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("str(data):")
    print(str(data))
    print("threshold:")
    print(threshold)
    print("round.to:")
    print(round.to)
    print("kit:")
    print(kit)
  }

  # Check data ----------------------------------------------------------------
  
  # Columns.  
  if(is.null(data$Sample.Name)){
    stop("'Sample.Name' does not exist!")
  }
  if(is.null(data$Marker)){
    stop("'Marker' does not exist!")
  }
  if(is.null(data$Size)){
    stop("'Size' does not exist!")
  }
  
  # Check if slim format.  
  if(sum(grepl("Size", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }

  # Check data type.  
  if(!is.numeric(data$Size)){
    data$Size <- as.numeric(data$Size)
    warning("'Size' not numeric! 'data' converted.")
  }
  
  # Prepare -------------------------------------------------------------------
  
  # Convert to data.table.  
  dtable <- data.table::data.table(data)

  if(is.null(kit)){

    # Detect kit.
    kit <- detectKit(data=dtable, index=FALSE, debug=debug)
    kit <- kit[1]
    message(paste("Using kit:", kit))
    
  }

  # Getdye channels.
  kitColors <- unique(getKit(kit=kit, what="Color", debug=debug)$Color)
  kitDyes <- addColor(data=kitColors, have="Color", need="Dye")
  
  if(is.null(threshold)){

    # Default to number of dyes minus one to allow for one unlabelled spike.
    threshold <- length(kitDyes) - 1
    message(paste("Using default spike threshold:", threshold))
    
  } else if (threshold > length(kitDyes)){
    
    # Threshold cannot be larger than the number of dyes in the kit.
    threshold <- length(kitDyes) - 1
    message(paste("'threshold' cannot be larger than the number of dyes in the kit"))
    message(paste("Using default spike threshold:", threshold))
    
  }
    
  # Add Dye if not present.
  if(!"Dye" %in% names(dtable)){
    dtable <- addColor(data=dtable, kit=kit, need="Dye")
  }
  
  # Add Id if not present.
  if(!"Id" %in% names(dtable)){
    dtable$Id <- paste(dtable$Sample.Name, dtable$File.Name)
  }
  
  # Round to nearest base pair
  dtable$Round <- round.to * round(dtable$Size / round.to) 
  
  # Remove NA's.
  if(any(is.na(dtable$Size))){
    dtable <- dtable[!is.na(dtable$Size), ]
  }

  # Analyse -------------------------------------------------------------------
  
  # Define column names (to avoid notes.)  
  col1 <- "File.Name"
  col2 <- "Sample.Name"
  col3 <- "Peaks"
  col4 <- "Size"
  
  # Count number of peaks of same size per sample.
  res <- dtable[, "Peaks":=.N, by=c("Id", "Round")]

  # Extract samples with at least 'threshold' peaks of the same size.
  # res <- res["Peaks" >= threshold,]
  res <- res[get(col3) >= threshold,]
  
  # Sort table.
  #res <- res[order("File.Name", "Sample.Name", -"Peaks", "Size")]
  res <- res[order(get(col1), get(col2), -get(col3), get(col4))]
  
  # Convert to data.frame.
  res <- as.data.frame(res)
  
  # Add attributes to result.
  attr(res, which="calculateSpike, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="calculateSpike, call") <- match.call()
  attr(res, which="calculateSpike, date") <- date()

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)
  

}