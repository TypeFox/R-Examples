################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added attributes to result.
# 28.08.2015: Added importFrom
# 26.08.2014: Fixed bug when scrambled markers (issue#5)
# 27.04.2014: Added option to ignore case in marker names.
# 01.03.2014: Added options 'bins' and calculation of size.
# 01.03.2014: Fixed bug kit always "ESX17".
# 11.09.2014: First version.

#' @title Add Size Information.
#'
#' @description
#' Add size information to alleles.
#'
#' @details
#' Adds a column 'Size' with the fragment size in base pair (bp) for each allele as
#' estimated from kit bins OR calculated from offset and repeat. The bins
#' option return NA for alleles not in bin. The calculate option handles
#' all named alleles including micro variants (e.g. '9.3').
#' Handles 'X' and 'Y' by replacing them with '1' and '2'.
#' 
#' @param data data.frame with at least columns 'Marker' and 'Allele'.
#' @param kit data.frame with columns 'Marker', 'Allele', and 'Size' (for bins=TRUE) or
#'  'Marker', 'Allele', 'Offset' and 'Repeat' (for bins=FALSE).
#' @param bins logical TRUE alleles get size from corresponding bin.
#'  If FALSE the size is calculated from the locus offset and repeat unit.
#' @param ignore.case logical TRUE case in marker names are ignored.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with additional columns for added size.
#' 
#' @export
#' 
#' @importFrom utils head str
#' 

addSize <- function(data, kit=NA, bins=TRUE, ignore.case=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("PARAMETERS:")
    print("data:")
    print(str(data))
    print(head(data))
    print("kit:")
    print(str(kit))
    print(head(kit))
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset.
  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'",
         call. = TRUE)
  }
  
  if(!"Allele" %in% names(data)){
    stop("'data' must contain a column 'Allele'",
         call. = TRUE)
  }
  
  # Check kit depending on 'bins'.
  if(bins){
    
    if(!"Size" %in% names(kit)){
      stop("'kit' must contain a column 'Size'",
           call. = TRUE)
    }
    
    if(!"Allele" %in% names(kit)){
      stop("'kit' must contain a column 'Allele'",
           call. = TRUE)
    }
    
  } else {
    
    if(!"Offset" %in% names(kit)){
      stop("'kit' must contain a column 'Offset'",
           call. = TRUE)
    }
    
    if(!"Repeat" %in% names(kit)){
      stop("'kit' must contain a column 'Repeat'",
           call. = TRUE)
    }
    
  }
  
  # Check kit.
  if(!"Marker" %in% names(kit)){
    stop("'kit' must contain a column 'Marker'",
         call. = TRUE)
  }
  
  # Check if character data.
  if(!is.character(data$Allele)){
    message("'Allele' must be character. 'data' converted")
    data$Allele <- as.character(data$Allele)
  }
  
  # PREPARE -----------------------------------------------------------------
  
  # Check for column 'Size'
  if("Size" %in% names(data)){
    
    message(paste("'data' already contain a column 'Size'\n",
                  "Size will be overwritten!"),
         call. = TRUE)
    
  }
  
  # Add a column 'Size'
  data$Size <- NA

  # Get markers in dataset.
  marker <- unique(data$Marker)
  
  # Check if case in marker names should be ignored.
  if(ignore.case){
    kit$Marker <- toupper(kit$Marker)
    marker <- toupper(marker)
  }
  
  # ADD SIZE ------------------------------------------------------------------

  if(debug){
    print("Markers:")
    print(marker)
  }

  # Loop over markers.
  for(m in seq(along=marker)){
    
    if(debug){
      print("Current marker:")
      print(marker[m])
    }
    
    if(ignore.case){

      # Select rows for current marker and ignore case in marker names.
      cMarker <- toupper(data$Marker) == toupper(marker[m])
      
    } else {
      
      # Select rows for current marker.
      cMarker <- data$Marker == marker[m]
      
    }

    # Get alleles for current marker.
    allele <- unique(data$Allele[cMarker])
    # Remove NAs.
    allele <- allele[!is.na(allele)]

    if(debug){
      print("Alleles:")
      print(allele)
    }
    
    # Loop over allele in current marker.
    for(a in seq(along=allele)){
      
      # Select rows for current allele.
      cAllele <- data$Allele == allele[a]
      
      # Combine selections.
      selection <- cMarker & cAllele
      
      if(bins){
        
        # Get size from matching bins.
        size <- kit$Size[kit$Marker == marker[m] & kit$Allele == allele[a]]
        
      } else {
        # Calculate size from 'offset' and 'repeat'.
        
        # Copy to temporary working variable.
        alleleTmp <- toupper(allele[a])
        
        # Check presence of X/Y.
        if ("X" %in% alleleTmp || "Y" %in% alleleTmp) {
          
          # Use 1 and 2 for X and Y.
          alleleTmp <- sub(pattern="X", replacement=1, x=alleleTmp)
          alleleTmp <- sub(pattern="Y", replacement=2, x=alleleTmp)
        }
        
        # Convert to numeric.
        alleleTmp <- as.numeric(alleleTmp)
        
        # Calculate estimated size.
        tmpOffset <- kit$Offset[kit$Marker == marker[m]]
        tmpRepeat <- kit$Repeat[kit$Marker == marker[m]]
        size <- tmpOffset + floor(alleleTmp) * tmpRepeat + (alleleTmp %% 1) * 10
        
      }
      
      # Store size for current allele in current marker.
      if(length(size) != 0){
        
        data$Size[selection] <- size
        
      } else {
        
        message(paste("Allele", allele[a],
                      "for marker", marker[m], "not in kit definition file."))
        
      }
      
    }
    
  }
  
  if(debug){
    print("Return:")
    print(str(data))
    print(head(data))
  }
  
  # Add attributes to result.
  attr(data, which="slim, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(data, which="slim, call") <- match.call()
  attr(data, which="slim, date") <- date()
  attr(data, which="slim, data") <- substitute(data)
  attr(data, which="slim, kit") <- substitute(kit)
  attr(data, which="slim, bins") <- bins
  attr(data, which="slim, ignore.case") <- ignore.case

  return(data)
  
}
