################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 02.12.2016: Possible to return multiple kits by specifying a vector.
# 29.08.2015: Added importFrom.
# 28.06.2015: Changed parameter names to format: lower.case
# 14.12.2014: what='Gender' changed to 'Sex.Marker' now return vector.
# 26.09.2014: Fixed error if kit=NULL and what!=NA.
# 26.08.2014: what=Offset/Repeat, now returns identical data frames.
# 03.08.2014: Added option to return kit index.
# 02.03.2014: Removed factor levels from 'Marker' before returning 'OFFSET'/'REPEAT'.
# 09.12.2013: Removed factor levels from 'Marker' before returning 'VIRTUAL'.
# 20.11.2013: Change parameter name 'kitNameOrIndex' to 'kit'.
# 10.11.2013: 'Marker' returns vector instead of factor.
# 24.10.2013: Fixed error when no matching kit and 'what'!=NA, return NA.
# 04.10.2013: Removed factor levels from 'Marker' before returning 'COLOR'.
# 17.09.2013: Added new parameter 'what' to specify return values.
# 16.09.2013: Changed to support new kits file structure.
# 05.06.2013: Added 'gender.marker'
# 19.05.2013: Re-written for reading data from text file.

#' @title Get Kit
#'
#' @description
#' Provides information about STR kits.
#'
#' @details
#' The function returns the following information for a kit specified in kits.txt:
#' Panel name, short kit name (unique, user defined), full kit name (user defined),
#' marker names, allele names, allele sizes (bp), 
#' minimum allele size, maximum allele size (bp), flag for virtual alleles,
#' marker color, marker repeat unit size (bp), minimum marker size, 
#' maximum marker, marker offset (bp), flag for sex markers (TRUE/FALSE).
#' 
#' If no matching kit or kit index is found NA is returned.
#' If 'NULL' or '0' a vector of available kits is printed and NA returned.
#' 
#' @param kit string or integer to specify the kit.
#' @param what string to specify which information to return. Default is 'NA' which return all info.
#' Not case sensitive.
#' @param show.messages logical, default TRUE for printing messages to the R promt.
#' @param .kit.info data frame, run function on a data frame instead of the kits.txt file.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with kit information.
#' 
#' @export
#' 
#' @importFrom utils read.delim
#' 
#' @examples
#' # Show all information stored for kit with short name 'ESX17'.
#' getKit("ESX17")

getKit<-function(kit=NULL, what=NA, show.messages=FALSE, .kit.info=NULL, debug=FALSE) {

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  .separator <- .Platform$file.sep # Platform dependent path separator.
  
  # LOAD KIT INFO  ############################################################

  if(is.null(.kit.info)){
    # Get package path.
    packagePath <- path.package("strvalidator", quiet = FALSE)
    subFolder <- "extdata"
    fileName <- "kit.txt"
    
    filePath <- paste(packagePath, subFolder, fileName, sep=.separator)
    
    .kit.info <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",
                           dec = ".", fill = TRUE, stringsAsFactors=FALSE)
    
  }

  # Available kits. Must match else if construct.
  kits<-unique(.kit.info$Short.Name)
  
	# Check if NULL
	if (is.null(kit)) {

		# Print available kits
		if (show.messages){
			message("Available kits:")
		}
		res<-kits

	# String provided.
	} else {

		# Check if number or string.
		if (is.numeric(kit)) {

			# Set index to number.
			index<-kit

		} else {

			# Find matching kit index (case insensitive)
			index<-match(toupper(kit),toupper(kits))

		}

		# No matching kit.
		if (any(is.na(index))) {
			
			# Print available kits
			if (show.messages){
				message(paste("No matching kit! \nAvailable kits:", 
                      paste(kits, collapse=", ")))
			}
			return(NA)

		# Assign matching kit information.
		} else {
		  
		  currentKit <- .kit.info[.kit.info$Short.Name %in% kits[index], ]
      
      res <- data.frame(Panel = currentKit$Panel,
                        Short.Name = currentKit$Short.Name,
                        Full.Name = currentKit$Full.Name,
                        Marker = currentKit$Marker,
                        Allele = currentKit$Allele,
                        Size = currentKit$Size,
                        Size.Min = currentKit$Size.Min,
                        Size.Max = currentKit$Size.Max,
                        Virtual = currentKit$Virtual,
                        Color = currentKit$Color,
                        Repeat = currentKit$Repeat,
                        Marker.Min = currentKit$Marker.Min,
                        Marker.Max = currentKit$Marker.Max,
                        Offset = currentKit$Offset,
                        Sex.Marker = currentKit$Sex.Marker,
                        stringsAsFactors = FALSE)
      
      # Create useful factors.
		  res$Marker <- factor(res$Marker, levels=unique(res$Marker))
		  
		} 

	}
  
  # Used in error message in 'else'.
  options <- paste("Index",
                   "Panel",
                   "Short.Name",
                   "Full.Name",
                   "Marker",
                   "Allele",
                   "Size",
                   "Virtual",
                   "Color",
                   "Repeat",
                   "Range",
                   "Offset",
                   "Sex.Marker", sep=", ")
  
  # WHAT ----------------------------------------------------------------------
  
  # Kit is required.
  if (!is.null(kit)) {
    
    if(is.na(what)){
      # Return all kit information.
      return(res)
      
    } else if (toupper(what) == "INDEX"){
      # Return kit index.
      return(index)
      
    } else if (toupper(what) == "PANEL"){
      # Return panel name.
      return(unique(res$Panel))
      
    } else if (toupper(what) == "SHORT.NAME"){
      # Return short name.
      return(unique(res$Short.Name))
      
    } else if (toupper(what) == "FULL.NAME"){
      # Return full name.
      return(unique(res$Full.Name))
      
    } else if (toupper(what) == "MARKER"){
      # Return all markers.
      return(as.vector(unique(res$Marker)))
      
    } else if (toupper(what) == "ALLELE"){
      # Return all alleles and markers.
      
      res <- data.frame(Marker=res$Marker, Allele=res$Allele)
      
      return(res)
      
    } else if (toupper(what) == "SIZE"){
      # Returns all alleles and their indicated normal size in base pair.
      # Their normal size range is idicated in min and max columns.
      # Grouped by marker.
      
      res <- data.frame(Marker=res$Marker,
                        Allele=res$Allele,
                        Size=res$Size,
                        Size.Min=res$Size.Min,
                        Size.Max=res$Size.Max,
                        stringsAsFactors=FALSE)
      
      return(res)
      
    } else if (toupper(what) == "VIRTUAL"){
      # Returns all alleles (bins) with a flag if it is virtual
      # 1 for virtual or 0 it it is a physical ladder fragment.
      # Grouped per marker.
      
      res <- data.frame(Marker=as.character(res$Marker),
                        Allele=res$Allele,
                        Virtual=res$Virtual,
                        stringsAsFactors=FALSE)
      
      return(res)
      
    } else if (toupper(what) == "COLOR"){
      # Return markers and their color as strings.
      
      marker <- getKit(kit, what="Marker")
      color <- NA
      
      for(m in seq(along=marker)){
        color[m] <- unique(res$Color[res$Marker == marker[m]])
      }
      
      res <- data.frame(Marker=marker,
                        Color=color,
                        stringsAsFactors=FALSE)
      
      return(res)
      
    } else if (toupper(what) == "REPEAT"){
      # Return markers and their repeat unit length in base pair.
      
      marker <- getKit(kit, what="Marker")
      offset <- NA
      repeatUnit <- NA
      
      for(m in seq(along=marker)){
        offset[m] <- unique(res$Offset[res$Marker == marker[m]])
        repeatUnit[m] <- unique(res$Repeat[res$Marker == marker[m]])
      }
      
      res <- data.frame(Marker=marker, Offset=offset, Repeat=repeatUnit,
                        stringsAsFactors=FALSE)
      
      return(res)
      
    } else if (toupper(what) == "RANGE"){
      # Return markers and their range (min and max) in base pair.
      
      marker <- getKit(kit, what="Marker")
      markerMin <- NA
      markerMax <- NA
      color <- NA
      
      for(m in seq(along=marker)){
        markerMin[m] <- unique(res$Marker.Min[res$Marker == marker[m]])
        markerMax[m] <- unique(res$Marker.Max[res$Marker == marker[m]])
        color[m] <- unique(res$Color[res$Marker == marker[m]])
      }
      
      res <- data.frame(Marker=marker,
                        Color=color,
                        Marker.Min=markerMin,
                        Marker.Max=markerMax,
                        stringsAsFactors=FALSE)
      
      # Create useful factors.
      res$Color <- factor(res$Color, levels=unique(res$Color))
      
      
      return(res)
      
    } else if (toupper(what) == "OFFSET"){
      # Return markers and their estimated offset in base pair.
      
      marker <- getKit(kit, what="Marker")
      offset <- NA
      repeatUnit <- NA
      
      for(m in seq(along=marker)){
        offset[m] <- unique(res$Offset[res$Marker == marker[m]])
        repeatUnit[m] <- unique(res$Repeat[res$Marker == marker[m]])
      }
      
      res <- data.frame(Marker=marker, Offset=offset, Repeat=repeatUnit,
                        stringsAsFactors=FALSE)
      
      return(res)
      
    } else if (toupper(what) == "SEX.MARKER"){
      # Return sex markers as vector.
      
      sexMarkers <- as.character(unique(res$Marker[res$Sex.Marker == TRUE]))
      
      return(sexMarkers)
      
    } else {
      
      warning(paste(what, "not supported! \nwhat = {", options,"}"))
      return(NA)
      
    }
    
  } else {
    # If kit is NULL return available kits.
    
    return(res)
    
  }

}
