################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 06.02.2015: Fixed error message 'not identical panels'.
# 15.09.2013: First version.

#' @title Combine Bins And Panels Files.
#'
#' @description
#' Combines useful information into one dataset.
#'
#' @details
#' Combines information from two sources ('Bins' and 'Panels' file) to create
#' a dataset containing information about panel name, marker name, alleles in the
#' allelic ladder, their size and size range, a flag indicating virtual alleles,
#' fluorophore color, repeat size, marker range. The short name, full name,
#'  and sex marker flag is populated through the \code{makeKit_gui} user interface.
#' In addition the function calculates an estimated offset for each marker, 
#' which can be used for creating epg-like plots.
#' Note: offset is estimated by taking the smallest physical ladder fragment
#' e.g. 98.28 for D3 in ESX17. Round this to an integer (98) and finally
#' subtract the number of base pair for that repeat i.e. 4*9=36,
#' which gives an offset of 98-36 = 62 bp.
#' Microvariants are handled by taking the decimal part multiplied with 10
#' and adding this to the number of base pair e.g. 9.3 = 4*9 + 0.3*10 = 39 bp.
#' 
#' @param bin data frame created from the 'bins' file.
#' @param panel data frame created from the 'panels' file.
#' 
#' @return data frame containing columns 'Panel', 'Marker', 'Allele',
#' 'Size', 'Size.Min', 'Size.Max', 'Virtual', 'Color', 'Repeat', 
#' 'Marker.Min', 'Marker.Max', 'Offset', 'Short.Name', 'Full.Name'
#' 

combineBinsAndPanels <- function(bin, panel){
  
  kit <- bin
  
  # Add new columns.
  kit$Color<-NA
  kit$Repeat <- NA
  kit$Marker.Min <- NA
  kit$Marker.Max <- NA
  kit$Offset <- NA

  # Get panels.
  binPanel <- unique(bin$Panel)
  binPanel2 <- unique(panel$Panel)
  
  if(!all(binPanel == binPanel2)){
    print(paste("bin panels:", paste(binPanel, collapse=",")))
    print(paste("panel panels:", paste(binPanel2, collapse=",")))
    stop("Panels in 'bin' and 'panel' files not identical")
  }
  
  # Loop over all panels.
  for(p in seq(along=binPanel)){

    # Get markers for current panel.
    binMarker <- unique(bin$Marker[bin$Panel == binPanel[p]])
    
    for(m in seq(along=binMarker)){
      
      # Add new info for current marker in current panel.
      
      # Color.
      kit$Color[kit$Panel == binPanel[p] & kit$Marker == binMarker[m]] <- 
        panel[panel$Panel==binPanel[p] & panel$Marker==binMarker[m],]$Color
      
      # Repeat unit size.
      kit$Repeat[kit$Panel == binPanel[p] & kit$Marker == binMarker[m]] <- 
        panel[panel$Panel==binPanel[p] & panel$Marker==binMarker[m],]$Repeat
      

      # Marker size range min.
      kit$Marker.Min[kit$Panel == binPanel[p] & kit$Marker == binMarker[m]] <- 
        panel[panel$Panel==binPanel[p] & panel$Marker==binMarker[m],]$Marker.Min
      
      # Marker size range max.
      kit$Marker.Max[kit$Panel == binPanel[p] & kit$Marker == binMarker[m]] <- 
        panel[panel$Panel==binPanel[p] & panel$Marker==binMarker[m],]$Marker.Max
      
    }
    
  }
  
  # Estimate marker offset by taking the smallest ladder fragment.
  # Round this to an integer.
  # Subtract the number of base pair for that repeat.

  # Get panels.
  panel <- unique(kit$Panel)
  
  # Loop over all panels.
  for(p in seq(along=panel)){
    
    # Select current panel.
    selPanel <- kit$Panel == panel[p]
    
    # Get markers for current panel.
    marker <- unique(kit$Marker[kit$Panel == panel[p]])

    # Loop over all markers.
    for(m in seq(along=marker)){

      # Select current marker.
      selMarker <- kit$Marker == marker[m]
      
      # Get smallest physical ladder fragment.
      fragments <- kit$Size[selPanel & selMarker & kit$Virtual == 0]
      minFragment <- min(fragments)
      
      # Get corresponsing allele and concert to numeric.
      minAllele <- kit$Allele[selPanel & selMarker & kit$Size == minFragment]
      if(minAllele == "X"){
        minAllele <- 1
      }
      minAllele <- as.numeric(minAllele)
      
      # Get the repeat unit.
      repeatUnit <- kit$Repeat[selPanel & selMarker & kit$Size == minFragment]

      # Calculate offset.
      minFragment <- round(minFragment)
      alleleSize <- floor(minAllele) * repeatUnit + ((minAllele %% 1) * 10)
      markerOffset <- minFragment - alleleSize

      # Add new info for current marker in current panel.
      kit$Offset[selPanel & selMarker] <- markerOffset
      
    }
    
  }
  
  return(kit)
  
}
