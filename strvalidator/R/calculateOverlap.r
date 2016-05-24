################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 06.01.2014: Fixed factor/character bug when using frequency database.
# 30.11.2013: Specified package for function in 'plyr' -> 'plyr::rbind.fill'
# 25.09.2013: First version.

#' @title Calculate Bins Overlap
#'
#' @description
#' Analyses the bins overlap between colors.
#'
#' @details By analysing the bins overlap between dye channels a measure of
#' the risk for spectral pull-up artefacts can be obtain. The default result
#' is a matrix with the total bins overlap in number of base pairs. If an allele
#' frequency database is provided the overlap at each bin is multiplied with the
#' frequency of the corresponding allele. If no frequence exist for that allele
#' a frequency of 5/2N will be used. X and Y alleles is given the frequency 1.
#' A penalty matrix can be supplied to reduce the effect by spectral distance, 
#' meaning that overlap with the neighbouring dye can be counted in full (100%)
#' while a non neighbour dye get its overlap reduced (to e.g. 10%).
#' 
#' @param data data frame providing kit information.
#' @param db data frame allele frequency database.
#' @param penalty vector with factors for reducing the impact from distant dye channels. 
#' NB! Length must equal number of dyes in kit minus one.
#' @param virtual logical default is TRUE meaning that overlap calculation includes virtual bins.
#' @param debug logical indicating printing debug information.
#' 
#' @importFrom plyr rbind.fill
#' @importFrom utils str
#' 
#' @export
#' 
#' @return data.frame with columns 'Kit', 'Color', [dyes], 'Sum', and 'Score'. 



calculateOverlap <- function (data, db=NULL, penalty=NULL, virtual=TRUE, debug=FALSE){
  
  # BIN OVERLAP.

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data:")
    print(str(data))
    print("db:")
    print(str(db))
    print("penalty:")
    print(penalty)
    print("virtual:")
    print(virtual)
  }
  
  # Check data ----------------------------------------------------------------

  # Check db if provided.
  if(!is.null(db)){
    
    # Check if all markers are in the provided frequency database.
    chkMarkers <- unique(data$Marker)
    if(!all(chkMarkers %in% names(db))){
      message(paste("Marker(s) ", paste(chkMarkers[!chkMarkers %in% names(db)], collapse=", "),
                    " not found in frequency database!",
                    "\nA frequency of 5/2N will be used as an estimate for missing alleles/markers.",
                    sep=""))
    }
    
    # Check if 'N' (size of frequency databse) is available.
    if(is.na(db$N[1])){
      message(paste("'N' not found in frequency database!",
                    "\nA frequency of 0 will be used for missing alleles/markers.",
                    sep=""))
      
    }
  }
  
  # Prepare -------------------------------------------------------------------
  
  # Initialise variables.
  bp<-vector()
  min<-vector()
  max<-vector()
  xLower <- vector()
  xUpper <- vector()
  allele <- vector()
  virtualAllele <- vector()
  mpenalty <- matrix()
  dfRes <- data.frame()

  # Create penalty matrix.
  if(!is.null(penalty)){

    if(debug){
      print("Creating penalty matrix...")
    }
    
    # Get size (number of colors)
    msize <- length(penalty) + 1
    mrow <- list()
    
    # Loop to construct matrix rows.
    for(i in 1:msize){
      
      # Calculate indexes:
      a <- ifelse (i == 1, 0, 1)
      b <- i - 1
      c <- ifelse (i %% msize == 0, 0, 1)
      d <- msize - i
      
      # Construct a matrix row.
      mrow[[i]] <- c(rev(penalty[a:b]), NA, penalty[c:d])
      
    }
    # Construct matrix.
    mpenalty <- matrix(unlist(mrow), nrow=msize)
    
    if(debug){
      print(mpenalty)
    }
    
  }
  
  # Get kits.
  kits <- unique(data$Panel)
  kitName <- unique(data$Full.Name)
  if(debug){
    print("Kits:")
    print(kits)
    print(kitName)
  }

  # Analyse -------------------------------------------------------------------
  
  # Loop over all kits.
  for(k in seq(along=kits)){

    # Initiate variable.
    totalOverlap<-0
    
    if(debug){
      print(paste("Analysing kit", kits[k]))
    }
    
    # Kit selection.
    selectedKit <- data$Panel==kits[k]
    
    # Get colors for current kit.
    colors <- unique(data$Color[selectedKit])

    # Create result matrix. (If use matrix better to only allow one kit)
    nColor <- length(colors)
    overlapMatrix <- matrix(rep(rep(NA, nColor), nColor), nrow=nColor)
    dimnames(overlapMatrix) <- list(colors, colors)

    if(debug){
      print("Result matrix created:")
      print(overlapMatrix)
    }
    
    # Get global min/max.
    xhigh <- max(data$Size.Max[selectedKit])
    xlow <- min(data$Size.Min[selectedKit])
    yhigh <- length(colors)
    
    # Loop over all colors. This loop is the color being compared.
    for(c in seq(along=colors)){

      if(debug){
        print(paste("Analysing color", colors[c]))
      }
      
      # Get current color.
      selectedColor <- data$Color==colors[c]

      # Combine current color and current kit selection.
      selectedData <- selectedColor & selectedKit
      
      # Get lower and upper bound of bins.
      xLower <- data$Size.Min[selectedData]
      xUpper <- data$Size.Max[selectedData]
      
      # Get markers and alleles.
      # Supress:  "Warning message: NAs introduced by coercion"
      allele <- suppressWarnings(as.numeric(data$Allele[selectedData]))
      marker <- data$Marker[selectedData]
      virtualAlleles <- data$Virtual[selectedData] == 1

      # Check if using virtual alleles.
      if(!virtual){

        # Remove all elements corresponding to virtual alleles.
        xLower <- xLower[!virtualAlleles] 
        xUpper <- xUpper[!virtualAlleles] 
        allele <- allele[!virtualAlleles] 
        marker <- marker[!virtualAlleles] 
      
        if(debug){
          print("Virtual alleles removed for current colour!")
        }
        
      }
      
      # c2<-2
      # Loop over all colors. This loop is colors to compare to.
      for(c2 in seq(along=colors)){

        # Do not compare to same color.
        if(colors[c] != colors[c2]){
          
          if(debug){
            print(paste("Comparing to", colors[c2]))
          }
          
          # Get current color2.
          selectedColor2 <- data$Color==colors[c2]
          
          # Combine current color2 and current kit selection.
          selectedData2 <- selectedColor2 & selectedKit
          
          # Get lower and upper bound of bins.
          xLower2 <- data$Size.Min[selectedData2]
          xUpper2 <- data$Size.Max[selectedData2]
          virtualAlleles2 <- data$Virtual[selectedData2] == 1

          # Check if using virtual alleles.
          if(!virtual){
            
            # Remove all elements corresponding to virtual alleles.
            xLower2 <- xLower2[!virtualAlleles2] 
            xUpper2 <- xUpper2[!virtualAlleles2]
            
            if(debug){
              print("Virtual alleles removed for comparing colour!")
            }
            
          }
          
          # Get maximum size bins to compare to.
          xMax2 <- max(xUpper2)
          
          # Loop over all bins in current color.
          for(x in seq(along=xLower)){
            
            # Get current values.
            cMarker <- as.character(marker[x])
            cAllele <- as.character(allele[x])

            if(debug){
              print(paste("Checking allele ", cAllele,
                          " in marker ", cMarker,".",
                          sep=""))
            }
            
            # Exit when bin in color is > max size in color2.
            if(xUpper[x] > xMax2){
#               if(debug){
#                 print(paste("Break at", xLower2[x2],
#                             ", current bin is larger than remaining bins."))
#               }
              break
            }
            
            # Loop over bins in current color2, starting from the first overlap.
            for(x2 in which(xUpper2 > xLower[x])){
              # Compare bin in color with all bins in color2.

#               if(debug){
#                 print(paste("Enter at", xLower2[x2]))
#               }

              # Exit when bin in color2 is > current bin size in color.
              if(xLower2[x2] > xUpper[x]){
#                 if(debug){
#                   print(paste("Break at", xLower2[x2],
#                               ", remaining bins are larger than current bin."))
#                 }
                break
              }

              # Make sequence out of bins with resolution 0.01 bp.
              currentBin <- seq(xLower[x], xUpper[x], by=0.01)
              currentBin2 <- seq(xLower2[x2], xUpper2[x2], by=0.01)
              
              # Find overlapping region.
              overlap <- intersect(currentBin, currentBin2)

              # Score overlap.
              if(length(overlap)!=0){
                
                # Calculate size of overlap in base pair.
                bp <- max(overlap) - min(overlap)
                
                # Adjust score by allele frequency.
                if(!is.null(db)){

                  # X and Y is NA
                  if(is.na(cAllele)){
                    freq <- 1
                  } else {
                    # Get frequency for current allele in current marker.
                    if(cMarker %in% names(db)){
                      freq <- db[cMarker][db$Allele == cAllele, ]
                    } else {
                      # Marker not in database.
                      freq <- NA
                    }
                  }
                  
                  # Check if frequence was found.
                  if(length(freq) == 0 || is.na(freq)){

                    # Use min freq (5/2N).
                    freq <- 5 / (2 * db$N[1])
                    
                    # Check if minimum frequency could be estimated.
                    if(is.na(freq)){
                      
                      freq <- 0
                      
                      warning(paste("Marker ", cMarker,": ", cAllele,
                                    " - 'N' not found in frequency database.",
                                    " Using a frequency of ", freq, sep=""))
                              
                    } else {
                      
                      warning(paste("Marker ", cMarker,": ", cAllele,
                                    " - No frequency found in frequency database.",
                                    " Using a frequency of ", freq,
                                    " (5/2N) as an estimate!", sep=""))
                      
                    }
                    
                  } else {

                    # Adjust overlap score.
                    bp <- bp * freq
                    
                  }
                  
                }
                
                # Sum total overlap score.
                totalOverlap <- totalOverlap + bp

                if(debug){
                  print(paste(bp, "basepair overlap in", cMarker,
                              "for allele", cAllele))
                }
                
              }
              
            }
            
          }
          
          if(debug){
            print(paste("TotalOverlap", totalOverlap))
          }
          
          # Same color.
          overlapMatrix[c, c2] <- totalOverlap
          
        } else {
          
          # Same color.
          overlapMatrix[c, c2] <- NA
          
        }
        
      } # Color loop 2 ends!
      
    } # Color loop 1 ends!

    if(debug){
      print(paste("Result matrix for kit", kits[k]))
      print(overlapMatrix)
    }
    
    # Apply scoring matrix to reduce impact by dye distance.
    if(!is.null(penalty)){
      overlapMatrix <- overlapMatrix * mpenalty
      
      if(debug){
        print("penalty matrix:")
        print(mpenalty)
      }
      
    }
    
    # Create data frame.
    dfKit <- data.frame(overlapMatrix)
    names(dfKit) <- rownames(overlapMatrix)
    
    # Add color (to the left).
    dfKit <- data.frame(Color = rownames(overlapMatrix), dfKit)
    
    # Add kit name (to the left).
    dfKit <- data.frame(Kit = kitName[k], dfKit)
    
    # Calculate and add score per dye.
    dfKit$Sum <- rowSums(overlapMatrix, na.rm=TRUE)
    
    # Calculate and add total score.
    dfKit$Score <- sum(overlapMatrix, na.rm=TRUE)
    
    # Combine with previous result.
    dfRes <- plyr::rbind.fill(dfRes, dfKit)
     
   } # Kit loop ends!
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(dfRes)
  
}
