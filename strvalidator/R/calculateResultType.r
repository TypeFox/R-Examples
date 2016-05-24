################################################################################
# TODO LIST
# TODO: use string constants instead of hard coded.

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 15.12.2014: Changed parameter names to format: lower.case
# 22.01.2014: Fixed bug by adding check that 'Height' is numeric and convert.
# 15.01.2014: Fixed NA's when 'mixture.limits' and 'partial.limits' is NULL.
# 15.01.2014: Added message to show progress.
# 03.11.2013: Added debug parameter and data check.
# <03.11.2013: Roxygenized.
# <03.11.2013: Added factors.
# <03.11.2013: Grouping of mixed results.
# <03.11.2013: Grouping of partial results.
# <03.11.2013: First version.

#' @title Calculate Result Type
#'
#' @description
#' Calculate the result type for samples.
#'
#' @details
#' Calculates result types for samples in 'data'.
#' Defined types are: 'No result', 'Mixture', 'Partial', and 'Complete'.
#' Subtypes can be defined by parameters.
#' An integer passed to 'threshold' defines a subtype of 'Complete' "Complete profile all peaks >threshold".
#' An integer or vector passed to 'mixture.limits' define subtypes of 'Mixture' "> [mixture.limits] markers".
#' An integer or vector passed to 'partial.limits' define subtypes of 'Partial' "> [partial.limits] peaks".
#' A string with marker names separated by pipe (|) passed to 'marker.subset' and
#'  a string 'subset.name' defines a subtype of 'Partial' "Complete [subset.name]".
#'  
#' @param data a data frame containing at least the column 'Sample.Name'.
#' @param kit character string or integer defining the kit.
#' @param add.missing.marker logical, defualt is TRUE which adds missing markers.
#' @param threshold integer indicating the dropout threshold.
#' @param mixture.limits integer or vector indicating subtypes of 'Mixture'.
#' @param partial.limits integer or vector indicating subtypes of 'Partial'.
#' @param subset.name string naming the subset of 'Complete'.
#' @param marker.subset string with marker names defining the subset of 'Complete'.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with columns 'Sample.Name','Type', and 'Subtype'.
#' 
#' @export
#' 
#' @importFrom utils head
#' 

calculateResultType <- function(data, kit=NULL, add.missing.marker=TRUE,
                                threshold=NULL, mixture.limits=NULL,
                                partial.limits=NULL, subset.name=NA,
                                marker.subset=NULL, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("head(data)")
    print(head(data))
    print("threshold")
    print(threshold)
    print("mixture.limits")
    print(mixture.limits)
    print("partial.limits")
    print(partial.limits)
    print("subset.name")
    print(subset.name)
    print("marker.subset")
    print(marker.subset)
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset.
  if(!any(grepl("Sample.Name", names(data)))){
    stop("'data' must contain a column 'Sample.Name'",
         call. = TRUE)
  }
  
  if(!any(grepl("Marker", names(data)))){
    stop("'data' must contain a column 'Marker'",
         call. = TRUE)
  }
  if(!any(grepl("Allele", names(data)))){
    stop("'data' must contain a column 'Allele'",
         call. = TRUE)
  }
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(add.missing.marker){
    if(is.null(kit)){
      stop("'kit' must be provided if 'add.missing.marker' is TRUE",
           call. = TRUE)
    } else {
      if(is.na(getKit(kit=kit, what="Short.Name"))){
        stop(paste("'kit' does not exist", "\nAvailable kits:",
                   paste(getKit(), collapse=", ")), call. = TRUE)
      }
    }
  }

  # PREPARE -------------------------------------------------------------------
  
  if(add.missing.marker){
    # Add missing markers to samples.
    markers <- getKit(kit=kit, what="Marker")
    data <- addMarker(data=data, marker=markers, ignore.case=TRUE, debug=debug)
  }
  
  if(!is.numeric(data$Height)){
    message("'Height' not numeric. Converting to numeric.")
    data$Height <- as.numeric(data$Height)
  }
  
  # CALCULATE -----------------------------------------------------------------
  # NB! Strings used for classification must be identical to the ones used to
  # create factors.
  
  # Get sample names.
	sampleNames <- unique(data$Sample.Name)

	# Create result data frame.
	res <- data.frame(matrix(NA, length(sampleNames), 3))
	# Add new column names.
	names(res) <- paste(c("Sample.Name","Type","Subtype"))  

	# Loop over all samples.
	for(s in seq(along=sampleNames)){

	  # Show progress.
	  message(paste("Calculate result type for sample (",
                  s, " of ", length(sampleNames), "): ", sampleNames[s], sep=""))
	  
		# Get current sample.
		sampleData <- data[data$Sample==sampleNames[s],]
    
    if(debug){
		  print("Current sample data:")
		  print(sampleData)
		}
		
		# Check result type.
		if(all(is.na(sampleData$Allele))){
			# No result.

			res[s, ] <- c(sampleNames[s], "No result", "No result")
		
		} else if(max(table(sampleData$Marker))>2){
			# Mixture.

			markers <- length(unique(sampleData$Marker[!is.na(sampleData$Allele)]))
			if(!is.null(mixture.limits)){
				for(t in rev(seq(along=mixture.limits))){
					if(markers<=mixture.limits[t]){
						subtype <- paste("<=", mixture.limits[t], "markers")
					} else if(markers > mixture.limits[length(mixture.limits)]){
						subtype <- paste(">", mixture.limits[length(mixture.limits)], "markers")
					}
				}
			} else {
					subtype <- paste("Mixture")
			}
			res[s, ] <- c(sampleNames[s], "Mixture", subtype)
			
		} else if(any(is.na(sampleData$Allele))){
			# Partial profile.

			alleles <- sum(!is.na(sampleData$Allele))
			if(!is.null(partial.limits)){
				for(t in rev(seq(along=partial.limits))){
					if(alleles<=partial.limits[t]){
						subtype <- paste("<=", partial.limits[t], "peaks")
					} else if(alleles > partial.limits[length(partial.limits)]){
						subtype <- paste(">", partial.limits[length(partial.limits)], "peaks")
					} 
				}
			} else {
					subtype <- paste("Partial")
			}
			res[s, ] <- c(sampleNames[s], "Partial", subtype)

			# Check for subset.
			if(!is.null(marker.subset)){
				# Subset data.
				selectedMarkers <- grepl(marker.subset, sampleData$Marker)		
				if(all(!is.na(sampleData$Allele[selectedMarkers]))){
					# Full subset profile.
					res[s, ] <- c(sampleNames[s], "Partial" , paste("Complete", subset.name))
				}
			}
			
		} else if(!any(is.na(sampleData$Allele))){
			# Complete profile.
			res[s, ] <- c(sampleNames[s], "Complete profile", "Complete profile")

			# Check against threshold.
			if(!is.null(threshold) && all(sampleData$Height>threshold)){
				# Complete profile, all peaks > T.
				res[s, ] <- c(sampleNames[s], "Complete profile", paste("all peaks >", threshold))
			}

		}
	}

  # FACTORS -------------------------------------------------------------------
  
	# Construct factor levels in correct order.
  # NB! Strings must be identical to the ones used in classification.
	factorLabels <- NULL
	blankLabels <- NULL
	mixtureLabels <- NULL
	partialLabels <- NULL
	completeLabels <- NULL

	factorLabelsSub <- NULL
	blankLabelsSub <- NULL
	mixtureLabelsSub <- NULL
	partialLabelsSub <- NULL
	completeLabelsSub <- NULL
	
	# Partial Labels.

	if(!is.null(marker.subset)){
		partialLabelsSub <- c(partialLabelsSub, paste("Complete", subset.name))
	}
	if(!is.null(partial.limits)){
		for(t in rev(seq(along=partial.limits))){
			if(t == length(partial.limits)){
				partialLabelsSub <- c(partialLabelsSub, paste(">", partial.limits[t], "peaks"))
			}
			partialLabelsSub <- c(partialLabelsSub, paste("<=", partial.limits[t], "peaks"))
		}
	}
	partialLabels <- "Partial"
	partialLabelsSub <- c("Partial", partialLabelsSub)

	# Mixture Labels.
	if(!is.null(mixture.limits)){
		for(t in rev(seq(along=mixture.limits))){
			if(t == length(mixture.limits)){
				mixtureLabelsSub <- c(mixtureLabelsSub, paste(">", mixture.limits[t], "markers"))
			}
			mixtureLabelsSub <- c(mixtureLabelsSub, paste("<=", mixture.limits[t], "markers"))
		}
	}
	mixtureLabels <- "Mixture"
	mixtureLabelsSub <- c("Mixture", mixtureLabelsSub)

	# Complete Labels.
	if(!is.null(threshold)){
		completeLabelsSub <- c(completeLabelsSub, paste("all peaks >", threshold))
	}
	completeLabels <- "Complete profile"
	completeLabelsSub <- c(completeLabelsSub, "Complete profile")

	# Blank Labels.
	blankLabels <- "No result"
	blankLabelsSub <- "No result"

	# All factor labels.
	factorLabels <- c(mixtureLabels, completeLabels, partialLabels, blankLabels)
	factorLabelsSub <- c(mixtureLabelsSub, completeLabelsSub, partialLabelsSub, blankLabelsSub)

  if(debug){
    print("factorLabels")
    print(factorLabels)
    print("factorLabelsSub")
    print(factorLabelsSub)
  }
  
	# Assign factors.
	res$Type <- factor(res$Type, levels = factorLabels)
	res$Subtype <- factor(res$Subtype, levels = factorLabelsSub)

  if(debug){
    print("head(res):")
    print(head(res))
    print(paste("EXIT:", match.call()[[1]]))
  }
  
	return(res)
}

