################################################################################
# TODO LIST
# TODO: Fix problem with overplotting resulting in 'invisible' peaks.

################################################################################
# CHANGE LOG
# 31.12.2015: Rewritten function.
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 31.05.2015: Added 'numbered=TRUE' to 'slim' function.
# 10.02.2015: Changed error message.
# 09.12.2014: Function moved from PCRsim package.

#' @title Generate EPG
#'
#' @description
#' Visualises an EPG from DNA profiling data.
#'
#' @details
#' Generates a electropherogram like plot from 'data' and 'kit'.
#' If 'Size' is not present it is estimated from kit information and allele values.
#' If 'Height' is not present a default of 1000 RFU is used.
#' Off-ladder alleles can be plotted if 'Size' is provided.
#' There are various options to customise the plot scale and labels.
#' It is also possible to plot 'distributions' of peak heights as boxplots.
#' 
#' @param data data frame containing at least columns 'Sample.Name', 'Allele', and 'Marker'.
#' @param kit string or integer representing the STR typing kit.
#' @param title string providing the title for the EPG.
#' @param wrap logical TRUE to wrap by dye.
#' @param peaks logical TRUE to plot peaks for distributions using mean peak height.
#' @param boxplot logical TRUE to plot distributions of peak heights as boxplots.
#' @param collapse logical TRUE to add the peak heights of identical alleles peaks within each marker.
#' NB! Removes off-ladder alleles.
#' @param silent logical FALSE to show plot.
#' @param ignore.case logical FALSE for case sensitive marker names.
#' @param at numeric analytical threshold (Height <= at will not be plotted).
#' @param scale character "free" free x and y scale, alternatively "free_y" or "free_x".
#' @param limit.x logical TRUE to fix x-axis to size range.
#' To get a common x scale set scale="free_y" and limit.x=TRUE.
#' @param label.size numeric for allele label text size.
#' @param label.angle numeric for allele label print angle.
#' @param label.vjust numeric for vertical justification of allele labels.
#' @param label.hjust numeric for horizontal justification of allele labels.
#' @param expand numeric for plot are expansion (to avoid clipping of labels).
#' @param debug logical for printing debug information to the console.
#' 
#' @return ggplot object.
#' 
#' @export
#' 
#' @importFrom utils str head tail
#' @importFrom stats as.formula
#' @importFrom ggplot2 geom_polygon ggplot aes_string scale_fill_manual
#' geom_boxplot scale_colour_manual geom_rect geom_text scale_y_continuous
#' facet_grid facet_wrap coord_cartesian theme element_blank labs xlab ylab
#' 

generateEPG <- function(data, kit, title = NULL, wrap = TRUE, boxplot = FALSE,
                        peaks = TRUE, collapse = TRUE, silent = FALSE,
                        ignore.case = TRUE, at = 0, scale = "free",
                        limit.x = TRUE, label.size = 3, label.angle = 0,
                        label.vjust = 1, label.hjust = 0.5, expand = 0.1,
                        debug = FALSE){

  # Debug info.
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data:")
    print(str(data))
    print(head(data))
    print(tail(data))
  }
  
  if(!"Sample.Name" %in% names(data)){
    stop("'data' must contain a column 'Sample.Name'", call. = TRUE)
  }
  
  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'", call. = TRUE)
  }
  
  if(length(grep("Allele", names(data))) == 0){
    stop("'data' must contain a column 'Allele'.", call. = TRUE)
  }

  if(!is.logical(wrap)){
    stop("'wrap' must be logical.", call. = TRUE)
  }
  
  if(!is.logical(boxplot)){
    stop("'boxplot' must be logical.", call. = TRUE)
  }
  
  if(!is.logical(peaks)){
    stop("'peaks' must be logical.", call. = TRUE)
  }
  
  if(!is.logical(collapse)){
    stop("'collapse' must be logical.", call. = TRUE)
  }
  
  if(!is.logical(silent)){
    stop("'silent' must be logical.", call. = TRUE)
  }
  
  if(!is.logical(ignore.case)){
    stop("'ignore.case' must be logical.", call. = TRUE)
  }
  
  if(!is.numeric(at)){
    stop("'at' must be numeric.", call. = TRUE)
  }
  
  if(!is.numeric(label.size)){
    stop("'label.size' must be numeric.", call. = TRUE)
  }
  
  if(!is.numeric(label.angle)){
    stop("'label.angle' must be numeric.", call. = TRUE)
  }
  
  if(!is.numeric(label.vjust)){
    stop("'label.vjust' must be numeric.", call. = TRUE)
  }
  
  if(!is.numeric(label.hjust)){
    stop("'label.hjust' must be numeric.", call. = TRUE)
  }
  
  if(!is.numeric(expand)){
    stop("'expand' must be numeric.", call. = TRUE)
  }
  
  # Prepare -------------------------------------------------------------------

  # Width of peaks in base pair.
  width = 1
  
  if(!collapse){
    
    if(boxplot){
      boxplot <- FALSE
      message("boxplot set to FALSE since collapse=FALSE")
    }
    
    if(peaks){
      peaks <- FALSE
      message("peaks set to FALSE since collapse=FALSE")
    }
    
  }
  
  # Add missing columns .......................................................

  # Check if height column exist.  
  if(!"Height" %in% names(data)){
    
    # Add Height if not present. 
    data$Height <- 1000
    
    message("'Height' is missing. Using default!")
    
  }

  # Check NA's.  
  if(any(is.na(data$Height))){
    
    tmp1 <- nrow(data)
    
    # Remove rows with zero height.
    data <- data[!is.na(data$Height), ]
    
    tmp2 <- nrow(data)
    
    message(tmp1 - tmp2, " peaks with Height = NA removed from data")
    
  }
  
  # Check height.  
  if(any(data$Height == 0)){
    
    tmp1 <- nrow(data)
    
    # Remove rows with zero height.
    data <- data[data$Height != 0, ]
    
    tmp2 <- nrow(data)
    
    message(tmp1 - tmp2, " peaks with Height = 0 removed from data")
    
  }

  # Check if dye column exist.
  if(!"Dye" %in% names(data)){
    
    message("'Dye' information not in 'data'.")
    
    # Add dye information.
    data <- addColor(data = data, kit = kit)
    
    message("Added 'Dye' information.")
    
  }

  # Check format ..............................................................
  
  if(!is.numeric(data$Height)){
    
    data$Height <- as.numeric(data$Height)
    
  }
  
  # Check if 'fat' format.
  if(length(grep("Allele", names(data))) > 1) {
    
    message("'fat' data format detected.")
    
    fixCol <- colNames(data = data, slim = TRUE, numbered = TRUE,
                       concatenate = NULL, debug = debug)
    
    stackCol <- colNames(data = data, slim = FALSE, numbered = TRUE,
                         concatenate = NULL, debug = debug)
    
    # Slim data frame.
    data <- slim(data = data, fix = fixCol, stack = stackCol, debug = debug)

    message("data converted to 'slim' format.")
    
  }

  # Filter data ...............................................................
  
  # Apply analytical threshold (AT).
  if(any(data$Height < at)){
    tmp1 <- nrow(data)
    data <- data[!(data$Height < at), ]
    tmp2 <- nrow(data)
    message("Removed", tmp1 - tmp2, "peaks below at =", at, "RFU")
  }
  

  # Get kit information .......................................................

  # Get kit markers, ranges, and colors.
  kitInfo <- getKit(kit=kit, what="Range")
  kitInfo <- addColor(kitInfo, have="Color", need="Dye")

  # Get unique Dyes.
  kitDye <- unique(kitInfo$Dye)
  
  # Get unique Colors.
  kitColors <- unique(kitInfo$Color)
  
  # Convert R colors.
  manualPlotColors <- addColor(kitColors, have="Color", need="R.Color")

    
  # Create EPG ----------------------------------------------------------------

  # Create id .................................................................
  
  # Add unique 'Id' column for grouping.
  if("Size" %in% names(data)){
    
    # Check if numeric.
    if(!is.numeric(data$Size)){
      
      # Convert to numeric.    
      data$Size <- as.numeric(data$Size)
      
      message("'Size' must be numeric. 'data' converted!")
      
    }  
    
    # Combine rounded 'Size' and 'Marker'.
    # This can preserve 'OL' peaks but may also result in two peaks for alleles.
    data$Id <- paste(round(data$Size, 0), data$Marker, sep="")
    
  } else {
    
    # Combine 'Allele' and 'Marker'.
    # This will add all 'OL' in one marker even if originally at different size.
    data$Id <- paste(data$Allele, data$Marker, sep="") 
    
  }
  
  # Collapse ..................................................................
  
  # 'Collapse' will add peak heights of identical alleles.
  # Not collapsing will 'overplot' samples on top of each other
  # without adding peak heights.
  if(collapse){

    message("Collapse dataset to mean peak heights over multiple samples.")

    # Convert to data.table for performance.
    DT <- data.table::data.table(data)
    
    # Calculate sum of peak heights for identical alleles in each sample.
    # to be used in boxplot.
    DT <- DT[, list(Marker = Marker, Dye = Dye, Allele = Allele,
                    Height = sum(Height)),
             by = list(Sample.Name, Id)]
    
    # Calculate mean peak height for each allele across all samples.
    # If plot peaks are true for boxplot.
    dataMean <- DT[, list(Sample.Name = "Mean", Height = mean(Height)),
                   by = list(Marker, Dye, Allele, Id)]

    # Add size.
    dataMean <- addSize(data = dataMean, kit = getKit(kit = kit, what = "Offset"),
                        bins = FALSE, ignore.case = ignore.case, debug = debug)

    # Check if distribution (boxplot).
    if(!boxplot){
        
      message("Collapse dataset by adding peak heights for all profiles.")
      
      # Calculate sum of peak heights for identical alleles across all samples.        
      DT <- DT[, list(Sample.Name = "Profile", Height = sum(Height)),
               by = list(Marker, Dye, Allele, Id)]
      
    } # end distribution.
    
    # Convert to data.frame to make sure strvalidator functions are working.
    data <- data.frame(DT)
      
  } # end collapse.
  
  # Check if size column exist.
  if(!"Size" %in% names(data)){
    
    message("'Size' information not in 'data'.")
    
    # Add estimated size.
    data <- addSize(data = data, kit = getKit(kit = kit, what = "Offset"),
                    bins = FALSE, ignore.case = ignore.case)
    
    message("Added estimated 'Size' information.")

    # Remove NA rows.
    if(any(is.na(data$Size))){
      tmp1 <- nrow(data)
      data <- data[!is.na(data$Size),]
      tmp2 <- nrow(data)
      message(tmp1 - tmp2, " rows with Size=NA removed.")
    }

  }
  
  # Check if numeric.
  if(!is.numeric(data$Size)){
    
    # Convert to numeric.    
    data$Size <- as.numeric(data$Size)
    
    message("'Size' must be numeric. 'data' converted!")
    
  }  

  # Sort 'Marker' and 'Dye' factors according 'kit'.
  data <- sortMarker(data = data, kit = kit)
  if(peaks){
    dataMean <- sortMarker(data = dataMean, kit = kit)
  }
  
  # Calculate coordinates for plotting peaks ..................................
  
  # Create new dataframe.  
  dataPeaks <- heightToPeak(data = data, width = width, debug = debug)
  if(peaks){
    dataMean <- heightToPeak(data = dataMean, width = width, debug = debug)
  }
  
  # Create allele labels ......................................................
  
  # Copy unique 'Marker'-'Allele' combinations in 'data'
  # to a new data frame for handling allele names.
  alleleInfo <- data[data$Height != 0,]
  # Remove duplicates.
  alleleInfo <- alleleInfo[!duplicated(alleleInfo[c("Marker", "Allele","Size")]), ]

  # Check if NA's in Size.
  if(any(is.na(alleleInfo$Size))){
    
    tmp1 <- sum(is.na(alleleInfo$Size))
    
    message("NA's in allele info Size")

    # Replace NA with the smallest size in kit (plot can't handle all NAs).
    alleleInfo$Size[is.na(alleleInfo$Size)] <- min(kitInfo$Marker.Min)

    tmp2 <- sum(is.na(alleleInfo$Size))
    
    message(tmp1 - tmp2, "NA's replaced by", min(kitInfo$Marker.Min))
    
  }
  
  # Create marker ranges ......................................................
  
  # Get information for annotation of markers.
  mDye <- kitInfo$Dye
  mXmin <- kitInfo$Marker.Min
  mXmax <- kitInfo$Marker.Max
  mText <- kitInfo$Marker
  mYmax <- vector()
  
  # Find max peak height by dye channel.
  DT <- data.table::data.table(data)
  
  # NB! Use keyby instead of key to sort the result.
  tmpYmax <- DT[, list(Max = max(Height)), keyby = Dye]
  mYtimes <- as.vector(table(kitInfo$Dye))

  # Check scale.  
  if(scale == "free_x"){
    
    # This means y max is equal for all dyes.
    mYmax <- rep(max(tmpYmax$Max), sum(mYtimes))
    
  } else{
    
    # Different y max for each dye.    
    mYmax <- rep(tmpYmax$Max, mYtimes)
    
  }
  
  # Create annotation data frame for loci.
  markerRanges <- data.frame(Dye = factor(mDye, levels = unique(mDye)), # Facet.
                             Color = mDye,    	         # Dye.
                             Xmin = mXmin,			         # Marker lower range.	
                             Xmax = mXmax,			         # Marker upper range.
                             Size = (mXmin + mXmax) / 2, # Midpoint of marker range.
                             Height = mYmax,		         # Lower edge of marker range.
                             Top = mYmax * 1.1,          # Upper edge of marker range.
                             Text = mText)			         # Marker names.

  # Create plot ...............................................................

  # Create plot.
  #gp <- ggplot(data = data, aes_string(x = "Size", y = "Height"))
  gp <- ggplot(data = dataPeaks, aes_string(x = "Size", y = "Height"))
  
  # Plot data.
  if(!boxplot){
    
    # Plot height as peaks.
    gp <- gp + geom_polygon(aes_string(group = "Id", fill = "Dye"), data = dataPeaks)

  } else {
    
    # Plot boxplots for distributions.
    gp <- gp + geom_boxplot(aes_string(group = "Id", color = "Dye"),
                            outlier.size = 1, data = data)

    if(peaks){
      
      # Plot mean peak height as peaks.
      gp <- gp + geom_polygon(aes_string(group = "Id", fill = "Dye"),
                              data = dataMean)
      
    }
    
  }

  # Add colours.
  gp <- gp + scale_fill_manual(values = manualPlotColors)

  if(wrap){
    # Add marker regions, names, and wrap by colour.
    
    # Add marker regions.
    gp <- gp + geom_rect(aes_string(xmin = "Xmin", xmax = "Xmax",
                                    ymin = "Height", ymax = "Top"),
                         alpha = .2, data = markerRanges,
                         fill = "blue", color = "red")
    
    # Add marker names.
    gp <- gp + geom_text(aes_string(label = "Text", y = "Top"), 
                         data = markerRanges, size = 3, vjust = 1)
    
    # Add allele names.
    gp <- gp + geom_text(aes_string(label = "Allele", x = "Size", y = 0),
                         data = alleleInfo, size = label.size,
                         angle = label.angle, vjust = label.vjust,
                         hjust = label.hjust)
    
    # expand plot area (to avoid clipping).
    gp <- gp + scale_y_continuous(expand=c(expand, 0))

    # NB! 'facet_wrap' does not seem to support strings.
    #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
    gp <- gp + facet_wrap(as.formula(paste("~", "Dye")), ncol = 1,
                          drop = FALSE, scales = scale)
    
  }
  
  # Set limits.
  if(limit.x){
    # Set x-axis limits to total marker range.
    gp <- gp + coord_cartesian(xlim = c(min(mXmin), max(mXmax))) 
  }
  
  # Strip facet labels and background.
  gp <- gp + theme(strip.text = element_blank())
  gp <- gp + theme(strip.background = element_blank())
  
  # Add title and axis labels.
  gp <- gp + labs(title=title)
  gp <- gp + xlab("Size (bp)")
  gp <- gp + ylab("Peak height (RFU)")
  
  # Show plot.
  if(!silent){
    print(gp)
  }
  
  # Debug info.
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return plot object.
  return(gp)
  
}
