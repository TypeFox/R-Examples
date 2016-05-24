intensity_plot <-
function(x, y, strata,
                           nbin = 20,
                           xmax = max(x), xmin = min(x), ymax = max(y), ymin = min(y),
                           strata_colours = c("black", "red", "turquoise3"), verbose = TRUE,
                           xlab = "x", ylab = "y", ...) { 
  
  if(length(x) == 0L) stop("x has length 0")
  if(length(x) != length(y)) stop("x & y are not of equal length")  
  if(missing(strata)) {
    strata <- !logical(length = length(x))
    useStrata <- FALSE
  } else {
    if(!is.logical(strata)) stop("strata is not a logical vector")
    if(length(x) != length(strata)) stop("x & strata are not of equal length")
    useStrata <- TRUE
  }
  
  removeNA <- is.na(x) | is.na(y)
  if(any(removeNA)) {
    if(all(removeNA)) stop("No non-missing pairs in vectors x and y")
    if(verbose) print(paste(" - excluding", sum(removeNA), "subjects with missing values"), quote = FALSE)
    x <- x[!removeNA]
    y <- y[!removeNA]
    strata <- strata[!removeNA]
  }
  removeNA <- sum(removeNA)
  
  if(min(x) < xmin) removeMinMax <- x < xmin else removeMinMax <- logical(length = length(x))
  if(max(x) > xmax) removeMinMax <- removeMinMax | x > xmax
  if(min(y) < ymin) removeMinMax <- removeMinMax | y < ymin
  if(max(y) > ymax) removeMinMax <- removeMinMax | y > ymax
  
  if(any(removeMinMax)) {
    if(all(removeMinMax)) stop("All samples are outside of predefined x/y range")
    if(verbose) print(paste(" - excluding", sum(removeMinMax), "Subjects outside of predefined x/y range"), quote = FALSE)
    x <- x[!removeMinMax]
    y <- y[!removeMinMax]
    strata <- strata[!removeMinMax]
  }
  removeMinMax <- sum(removeMinMax)
  
  if(any(is.na(strata))) stop("Missing values in strata vector")
  
  if(useStrata) {
    useHQ <- any( strata)
    useLQ <- any(!strata)
  } else {
    useHQ <- TRUE
    useLQ <- FALSE
  }
  
  xbin_width <- (xmax - xmin) / nbin
  xbin <- ifelse(x == xmax,
                 (nbin - 0.5) * xbin_width + xmin,
                 (floor((x - xmin)/ xbin_width) + 0.5) * xbin_width + xmin)
  xbin_values <- xmin + (1:nbin - 0.5) * xbin_width
  
  ybin_width <- (ymax - ymin) / nbin
  ybin <- ifelse(y == ymax,
                 (nbin - 0.5) * ybin_width + ymin,
                 (floor((y - ymin)/ ybin_width) + 0.5) * ybin_width + ymin)
  ybin_values <- ymin + (1:nbin - 0.5) * ybin_width
  
  table_counts_all <- table(xbin, ybin)
  table_counts <- as.list(rep(NA,2))
  table_freqs <- as.list(rep(NA,2))
  
  if(useHQ) {
    table_counts[[1]] <- table(c(rep(xbin_values, nbin), xbin[ strata]), c(sort(rep(ybin_values, nbin)), ybin[ strata])) - 1
    table_freqs[[1]] <- table_counts[[1]]/sum(table_counts_all)
    cexHQ <- table_freqs[[1]]^0.25*8
  }
  if(useLQ) {
    table_counts[[2]] <- table(c(rep(xbin_values, nbin), xbin[!strata]), c(sort(rep(ybin_values, nbin)), ybin[!strata])) - 1
    table_freqs[[2]] <- table_counts[[2]]/sum(table_counts_all)
    cexLQ <- table_freqs[[2]]^0.25*8
  }
  
  plot(x = c(xmin, xmin, xmax, xmax), y = c(ymin, ymax, ymin, ymax),
       col = "transparent", xlab = xlab, ylab = ylab, ...)
  if(useHQ & useLQ) {
    plotX <- rep(xbin_values, nbin)
    plotY <- sort(rep(ybin_values, nbin))
    
    strata_overlap <- abs(cexLQ - cexHQ) < 0.8 & cexHQ != 0 & cexLQ != 0
    # strata overlap selects  bins that are very similar in size.
    # Obviously, if either of the bins is empty, there is no overlap
    
    points(plotX, plotY,
           cex = ifelse(cexHQ - cexLQ >= 0.8 | cexLQ == 0, cexHQ, NA),
           pch = 19, col = strata_colours[1])
    points(plotX, plotY,
           cex = ifelse(strata_overlap,
                        ifelse(cexLQ > cexHQ, cexLQ, cexHQ),
                        cexLQ),
           pch = 19,
           col = ifelse(strata_overlap, strata_colours[3], strata_colours[2]))
    points(plotX, plotY,
           cex = ifelse(cexLQ - cexHQ >= 0.8 | cexLQ == 0, cexHQ, NA),
           pch = 19, col = strata_colours[1])
  } else {
    points(rep(xbin_values, nbin), sort(rep(ybin_values, nbin)),
           cex = if(useHQ) cexHQ else cexLQ,
           pch = 19, col = strata_colours[if(useHQ) 1 else 2])
  }
  return(invisible(list(NA_removed = removeNA, outliers_removed = removeMinMax)))  
}
