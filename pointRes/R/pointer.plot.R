#' Plot pointer years for multiple sites
#'
#' @description The function creates a dot plot showing positive and (or) negative pointer years from \code{lists} of the type as produced by either \code{\link{pointer.norm}} or \code{\link{pointer.rgc}}.
#' 
#' @usage pointer.plot(list.sites, sign = c("both", "pos", "neg"),
#'              start.yr = NULL, end.yr = NULL, labels = NULL,
#'              x.tick.major = 10, x.tick.minor = 5) 
#'
#' @param list.sites a \code{list} with \code{lists} as produced by either \code{\link{pointer.norm}} or \code{\link{pointer.rgc}} for individual sites (created using list(site1, site2,..)).
#' @param sign a \code{character} string specifying whether both positive and negative (\code{"both"}), or only positive (\code{"pos"}) or negative (\code{"neg"}) pointer years should be displayed. Defaults to \code{"both"}.
#' @param start.yr an \code{integer} specifying the first year to be plotted. Defaults to the first year with data if \code{\var{start.yr}} is \code{NULL}.
#' @param end.yr an \code{integer} specifying the last year to be plotted. Defaults to the last year with data if \code{\var{end.yr}} is \code{NULL}.
#' @param labels a \code{character} vector with labels for the sites. Defaults to 'site 1, 2, .., \code{\var{i}}'.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a dot plot showing pointer years for multiple sites. Positive and negative pointer years are indicated with different symbols. If event years were defined using \code{method.thresh "Neuwirth"} (\code{\link{pointer.norm}}), different tones of gray indicate weak, strong and extreme pointer years, based on the most common event year class.
#' 
#' @return 
#' Dot plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot negative pointer years for multiple sites from pointer.rgc output
#' data(s033)
#' site1 <- pointer.rgc(s033, nb.yrs = 4)
#' site2 <- pointer.rgc(s033, nb.yrs = 6)
#' sites <- list(site1, site2)
#' pointer.plot(sites, sign = "neg", start.yr = 1950, end.yr = NULL) 
#'
#' ## Plot pointer years for multiple sites from pointer.norm output (method "Neuwirth")
#' data(s033)
#' site1 <- pointer.norm(s033, window = 5, method.thresh = "Neuwirth")
#' site2 <- pointer.norm(s033, window = 11, method.thresh = "Neuwirth")
#' sites <- list(site1, site2)
#' site.names <- c("schneetal5", "schneetal11")
#' pointer.plot(sites, start.yr = 1950, end.yr = NULL, labels = site.names) 
#'            
#' @import ggplot2
#' @import stats
#' @importFrom plyr round_any
#' @importFrom TripleR matrix2long
#'       
#' @export pointer.plot
#' 
pointer.plot <- function(list.sites, sign = c("both", "pos", "neg"), start.yr = NULL, end.yr = NULL, labels = NULL, x.tick.major = 10, x.tick.minor = 5) 
{
  stopifnot(is.list(list.sites))
  for(i in 1:length(list.sites)) {
    if(FALSE %in% (class(list.sites[[i]])[1] == "pointer.rgc") & FALSE %in% (class(list.sites[[i]])[1] == "pointer.norm")){
      stop("'list.sites' contains list(s) that are no output of function pointer.rgc or pointer.norm")
    }
  }
  check2 <- vector()
  for(i in 1:length(list.sites)) {
    check2[i] <- class(list.sites[[i]])[1] == "pointer.norm"
  }
  if((TRUE %in% check2) && (FALSE %in% check2)){
    stop("'list.sites' contains output of both function pointer.rgc and pointer.norm")
  }
  if(check2[1] == TRUE) {
    check3 <- vector()
    for(i in 1:length(list.sites)){
      check3[i] <- class(list.sites[[i]])[2] == "Neuwirth"
    }
    if((TRUE %in% check3) && (FALSE %in% check3)) {
      stop("'list.sites' contains output of function pointer.norm for both method Cropper and Neuwirth")
    }
  } else {
      check3 <- vector()
      check3[1] <- FALSE
  }
  if (x.tick.minor > x.tick.major){
    stop("'x.tick.minor' should be smaller then 'x.tick.major'")
  }
  
  sign2 <- match.arg(sign, c("both", "pos", "neg"))
  
  year <- site <- PYvalues <- NULL
  
  if(length(labels) == 0) {
    labels2 <- vector()
    for(i in 1:length(list.sites)) {
      labels2[i] <- paste("site", i, sep = " ")
    }
  } else {
      labels2 <- labels
  }
  
  vec.min <- vec.max <- vector()
  for(i in 1:length(list.sites)) {
    vec.min[i] <- min(list.sites[[i]]$out[, "year"])
    vec.max[i] <- max(list.sites[[i]]$out[, "year"])
  }
  min.yrs <- min(vec.min)
  max.yrs <- max(vec.max)
  yrs <- seq(min.yrs, max.yrs, 1)
  
  if(!is.null(start.yr) && start.yr < min.yrs) {
    stop("'start.yr' is out of bounds. By default (start.yr = NULL) the first year is displayed")
  }
  if(!is.null(end.yr) && end.yr > max.yrs) {
    stop("'end.yr' is out of bounds. By default (end.yr = NULL) the last year is displayed")
  }
  
  nature <- matrix(nrow = length(yrs), ncol = length(list.sites))
  rownames(nature) <- yrs
  for(i in 1:length(list.sites)) {
    nature[as.character(list.sites[[i]]$out[, "year"]),i] <- list.sites[[i]]$out[, "nature"]
  }
  colnames(nature) <- labels2
  
  start.yr2 <- ifelse(length(start.yr) != 0, start.yr, min.yrs)
  end.yr2 <- ifelse(length(end.yr) != 0, end.yr, max.yrs)
  start.yr3 <- round_any(start.yr2, 10, f = floor)
  end.yr3 <- round_any(end.yr2, 5, f = ceiling)
  
  if(check3[1] == TRUE) { 
    int.class <- matrix(nrow = length(yrs), ncol = length(list.sites))
    rownames(int.class) <- yrs
    for(i in 1:length(list.sites)) {
      int.class[as.character(list.sites[[i]]$out[, "year"]),i] <- ifelse(list.sites[[i]]$out[, "nature"] == (-1), 
             max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,6:8], ties.method = "first"), 
             max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,3:5], ties.method = "first"))
    }
    colnames(int.class) <- labels2
    
    input <- matrix2long(t(nature), new.ids = FALSE)
    input.int <- matrix2long(t(int.class), new.ids = FALSE)
    input[,4] <- input.int[,3]
    input[,4] <- ifelse(input[,3] == (-1), paste("-", input[,4], sep = ''), input[,4])
    input[,4] <- ifelse(input[,3] == 0, 0, input[, 4])
    input2 <- na.omit(input)
    rownames(input2) <- NULL
    colnames(input2) <-c ("site", "year", "PYvalues", "int.class")
    input3 <- subset(input2, input2[, "year"] >= start.yr2 & input2[, "year"] <= end.yr2)
    rownames(input3) <- NULL
    
    if(sign2 == "both") {
      int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
      fill.levels <- c("black", "#bdbdbd", "white" , "#bdbdbd", "white", "#bdbdbd", "black")
      label.levels <- c("negative extreme", "negative strong", "negative weak", "none", 
                        "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(25, 25, 25, 95, 24, 24, 24)
    }
    if(sign2 == "pos") {
      input3[input3$int.class < 0, "int.class"] <- 0
      int.levels <- c(0, 1, 2, 3)
      fill.levels <- c("#bdbdbd", "white", "#bdbdbd", "black")
      label.levels <- c("other", "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(95, 24, 24, 24)
    }
    if(sign2 == "neg") {
      input3[input3$int.class > 0, "int.class"] <- 0
      int.levels <- c(-3, -2, -1, 0)
      fill.levels <- c("black", "#bdbdbd", "white" , "#bdbdbd")
      label.levels <- c("negative extreme", "negative strong", "negative weak", "other")
      shape.levels <- c(25, 25, 25, 95)
    }
    
    pl <- ggplot(input3, aes(x = year, y = site, shape = factor(int.class), 
                             fill = factor(int.class)))
    pl + geom_point(size = 2, colour = "black") +
      scale_shape_manual(name = "pointer year class", limits = int.levels,
                         labels = label.levels, values = shape.levels) +
      scale_fill_manual(name = "pointer year class", limits = int.levels, 
                        labels = label.levels, values = fill.levels) +
      scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                         minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                         limits = c(start.yr3, end.yr3)) +
      theme_bw() + theme(legend.key = element_blank())
  }
  else {
    nature2 <- nature[as.character(start.yr2:end.yr2),]
    input <- matrix2long(t(nature2), new.ids = FALSE)
    input2 <- na.omit(input)
    rownames(input2) <- NULL
    colnames(input2) <-c ("site", "year", "PYvalues")
      
    if(sign2 == "both") {
      int.levels <- c(-1, 0, 1)
      label.levels <- c("negative", "none", "positive")
      shape.levels <- c(25, 95, 24)
    }
    if(sign2 == "pos") {
      input2[input2$PYvalues < 0, "PYvalues"] <- 0
      int.levels <- c(0, 1)
      label.levels <- c("other", "positive")
      shape.levels <- c(95, 24)
    }
    if(sign2 == "neg") {
      input2[input2$PYvalues > 0, "PYvalues"] <- 0
      int.levels <- c(-1, 0)
      label.levels <- c("negative", "other")
      shape.levels <- c(25, 95)
    }
    
    pl <- ggplot(input2, aes(x = year, y = site, shape = factor(PYvalues))) 
    pl + geom_point(size = 2, colour = "black", fill = "#bdbdbd") +
      scale_shape_manual(name = "event year", limits = int.levels,
                         labels = label.levels, values = shape.levels) + 
      scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                         minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                         limits = c(start.yr3, end.yr3)) +
      theme_bw() + theme(legend.key = element_blank())
  }
}






