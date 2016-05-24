#' Plot mean relative growth changes and pointer years
#'
#' @description The function creates a bar plot of mean relative growth changes from a \code{list} of the type as produced by \code{\link{pointer.rgc}} and highlights years identified as pointer years.
#' 
#' @usage rgc.plot(list.name, start.yr = NULL, end.yr = NULL, 
#'          sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
#'
#' @param list.name a \code{list} as produced by \code{\link{pointer.rgc}}
#' @param start.yr an \code{integer} specifying the first year to be plotted. Defaults to the first year with data if \code{\var{start.yr}} is \code{NULL}.
#' @param end.yr an \code{integer} specifying the last year to be plotted. Defaults to the last year with data if \code{\var{end.yr}} is \code{NULL}.
#' @param sd.disp a \code{logical} specifying whether error bars (stdev) should be displayed. Defaults to FALSE.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a plot showing mean relative growth changes; pointer years are indicated with dark-gray bars. Error bars can be set.
#'
#' @return 
#' Bar plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot mean relative growth changes and pointer years
#' data(s033)
#' py <- pointer.rgc(s033, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, 
#'                   series.thresh = 75)
#' rgc.plot(py, start.yr = 1950, end.yr = NULL,  
#'          sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
#' 
#' @import ggplot2
#' @importFrom plyr round_any
#' 
#' @export rgc.plot
#' 
rgc.plot <- function(list.name, start.yr = NULL, end.yr = NULL,
                     sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
{
  stopifnot(is.list(list.name))
  if(class(list.name)[1] != "pointer.rgc") {
    stop("'list.name' is no list output of function pointer.rgc")
  }
  if(is.data.frame(list.name$out) == FALSE) {
    stop("'list.name' is no list output of function pointer.rgc")
  }
  if("dev_mean" %in% colnames(list.name$out) == FALSE) {
    stop("'list.name' is no list output of function pointer.rgc")
  }
  if(nrow(list.name$out) < 2) {
    stop("'list.name'$out contains < 2 years and no plot is created")
  }
  if(!is.null(start.yr) && start.yr < min(list.name$out[, "year"])) {
    stop("'start.yr' is out of bounds. By default (start.yr = NULL) the first year is displayed")
  }
  if(!is.null(end.yr) && end.yr > max(list.name$out[, "year"])) {
    stop("'end.yr' is out of bounds. By default (end.yr = NULL) the last year is displayed")
  }
  if(x.tick.minor > x.tick.major) {
    stop("'x.tick.minor' should be smaller then 'x.tick.major'")
  }
  
  start.yr2 <- ifelse(length(start.yr) != 0, start.yr, min(list.name$out[, "year"])) 
  end.yr2 <- ifelse(length(end.yr) != 0, end.yr, max(list.name$out[, "year"]))
  start.yr3 <- round_any(start.yr2, 10, f = floor)
  end.yr3 <- round_any(end.yr2, 5, f = ceiling)
  
  data2 <- list.name$out[which(list.name$out[, "year"] == start.yr2):which(list.name$out[, "year"] == end.yr2),]
  data3 <- as.data.frame(data2)
  
  year <- nature <- dev_mean <- dev_sd <- NULL

  nat.levels <- c(-1, 0, 1)
  fill.levels <- c("#636363", "#f0f0f0", "#636363")
  
  if(sd.disp) {
  limits <- aes(ymax = dev_mean + dev_sd, ymin = dev_mean - dev_sd)
  
  pl <- ggplot(data3, aes(x = year, y = dev_mean, fill = factor(nature))) 
  pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
    scale_fill_manual(limits = nat.levels, values = fill.levels) +
    guides(fill = FALSE) +
    scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                       minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                       limits = c(start.yr3-1, end.yr3+1)) +
    ylab("mean growth deviation (%)") + theme_bw() +
    geom_errorbar(limits, width=0.25, colour = "gray60")
  }
  else {
    pl <- ggplot(data3, aes(x = year, y = dev_mean, fill = factor(nature))) 
    pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
      scale_fill_manual(limits = nat.levels, values = fill.levels) +
      guides(fill = FALSE) +
      scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                         minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                         limits = c(start.yr3-1, end.yr3+1)) +
      ylab("mean growth deviation (%)") + theme_bw()
  }
}


