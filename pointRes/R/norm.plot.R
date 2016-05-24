#' Plot mean Cropper values and pointer years
#'
#' @description The function creates a bar plot of mean Cropper values from a \code{list} of the type as produced by \code{\link{pointer.norm}} and highlights years identified as pointer years.
#' 
#' @usage norm.plot(list.name, start.yr = NULL, end.yr = NULL, 
#'           sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
#'
#' @param list.name a \code{list} as produced by \code{\link{pointer.norm}}
#' @param start.yr an \code{integer} specifying the first year to be plotted. Defaults to the first year included in the \code{out} component of the \code{list} if \code{\var{start.yr}} is \code{NULL}.
#' @param end.yr an \code{integer} specifying the last year to be plotted. Defaults to the last year included in the \code{out} component of the \code{list} if \code{\var{end.yr}} is \code{NULL}.
#' @param sd.disp a \code{logical} specifying whether error bars (stdev) should be displayed. Defaults to FALSE.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a plot showing mean Cropper values; pointer years are indicated with dark-gray bars. If event years were defined using \code{method.thresh "Neuwirth"} (\code{\link{pointer.norm}}), different tones of gray indicate weak, strong and extreme pointer years, based on the most common event year class. Error bars can be set.
#'
#' @return 
#' Bar plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#' 
#' @examples ## Plot mean Cropper values and pointer years (method "Cropper")
#' data(s033)
#' py_c <- pointer.norm(s033, window = 5, method.thresh = "Cropper", 
#'                      series.thresh = 75)
#' norm.plot(py_c, start.yr = 1950, end.yr = NULL, 
#'           sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
#'
#' ## Plot mean Cropper values and pointer years (method "Neuwirth")
#' data(s033)
#' py_n <- pointer.norm(s033, window = 5, method.thresh = "Neuwirth",
#'                      series.thresh = 75)
#' norm.plot(py_n, start.yr = 1950, end.yr = NULL, 
#'           sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
#'           
#' @import ggplot2
#' @importFrom plyr round_any
#' 
#' @export norm.plot
#' 
norm.plot <- function(list.name, start.yr = NULL, end.yr = NULL,
                      sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
{
  stopifnot(is.list(list.name))
  if(class(list.name)[1] != "pointer.norm") {
    stop("'list.name' is no list output of function pointer.norm")
  }
  if(is.data.frame(list.name$out) == FALSE) {
    stop("'list.name' is no list output of function pointer.norm")
  }
  if("Cvalues_mean" %in% colnames(list.name$out) == FALSE) {
    stop("'list.name' is no list output of function pointer.norm")
  }
  if(nrow(list.name$out) < 2){
    stop("'list.name'$out contains < 2 years and is not displayed")
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
  
  year <- nature <- Cvalues_mean <- Cvalues_sd <- int.class <- NULL
  limits <- aes(ymax = Cvalues_mean + Cvalues_sd, ymin = Cvalues_mean - Cvalues_sd)
  
  if(colnames(data3)[3] == "perc.pos") {
    nat.levels <- c(-1, 0, 1)
    fill.levels <- c("#636363", "#f0f0f0", "#636363")
    
    if(sd.disp) {
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(nature))) 
      pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = nat.levels, values = fill.levels) +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                           minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                           limits = c(start.yr3-1, end.yr3+1)) +
        ylab("mean Cropper value") + theme_bw() + 
        geom_errorbar(limits, width=0.25, colour = "gray60")
    }
    else {
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(nature))) 
      pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = nat.levels, values = fill.levels) +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                           minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                           limits = c(start.yr3-1, end.yr3+1)) +
        ylab("mean Cropper value") + theme_bw()
    }
  }
  else {
    data3[,12] <- ifelse(data2[, "nature"] == (-1), 
                         max.col(data2[,c(1, 2, 5, 4, 3, 6:11)][, 6:8], ties.method = "first"), 
                         ifelse(data2[, "nature"] == 1,
                                max.col(data2[,c(1, 2, 5, 4, 3, 6:11)][, 3:5], ties.method = "first"), 0))
    data3[,12] <- ifelse(data2[, "nature"] == (-1), paste("-", data3[, 12], sep = ''), data3[, 12])
    colnames(data3)[12] <- "int.class"
    
    int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
    fill.levels <- c("black", "#636363","#bdbdbd", "#f0f0f0", "#bdbdbd", "#636363", "black")
    
    if(sd.disp) {
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(int.class))) 
      pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = int.levels, values = fill.levels) +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                           minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                           limits = c(start.yr3-1, end.yr3+1)) +
        ylab("mean Cropper value") + theme_bw() + 
        geom_errorbar(limits, width=0.25, colour = "gray60")
    }
    else {
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(int.class))) 
      pl + geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = int.levels, values = fill.levels) +
        guides(fill = FALSE) +
        scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                           minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                           limits = c(start.yr3-1, end.yr3+1)) +
        ylab("mean Cropper value") +  theme_bw()
    }
  }
}