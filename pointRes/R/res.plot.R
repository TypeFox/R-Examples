#' Plot resilience components
#'
#' @description The function creates box plots of the resilience components resistance, recovery, resilience and relative resilience as produced by \code{\link{res.comp}} for years identified as negative pointer years, as well as for selected years. 
#' 
#' @usage res.plot(list.name, select.yr = NULL, multi.panel = TRUE)
#'
#' @param list.name a \code{list} as produced by \code{\link{res.comp}}.
#' @param select.yr an \code{integer} specifying the (pointer) years to be plotted (e.g., c(1948, 1992)). Defaults to all years defined as negative pointer year with \code{\var{nb.series}} >= 5 in the \code{list} component \code{out.select}.
#' @param multi.panel a \code{logical} specifying whether box plots should be plotted in a 2x2 grid. Defaults to TRUE.
#' 
#' @details The function makes a box plot for each resilience component showing the full range of variation for individual trees in negative pointer years (or selected years). Box plots are only created for years with \code{\var{nb.series}} >= 5, as this value represents the number of statistics that a box plot represents in its' simplest form.
#'
#' @return 
#' Four box plots.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot resilience components for all defined pointer years
#' # note: pointer years with < 5 series (here 1882) are not displayed (warning)
#' data(s033)
#' res <- res.comp(s033, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)
#' res.plot(res, select.yr = NULL, multi.panel = TRUE)
#'
#' ## Plot resilience components for selected years
#' # note: inclusion of non-pointer years (here 2002) results in a warning
#' data(s033)
#' res <- res.comp(s033, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)
#' res.plot(res, select.yr = c(1948, 1992, 2002), multi.panel = TRUE)
#' 
#' @import ggplot2
#' @import stats
#' @importFrom gridExtra grid.arrange
#' @importFrom TripleR matrix2long
#' @import graphics
#' 
#' @export res.plot
#' 
res.plot <- function(list.name, select.yr = NULL, multi.panel = TRUE)
{
  stopifnot(is.list(list.name))
  if(class(list.name) != "res.comp") {
    stop("'list.name' is no list output of function res.comp")
  }
  if(is.data.frame(list.name$out.select) == FALSE) {
    stop("'list.name' is no list output of function res.comp")
  }
  if(is.null(select.yr) && nrow(list.name$out.select) == 0) {
    stop("no pointer years to display here")
  }
  if(is.null(select.yr) && all(list.name$out.select[,"nb.series"] < 5)) {
    stop("all pointer years have < 5 series and are not displayed")
  }
  if(is.null(select.yr) && TRUE %in% (list.name$out.select[,"nb.series"] < 5)) {
    warning("pointer years with < 5 series are not displayed")
  }
  if(!is.null(select.yr) && all(subset(list.name$out, list.name$out[,"year"] %in% select.yr)[,"nb.series"] < 5)) {
    stop("all selected years have < 5 series and are not displayed")
  }
  if(FALSE %in% (select.yr %in% list.name$out[,"year"])) {
    stop("the selection under 'select.yr' contains year(s) that are not in the dataset")
  }
  if(TRUE %in% ((subset(list.name$out,
                         list.name$out[,"year"] %in% select.yr)[,"nb.series"] < 5))) {
    warning("selected years with < 5 series are not displayed")
  }
  if(FALSE %in% (select.yr %in% list.name$out.select[,"year"])) {
    warning("the selection under 'select.yr' contains year(s) not identified as pointer year(s)")
  }

  year <- value <- NULL
  
  if(length(select.yr) == 0) {
    data2 <- subset(list.name$out.select,list.name$out.select[, "nb.series"] >= 5)
    if(nrow(data2) == 1) {
      yrs <- as.character(data2[, "year"])
      resist2 <- matrix2long(list.name$resist[yrs,], new.ids = FALSE)
      recov2 <- matrix2long(list.name$recov[yrs,], new.ids = FALSE)
      resil2 <- matrix2long(list.name$resil[yrs,], new.ids = FALSE)
      rel.resil2 <- matrix2long(list.name$rel.resil[yrs,], new.ids = FALSE)
      colnames(resist2) <- colnames(recov2) <- colnames(resil2) <- colnames(rel.resil2) <- c("tree","year","value")
    }
    else {
      yrs <- as.character(data2[, "year"])
      resist2 <- matrix2long(t(list.name$resist[yrs,]), new.ids = FALSE)
      recov2 <- matrix2long(t(list.name$recov[yrs,]), new.ids = FALSE)
      resil2 <- matrix2long(t(list.name$resil[yrs,]), new.ids = FALSE)
      rel.resil2 <- matrix2long(t(list.name$rel.resil[yrs,]), new.ids = FALSE)
      colnames(resist2) <- colnames(recov2) <- colnames(resil2) <- colnames(rel.resil2) <- c("tree","year","value")
    }
  }
  else {
    if(length(select.yr) == 1) {
      data2 <- subset(list.name$out,list.name$out[,"year"] %in% select.yr & list.name$out[, "nb.series"] >= 5)
      yrs <- as.character(data2[, "year"])
      resist2 <- matrix2long(list.name$resist[yrs,], new.ids = FALSE)
      recov2 <- matrix2long(list.name$recov[yrs,], new.ids = FALSE)
      resil2 <- matrix2long(list.name$resil[yrs,], new.ids = FALSE)
      rel.resil2 <- matrix2long(list.name$rel.resil[yrs,], new.ids = FALSE)
      colnames(resist2) <- colnames(recov2) <- colnames(resil2) <- colnames(rel.resil2) <- c("tree","year","value")  
    }
    else {
      data2 <- subset(list.name$out,list.name$out[,"year"] %in% select.yr & list.name$out[, "nb.series"] >= 5)
      yrs <- as.character(data2[, "year"])
      resist2 <- matrix2long(t(list.name$resist[yrs,]), new.ids = FALSE)
      recov2 <- matrix2long(t(list.name$recov[yrs,]), new.ids = FALSE)
      resil2 <- matrix2long(t(list.name$resil[yrs,]), new.ids = FALSE)
      rel.resil2 <- matrix2long(t(list.name$rel.resil[yrs,]), new.ids = FALSE)
      colnames(resist2) <- colnames(recov2) <- colnames(resil2) <- colnames(rel.resil2) <- c("tree","year","value")  
    }
  }

    res1 <- ggplot(na.omit(resist2), aes(x = factor(year), y = value)) + 
      geom_boxplot(fill = "#f0f0f0") + theme_bw() + guides(fill = FALSE) +
      xlab("year") + ylab("resistance index") + scale_x_discrete(labels = yrs)
    res2 <- ggplot(na.omit(recov2), aes(x = factor(year), y = value)) + 
      geom_boxplot(fill = "#f0f0f0") + theme_bw() + guides(fill = FALSE) +
      xlab("year") + ylab("recovery index") + scale_x_discrete(labels = yrs)
    res3 <- ggplot(na.omit(resil2), aes(x = factor(year), y = value)) + 
      geom_boxplot(fill = "#f0f0f0") + theme_bw() + guides(fill = FALSE) +
      xlab("year") + ylab("resilience index") + scale_x_discrete(labels = yrs)
    res4 <- ggplot(na.omit(rel.resil2), aes(x = factor(year), y = value)) + 
      geom_boxplot(fill = "#f0f0f0") + theme_bw() + guides(fill = FALSE) +
      xlab("year") + ylab("relative resilience index") + scale_x_discrete(labels = yrs)
    
  if(multi.panel) {
    grid.arrange(res1, res2, res3, res4)
    }
    else {
      plot(res1)
      plot(res2)
      plot(res3)
      plot(res4)
    }
}




