#' @name segment.plot
#' @title Plot Segmented Test
#' @description
#' Plots the segmented function, if successful.
#' @param mod  Segmented model object
#' @param dosecolumn   Name of dose column.
#' @param targetcolumn   Name of response column.
#' @param data   Input dataframe.
#' @keywords internal

segment.plot <- function (mod=NA, dosecolumn="", targetcolumn="", data=NA) {

    graphics::plot(data[,dosecolumn], data[,targetcolumn], type="p", xlab = dosecolumn, ylab = targetcolumn, pch=19)
    segmented::plot.segmented(mod, add=TRUE)

}