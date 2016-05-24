#' @name segment
#' @title Segmented Hockey Stick Test
#' @usage
#' segment(dosecolumn = "", targetcolumn = "", data = NA)
#' @description
#' This function returns a two-segment linear dose-response model, by
#' imposing a zero slope for the initial (left) segment,
#' detecting one breakpoint where the dose-response relation changes 
#' to a positive slope (if such a breakpoint exists), then reporting 
#' the breakpoint dose, its standard error and p-value, and plotting the model.
#'
#' @details
#' This function:
#'
#' 1) Attempts to identify one breakpoint using the 'segmented' function,
#' starting the search at the median of the input dose variable,
#' and imposing a zero-slope left-hand segment before any identified breakpoint.
#'
#' 2) If a breakpoint has been identified using this iterative approach, the p-value
#' is returned and model plotted; otherwise the function returns no breakpoint.
#' 
#' This function is currently only intended for use on continuous outcome data.
#' @param dosecolumn   Name of dose column of interest in dataframe.
#' @param targetcolumn  Name of response column of interest in dataframe.
#' @param data   Input dataframe.
#' @return
#' Returns the estimated breakpoint, standard error, and 90 percent confidence limits on the breakpoint,
#' as well as a plot of the estimated two-segment dose-response relationship.
#' @examples
#' # Prints the breakpoint, its standard error, and 95% confidence limits
#' # along with a plot of the estimated two-segment linear relationship:
#' segment("dose", "MF_Log", data=DRdata)
#' @export

segment <- function (dosecolumn="", targetcolumn="", data=NA) {

     target<- data[,targetcolumn]
     fit.lm <- stats::lm(target ~ 1, data=data)
     dose <- data[,dosecolumn]
    
    out <- tryCatch(
        out <- fit.seg <- segmented::segmented(fit.lm, seg.Z = ~ dose, psi = stats::median(dose), control=segmented::seg.control(it.max = 100)),
    
        error = function(cond) {
        message("No breakpoint could be estimated from these data within the dose range.")
        return(NA)
        }
    )
    if ("segmented" %in% class(out))
    {
    Plot <- get("segment.plot", envir = environment(drsmooth))
    Plot(mod = fit.seg, dosecolumn=dosecolumn, targetcolumn=targetcolumn, data=data)
    Print <- get("segment.print", envir = environment(drsmooth))
    Print(fit.seg)
    }
}
    



