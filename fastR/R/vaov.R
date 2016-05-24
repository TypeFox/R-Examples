#' ANOVA vectors
#' 
#' Compute vectors associated with 1-way ANOVA
#' 
#' This is primarily designed for demonstration purposes to show how 1-way
#' ANOVA models partition variance.  It may not work properly for more
#' complicated models.
#' 
#' @aliases vaov vaov.formula
#' @param x a formula.
#' @param data a data frame.
#' @param \dots additional arguments.
#' @return A data frame with variables including \code{grandMean},
#' \code{groupMean}, \code{ObsVsGrand}, \code{STotal}, \code{ObsVsGroup},
#' \code{SError}, \code{GroupVsGrand}, and \code{STreatment}. The usual SS
#' terms can be computed from these by summing.
#' @author Randall Pruim
#' @keywords stats
#' @export
#' @examples
#' 
#' aov(pollution ~ location, data=airpollution)
#' vaov(pollution ~ location, data=airpollution)
#' 
vaov <-
function (x, ...) 
{
    UseMethod("vaov", x)
}


#' @rdname vaov
#' @method vaov formula
#' @export
vaov.formula <-
function (x, data = parent.frame(), ...) 
{
    groupMeans <- funvec(x, data, mean)
    form <- latticeParseFormula(x, data, ...)
    overallMeans <- rep(mean(form$left), length(form$left))
    df <- data.frame(form$right, form$left, overallMeans, groupMeans, 
        (form$left - overallMeans), (form$left - overallMeans)^2, 
        (form$left - groupMeans), (form$left - groupMeans)^2, 
        (groupMeans - overallMeans), (groupMeans - overallMeans)^2)
    names(df) = c(form$right.name, form$left.name, "GrandMean", 
        "GroupMean", "ObsVsGrand", "STotal", "ObsVsGroup", "SError", 
        "GroupVsGrand", "STreatment")
    return(df)
}
