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
