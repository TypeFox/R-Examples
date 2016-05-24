"replaceControls" <-
function (control, frameModels = NULL, frameViews = frameModels@models[[modelIndex]], 
    modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
    viewIndex = 1, ...) 
{
    dots <- list(...)
    localArguments <- dots$Arguments
    redrawView <- .get.redrawView(frameViews, localArguments)
    dg <- graphWindow@dg
    if (length(formals(redrawView)) == 0) {
        warning("Invalid 'redrawView' of argument")
    }
    else redrawView(frameModels = frameModels, frameViews = frameViews, 
        graphWindow = graphWindow, dg = dg, control = control, 
        Arguments = localArguments)
}
