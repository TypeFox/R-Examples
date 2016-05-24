"replaceBlockList" <-
function (blockList, frameModels = NULL, frameViews = frameModels@models[[modelIndex]], 
    modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
    viewIndex = 1, ...) 
{
    dots <- list(...)
    localArguments <- dots$Arguments
    redrawView <- .get.redrawView(frameViews, localArguments)
    dg <- graphWindow@dg
    control <- frameModels@control
    if (length(formals(redrawView)) == 0) {
        warning("Invalid 'redrawView' of argument")
    }
    else redrawView(frameModels = frameModels, frameViews = frameViews, 
        graphWindow = graphWindow, blockList = blockList, dg = dg, 
        control = control, Arguments = localArguments)
}
