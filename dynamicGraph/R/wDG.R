"wDG" <-
function (sdg = NULL, object = NULL, frameModels = NULL, frameViews = NULL, 
    graphWindow = NULL, dg = NULL, addModel = FALSE, addView = FALSE, 
    overwrite = FALSE, returnNewMaster = FALSE, redraw = FALSE, 
    control = dg.control(...), ...) 
{
    Arguments <- list(...)
    ".wDG" <- function(localVertices = NULL, localBlockList = NULL, 
        dg = NULL, addArguments = FALSE, ...) {
        if (addModel) {
            if (is.null(frameModels)) 
                warning("'frameModels' should be set;")
            if (is.object(frameModels)) {
                frameModels.env <- .get.env.frameModels(frameModels = frameModels)
                drawModel <- frameModels.env$env$DrawModel
                if (length(formals(drawModel)) == 0) 
                  drawModel <- Arguments$Arguments$drawModel
            }
            else drawModel <- frameModels$drawModel
            if (!overwrite) {
                frameViews <- NULL
                graphWindow <- NULL
            }
            else {
                if (addArguments && !is.null(frameViews)) 
                  frameViews@model <- list(object)
            }
            if (length(formals(drawModel)) == 0) {
                warning("Invalid 'drawModel' of argument")
            }
            else drawModel(frameModels = frameModels, frameViews = frameViews, 
                graphWindow = graphWindow, dg = dg, object = object, 
                control = control, ...)
        }
        else if (addView) {
            if (is.null(frameModels)) 
                warning("'frameModels' should be set;")
            if (addArguments && !is.null(object)) 
                warning("'object' ignored!  Did you mean 'addModel'? ; ")
            if (is.object(frameViews)) {
                frameViews.env <- .get.env.frameViews(frameViews = frameViews, 
                  frameModels = frameModels)
                redrawView <- frameViews.env$env$RedrawView
                if (length(formals(redrawView)) == 0) 
                  redrawView <- Arguments$Arguments$redrawView
            }
            else redrawView <- frameViews$redrawView
            if (addArguments && !overwrite) {
                graphWindow <- NULL
            }
            else {
            }
            if (length(formals(redrawView)) == 0) {
                warning("Invalid 'redrawView' of argument")
            }
            else redrawView(frameModels = frameModels, frameViews = frameViews, 
                dg = dg, object = object, control = control, 
                ...)
        }
        else {
            if (!redraw && !returnNewMaster) 
                frameModels <- NULL
            dynamicGraphMain(vertexList = localVertices, blockList = localBlockList, 
                dg = dg, object = object, frameModels = frameModels, 
                returnNewMaster = returnNewMaster, redraw = redraw, 
                control = control, ...)
        }
    }
    X <- c("edgeList", "blockEdgeList", "factorVertexList", "factorEdgeList", 
        "extraList", "extraEdgeList")
    if (any(X %in% names(list(...)))) {
        if (!is.null(sdg)) {
            if (!.IsEmpty(sdg@from)) 
                warning("Argument 'from      ' ignored")
            if (!.IsEmpty(sdg@to)) 
                warning("Argument 'to        ' ignored")
            if (!.IsEmpty(sdg@edge.types)) 
                warning("Argument 'edge.types' ignored")
            if (!.IsEmpty(sdg@edge.list)) 
                warning("Argument 'edge.list ' ignored")
            if (!.IsEmpty(sdg@blocks)) 
                warning("Argument 'blocks    ' ignored")
            if (!.IsEmpty(sdg@block.tree)) 
                warning("Argument 'block.tree' ignored")
            if (!.IsEmpty(sdg@factors)) 
                warning("Argument 'factors   ' ignored")
            if (!.IsEmpty(sdg@texts)) 
                warning("Argument 'texts     ' ignored")
        }
        .wDG(dg = NULL, addArguments = TRUE, ...)
    }
    else {
        dg <- simpleGraphToGraph(sdg = sdg, frameModels = frameModels, 
            dg = if (!is.null(dg)) 
                dg
            else if (!addModel && !addView) 
                if (!is.null(graphWindow)) 
                  graphWindow@dg, ...)
        .wDG(localVertices = dg@vertexList, localBlockList = dg@blockList, 
            dg = dg, ...)
    }
}
