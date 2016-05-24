"dg.control" <-
function (label = "dynamicGraph", width = 400, height = 400, 
    w = 6, margin = 100, closeenough = 2, background = "white", 
    transformation = NULL, permitZoom = TRUE, UserMenus = NULL, 
    constrained = FALSE, vertexColor = "red", extraVertexColor = "white", 
    edgeColor = "black", factorVertexColor = "default", factorEdgeColor = "brown", 
    blockEdgeColor = "default", blockColors = NULL, extraEdgeColor = "peru", 
    drawblocks = TRUE, right.to.left = FALSE, nested.blocks = FALSE, 
    overlaying = TRUE, fixedFactorPositions = FALSE, diagonal = TRUE, 
    N = 3, vertexClasses = validVertexClasses(), factorClasses = validFactorClasses(), 
    edgeClasses = validEdgeClasses(), viewClasses = validViewClasses(), 
    drawBlockFrame = TRUE, drawBlockBackground = FALSE, useNamesForLabels = TRUE, 
    namesOnEdges = TRUE, updateEdgeLabels = TRUE, enterLeaveUpdate = TRUE, 
    updateAllViews = TRUE, saveTkReferences = TRUE, saveFunctions = TRUE, 
    returnNull = FALSE, hasMethods = TRUE, variableFrame = TRUE, 
    debug.strata = FALSE, debug.edges = FALSE, debug.position = FALSE, 
    debug.update = FALSE, ...) 
{
    args <- list(...)
    if (!is.null(args$returnLink)) {
        if (!saveTkReferences) 
            saveTkReferences <- args$returnLink
        if (!saveFunctions) 
            saveFunctions <- as.logical(args$returnLink) && (as.numeric(args$returnLink) < 
                2)
    }
    if (!is.null(args$title) && (label == "dynamicGraph")) {
        label <- args$title
    }
    list(label = label, width = width, height = height, w = w, 
        margin = margin, closeenough = closeenough, background = background, 
        transformation = transformation, permitZoom = permitZoom, 
        UserMenus = UserMenus, constrained = constrained, vertexColor = vertexColor, 
        extraVertexColor = extraVertexColor, edgeColor = edgeColor, 
        factorVertexColor = factorVertexColor, factorEdgeColor = factorEdgeColor, 
        blockEdgeColor = blockEdgeColor, blockColors = blockColors, 
        extraEdgeColor = extraEdgeColor, drawblocks = drawblocks, 
        right.to.left = right.to.left, nested.blocks = nested.blocks, 
        overlaying = overlaying, fixedFactorPositions = fixedFactorPositions, 
        diagonal = diagonal, N = N, vertexClasses = vertexClasses, 
        factorClasses = factorClasses, edgeClasses = edgeClasses, 
        viewClasses = viewClasses, drawBlockFrame = drawBlockFrame, 
        drawBlockBackground = drawBlockBackground, useNamesForLabels = useNamesForLabels, 
        namesOnEdges = namesOnEdges, updateEdgeLabels = updateEdgeLabels, 
        enterLeaveUpdate = enterLeaveUpdate, updateAllViews = updateAllViews, 
        saveTkReferences = saveTkReferences, saveFunctions = saveFunctions, 
        returnNull = returnNull, variableFrame = variableFrame, 
        hasMethods = hasMethods, debug.strata = debug.strata, 
        debug.edges = debug.edges, debug.position = debug.position, 
        debug.update = debug.update)
}
