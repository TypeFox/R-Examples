"DynamicGraph" <-
function (names = character(), types = character(), from = vector(), 
    to = vector(), edge.list = list(NULL), labels = names, edge.types = character(), 
    blocks = list(NULL), block.tree = list(NULL), oriented = NA, 
    factors = list(NULL), texts = character(), extra.from = vector(), 
    extra.to = vector(), extra.edge.list = list(NULL), object = NULL, 
    viewType = "Simple", frameModels = NULL, frameViews = NULL, 
    graphWindow = NULL, addModel = FALSE, addView = FALSE, overwrite = FALSE, 
    returnNewMaster = FALSE, redraw = FALSE, control = dg.control(...), 
    ...) 
{
    wDG(vertex.names = names, types = types, from = from, to = to, 
        edge.types = edge.types, edge.list = edge.list, labels = labels, 
        blocks = blocks, block.tree = block.tree, oriented = oriented, 
        factors = factors, texts = texts, extra.from = extra.from, 
        extra.to = extra.to, extra.edge.list = extra.edge.list, 
        object = object, viewType = viewType, frameModels = frameModels, 
        frameViews = frameViews, graphWindow = graphWindow, addModel = addModel, 
        addView = addView, overwrite = overwrite, returnNewMaster = returnNewMaster, 
        redraw = redraw, control = control, ...)
}
