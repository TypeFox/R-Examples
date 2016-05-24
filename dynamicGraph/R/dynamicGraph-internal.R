".addModel" <-
function (dg, modelObject = NULL, modelObjectName = NULL, frameModels = NULL, 
    frameViews = frameModels@models[[modelIndex]], modelIndex = 1, 
    graphWindow = frameViews@graphs[[viewIndex]], viewIndex = 1, 
    overwrite = FALSE, control = NULL, ...) 
{
    dots <- list(...)
    localArguments <- dots$Arguments
    if (is.null(frameModels)) 
        warning("'frameModels' should be set;")
    if (is.object(frameModels)) {
        frameModels.env <- .get.env.frameModels(frameModels = frameModels)
        drawModel <- frameModels.env$env$DrawModel
    }
    if (length(formals(drawModel)) == 0) 
        drawModel <- localArguments$drawModel
    if (overwrite) {
        if (!is.null(frameViews)) 
            frameViews@model <- list(modelObject)
    }
    else {
        frameViews <- NULL
        graphWindow <- NULL
    }
    if (is.null(control)) 
        control <- frameModels@control
    label <- control$label
    if (any(slotNames(modelObject) == ".title")) {
      label <- modelObject@.title
    }
    if (length(formals(drawModel)) == 0) {
        warning("Invalid function 'drawModel' of arguments")
    }
    else drawModel(frameModels = frameModels, frameViews = frameViews, 
        graphWindow = graphWindow, dg = dg, object = modelObject, 
        objectName = modelObjectName, control = control, title = label, ...)
}

".addView" <-
function (dg, frameModels = NULL, 
    frameViews = frameModels@models[[modelIndex]], 
    modelIndex = 1, graphWindow = frameViews@graphs[[viewIndex]], 
    viewIndex = 1, overwrite = FALSE, control = NULL, ...) 
{
    dots <- list(...)
    localArguments <- dots$Arguments
    redrawView <- .get.redrawView(frameViews, localArguments)
    if (overwrite) {
    }
    else {
        graphWindow <- NULL
    }
    if (is.null(control)) 
        control <- frameModels@control
    if (length(formals(redrawView)) == 0) {
        warning("Invalid function 'redrawView' of arguments")
    }
    else redrawView(frameModels = frameModels, frameViews = frameViews, 
        graphWindow = graphWindow, dg = dg, control = control, 
        Arguments = localArguments)
}

".asDataFrame" <-
function (v, setRowLabels = FALSE, shortNames = TRUE, excludeSlots = FALSE) 
{
    ".removeNull" <- function(x) {
        result <- NULL
        for (i in seq(along = x)) 
            if (!is.null(x[[i]]) && (length(x[[i]] > 0))) 
            result <- append(result, list(x[[i]]))
        return(result)
    }
    "f" <- function(name, v) {
        if (name == "class") 
            x <- lapply(v, function(i) slot(i, name))
        else x <- lapply(v, function(i) if (name %in% slotNames(i)) 
            slot(i, name))
        w <- unique(range(unlist(lapply(x, length))))
        if (length(w) > 1) 
            x <- lapply(x, function(i) c(i, rep(NA, max(w) - 
                length(i))))
        w <- max(w)
        X <- unlist(x)
        if (shortNames) {
            if (name == "vertex.indices") 
                name <- "v..i."
            if (name == "label.position") 
                name <- "l..pos."
            if (name == "position") 
                name <- "p."
        }
        if (name == "dist") {
            X <- matrix(paste(X), ncol = 1)
            dimnames(X) <- list(NULL, name)
        }
        else if (w > 1) {
            X <- matrix(X, ncol = w, byrow = TRUE)
            dimnames(X) <- list(NULL, paste(name, 1:w, sep = "."))
        }
        else {
            X <- matrix(X, ncol = 1)
            dimnames(X) <- list(NULL, name)
        }
        names(X) <- NULL
        return(X)
    }
    names.all <- unique(sort(unlist(lapply(v, function(i) slotNames(i)))))
    names.all <- c(names.all, "class")
    vertex.names <- slotNames(v[[1]])
    vertex.names <- c(vertex.names, names.all[is.na(match(names.all, 
        vertex.names))])
    if (is.character(excludeSlots) || excludeSlots) {
        if (is.logical(excludeSlots)) 
            excludeSlots <- c("class", "name", "label.position", 
                "descendants", "ancestors")
        if (("label" %in% excludeSlots) && setRowLabels) 
            excludeSlots <- excludeSlots[excludeSlots != "label"]
        vertex.names <- vertex.names[is.na(match(vertex.names, 
            excludeSlots))]
    }
    Y <- data.frame(.removeNull(lapply(vertex.names, function(name) f(name, 
        v))))
    if (setRowLabels && ("label" %in% vertex.names)) {
        labels <- t(Y["label"])
        if (length(unique(match(labels, labels))) < length(labels)) 
            labels <- paste(labels, 1:length(labels), sep = ".")
        dimnames(Y)[[1]] <- labels
        if (is.character(excludeSlots) || excludeSlots) 
            Y <- Y[!(dimnames(Y)[[2]] == "label")]
    }
    return(Y)
}

".asRow" <-
function (positions) 
if (is.null(dim(positions))) return(positions) else return(t(positions))
".but.last" <-
function (x, n = nchar(x), sep = ".") 
{
    p <- (1:n)[substring(x, 1:n, 1:n) == sep]
    b <- p[length(p)] - 1
    return(paste(substring(x, 1:b, 1:b), collapse = ""))
}

".cliquesFromEdges" <-
function (Edges, Vertices, VisibleVertices) 
{
    require(ggm)
    e <- NodeIndices(Edges)
    if (length(e) > 0) {
        e <- lapply(e, function(egde) 
            if (sum(abs(egde)) > 0) egde)
        e <- .removeNull(e)
    }
    else e <- NULL
    factors <- NULL
    if (length(e) < 2) {
        if (length(e) == 1) 
            factors <- append(e, as.list(VisibleVertices))
        else if (length(VisibleVertices) > 0) 
            factors <- as.list(VisibleVertices)
    }
    else {
        n <- Names(Vertices)
        X <- matrix(rep(0, length(n)^2), ncol = length(n))
        lapply(e, function(i) X[i[1], i[2]] <<- 1)
        dimnames(X) <- list(n, n)
        X <- X[VisibleVertices, VisibleVertices]
        factors <- cliques(X + t(X)) # 'cliques' from 'ggm' ???
    }
    return(factors)
}

".dashReplaceMethod" <-
function (x, value) 
{
    "abort" <- function() {
        message("Invalid DASH PATTERN, see the Tcl/tk Reference Manual;")
        value <<- x@dash
    }
    if (!is.na(as.logical(value))) {
        if (as.logical(value)) 
            value <- "."
        else value <- " "
    }
    else if (is.character(value)) {
        a <- substring(value, 1:nchar(value), 1:nchar(value))
        A <- c(".", ",", "-", "_")
        if (any(a %in% A)) {
            A <- c(".", ",", "-", "_", " ")
            if (!all(a %in% A)) 
                abort()
        }
        else {
            A <- c(paste(0:9))
            if (any(a %in% A)) {
                A <- c(paste(0:9), " ")
                if (!all(a %in% A)) 
                  abort()
            }
            else abort()
        }
    }
    else abort()
    x@dash <- value
    x
}

".Dg.newenv" <-
function (ID) 
{
    win <- list(ID = ID, env = evalq(new.env(), .GlobalEnv))
    evalq(num.subenv <- 0, win$env)
    class(win) <- "dgenv"
    win
}

# ".DgRoot" <- structure(list(ID = "", env = <environment>), .Names = c("ID", "env"), class = "dgenv")

.DgRoot <- .Dg.newenv("")

".Dg.toplevel" <-
function (parent = .DgRoot, ...) 
{
    ".Dg.subenv" <- function(parent) {
        ID <- paste(parent$ID, evalq(num.subenv <- num.subenv + 
            1, parent$env), sep = ".")
        win <- .Dg.newenv(ID)
        assign(ID, win, envir = parent$env)
        assign("parent", parent, envir = win$env)
        win
    }
    ".Dg.ID" <- function(win) win$ID
    w <- .Dg.subenv(parent)
    ID <- .Dg.ID(w)
    w
}

".emptyDgList" <-
function (type = "dg.list", n = 0) 
new(type, vector("list", n))

".emptyToNull" <-
function (x) 
{
    if ((length(x) == 1) && is.null(x[[1]])) 
        return(NULL)
    else return(x)
}

".First.lib" <-
function (lib, pkg) 
{
    # .onLoad.dynamicGraph()
}

".get.env.frameModels" <-
function (id = frameModels@id.env, frameModels = NULL, env = .DgRoot$env) 
{
    if (is.null(env)) {
        message("NULL frameModels")
        NULL
    }
    else if (is.element(id, ls(env, all.names = TRUE))) {
        get(id, env)
    }
    else {
        message("Invalid frameModels")
        NULL
    }
}

".get.env.frameViews" <-
function (id = frameViews@id.env, frameViews = NULL, 
    env = .get.env.frameModels(id = id.fm, 
    env = .DgRoot$env)$env, 
    id.fm = if (is.null(frameModels)) .but.last(id) else frameModels@id.env, 
    frameModels = NULL) 
{
    if (is.null(env)) {
        message("NULL frameViews")
        NULL
    }
    else if (is.element(id, ls(env, all.names = TRUE))) {
        get(id, env)
    }
    else {
        message("Invalid frameViews")
        NULL
    }
}

".get.env.graphWindow" <-
function (id = graphWindow@id.env, graphWindow = NULL, 
    env = .get.env.frameViews(id = id.fv, 
    env = env.fm)$env, 
    id.fv = if (is.null(frameViews)) .but.last(id) else frameViews@id.env, 
    frameViews = NULL, 
    env.fm = .get.env.frameModels(id = id.fm, env = .DgRoot$env)$env, 
    id.fm = if (is.null(frameModels)) .but.last(.but.last(id)) 
            else frameModels@id.env, 
    frameModels = NULL) 
{
    if (is.null(env)) {
        message("NULL graphWindow")
        NULL
    }
    else if (is.element(id, ls(env, all.names = TRUE))) {
        get(id, env)
    }
    else {
        message("Invalid graphWindow")
        NULL
    }
}

".get.redrawView" <-
function (frameViews, localArguments) 
{
    redrawView <- NULL
    if (is.null(frameViews)) 
        warning("'frameViews' should be set;")
    if (is.object(frameViews)) {
        frameViews.env <- .get.env.frameViews(frameViews = frameViews)
        redrawView <- frameViews.env$env$RedrawView
    }
    if (length(formals(redrawView)) == 0) 
        redrawView <- localArguments$redrawView
    return(redrawView)
}

".is.dgenv" <-
function (x) 
inherits(x, "dgenv")

".IsEmpty" <-
function (x) 
{
    if (is.null(x) || (length(x) == 0) || (length(x) == 1) && 
        is.null(x[[1]])) 
        return(TRUE)
    else return(FALSE)
}

".newDgGraph" <-
function (viewType = "Simple", vertexList = NULL, 
    visibleVertices = 1:length(vertexList), 
    visibleBlocks = 1:length(blockList), oriented = NA, edgeList = NULL, 
    blockList = NULL, blockEdgeList = NULL, factorVertexList = NULL, 
    factorEdgeList = NULL, extraList = NULL, extraEdgeList = NULL, 
    ...) 
{
    new("dg.graph", viewType = viewType, 
        vertexList = .nullToEmpty(vertexList), 
        visibleVertices = .nullToEmpty(visibleVertices), 
        visibleBlocks = .nullToEmpty(visibleBlocks), 
        edgeList = .nullToList(edgeList,  
            type = "dg.VertexEdgeList"), 
        oriented = oriented, 
        blockList = .nullToList(blockList,  
            type = "dg.BlockList"), 
        blockEdgeList = .nullToList(blockEdgeList,  
            type = "dg.BlockEdgeList"), 
        factorVertexList = .nullToList(factorVertexList,  
            type = "dg.FactorVertexList"), 
        factorEdgeList = .nullToList(factorEdgeList, 
            type = "dg.FactorEdgeList"), 
        extraList = .nullToList(extraList, 
            type = "dg.VertexList"), 
        extraEdgeList = .nullToList(extraEdgeList, 
            type = "dg.ExtraEdgeList"))
}

".newDgGraphEdges" <-
function (viewType = "Simple", vertexList = NULL,  
   visibleVertices = 1:length(vertexList), 
    visibleBlocks = 1:length(blockList), oriented = NA, edgeList = NULL, 
    blockList = NULL, blockEdgeList = NULL, factorVertexList = NULL, 
    factorEdgeList = NULL, extraList = NULL, extraEdgeList = NULL, 
    ...) 
{
    new("dg.graphedges", viewType = viewType, 
        visibleVertices = .nullToEmpty(visibleVertices), 
        visibleBlocks = .nullToEmpty(visibleBlocks), 
        edgeList = .nullToList(edgeList, type = "dg.VertexEdgeList"), 
        oriented = oriented, 
        blockEdgeList = .nullToList(blockEdgeList, type = "dg.BlockEdgeList"), 
        factorVertexList = .nullToList(factorVertexList, type = "dg.FactorVertexList"), 
        factorEdgeList = .nullToList(factorEdgeList, type = "dg.FactorEdgeList"), 
        extraList = .nullToList(extraList, type = "dg.VertexList"), 
        extraEdgeList = .nullToList(extraEdgeList, type = "dg.ExtraEdgeList"))
}

".newDgSimpleGraph" <-
function (vertex.names = character(), labels = vertex.names, 
    types = character(), from = vector(), to = vector(), 
    edge.list = list(NULL), 
    edge.types = character(), blocks = list(NULL), block.tree = list(NULL), 
    oriented = NA, factors = list(NULL), texts = character(), 
    extra.from = vector(), extra.to = vector(), extra.edge.list = list(NULL), 
    viewType = "Simple", ...) 
{
    new("dg.simple.graph", vertex.names = vertex.names, types = types, 
        labels = labels, from = if (is.null(from)) 
            vector()
        else from, to = if (is.null(to)) 
            vector()
        else to, edge.list = if (is.null(edge.list)) 
            list(NULL)
        else edge.list, edge.types = edge.types, blocks = if (is.null(blocks)) 
            list(NULL)
        else blocks, block.tree = if (is.null(block.tree)) 
            list(NULL)
        else block.tree, oriented = oriented, factors = if (is.null(factors)) 
            list(NULL)
        else factors, texts = texts, extra.from = if (is.null(extra.from)) 
            vector()
        else extra.from, extra.to = if (is.null(extra.to)) 
            vector()
        else extra.to, extra.edge.list = if (is.null(extra.edge.list)) 
            list(NULL)
        else extra.edge.list, viewType = viewType)
}

".newDynamicGraphModelObject" <-
function (object, model = list(object), graphs = list(NULL), 
    label = "", index = 0, parent = "", env = .Dg.toplevel(parent)) 
{
    return(new("DynamicGraphModel", id.env = env$ID, label = label, 
        index = index, model = model, graphs = graphs))
}

".newDynamicGraphObject" <-
function (vertices, blocks = list(NULL), blockTree = list(NULL), 
    models = list(NULL), label = "", control = dg.control(), 
    parent = "", env = .Dg.toplevel(parent)) 
{
    ".nullToTree" <- function(x) if (is.null(x)) 
        list()
    else return(x)
    return(new("DynamicGraph", id.env = env$ID, label = label, 
        vertices = vertices, 
        blocks = .nullToList(blocks, type = "dg.BlockList"), 
        control = control, models = models))
}

".nullToEmpty" <-
function (x) 
if (is.null(x)) return(numeric()) else return(x)

".nullToList" <-
function (x, type = "dg.list") 
if (is.null(x)) .emptyDgList(type) else return(x)

".onAttach" <-
function (lib, pkg) 
{
    require(tcltk)
}

".onLoad" <-
function (lib, pkg) 
{
    # .onLoad.dynamicGraph()
}

".onLoad.dynamicGraph" <-
function () 
{
}

".onLoad.dynamicGraphInterface" <-
function () 
{
}

".propertyDialog" <-
function (object, classes = NULL, title = class(object), 
    sub.title = label(object), 
    name.object = name(object), okReturn = TRUE, fixedSlots = NULL, 
    difficultSlots = NULL, top = NULL, entryWidth = 20, do.grab = FALSE) 
{
    "subSelectDialog" <- function(dlg, itemNames, init = 0, title = title, 
        background = "white", subtitle = sub.title) {
        scr <- tkscrollbar(dlg, repeatinterval = 5, 
                           command = function(...) tkyview(tl, ...))
        tl <- tklistbox(dlg, height = length(itemNames), 
                        selectmode = "single", 
                        yscrollcommand = function(...) 
                                 tkset(scr, ...), background = background)
        tkgrid(tklabel(dlg, text = subtitle), tl, scr)
        tkgrid.configure(scr, rowspan = length(itemNames), sticky = "nsw")
        for (i in (1:length(itemNames))) tkinsert(tl, "end", 
            itemNames[i])
        if (!is.null(init)) 
            tkselection.set(tl, init)
        return(tl)
    }
    "subTextDialog" <- function(dlg, subtitle, entryInit, entryWidth = 20, 
        background = "white") {
        textEntryVarTcl <- tclVar(paste(entryInit))
        textEntryWidget <- tkentry(dlg, width = paste(entryWidth), 
            textvariable = textEntryVarTcl, background = background)
        tkgrid(tklabel(dlg, text = subtitle), textEntryWidget)
        return(textEntryVarTcl)
    }
    "subListDialog" <- function(dlg, subtitle, entryInit, entryWidth = 20, 
        background = "white") {
        r <- list()
        for (i in 1:length(entryInit)) {
            t <- subTextDialog(dlg, subtitle = paste(subtitle, 
                names(entryInit[i]), i, ": "), entryInit = entryInit[i], 
                entryWidth = entryWidth, background = background)
            r <- append(r, list(t))
        }
        return(r)
    }
    "onCancel" <- function() {
        tkgrab.release(dlg)
        tkdestroy(dlg)
        if (!is.null(top)) 
            tkfocus(top)
    }
    "onOK" <- function() {
        "getList" <- function(t.list) {
            r <- list()
            for (i in 1:length(t.list)) {
                t <- (tclvalue(t.list[[i]]))
                r <- append(r, list(t))
            }
            return(unlist(r))
        }
        result <<- lapply(t.list, function(t) {
            if (length(t) > 1) 
                r <- getList(t)
            else r <- tclvalue(t)
        })
        if (!is.null(classes)) 
            r.class <<- as.numeric(tkcurselection(t.class)) + 1
        if ((length(r.class) == 0) && !(length(init.class) == 0)) 
            r.class <<- init.class + 1
        tkgrab.release(dlg)
        tkdestroy(dlg)
        if (!is.null(top)) 
            tkfocus(top)
    }
    dlg <- tktoplevel()
    tkwm.deiconify(dlg)
    if (do.grab) 
        tkgrab.set(dlg)
    tkfocus(dlg)
    tkwm.title(dlg, paste(title, ": ", name.object))
    idCancel <- NULL
    r.class <- idCancel
    result <- idCancel
    tkgrid(tklabel(dlg, text = "    "))
    tkgrid(tklabel(dlg, text = "Name: "), tklabel(dlg, text = name.object))
    tkgrid(tklabel(dlg, text = "    "))
    t.list <- lapply(slotNames(object), function(slot.name) {
        slot.value <- slot(object, slot.name)
        if (is.language(slot.value)) 
            size <- 0
        else size <- length(c(unlist(slot.value)))
        color <- "white"
        if (slot.name %in% difficultSlots) 
            color <- "grey"
        if (slot.name %in% fixedSlots) 
            color <- "DarkGrey"
        if (is.language(slot.value)) 
            t <- subTextDialog(dlg, subtitle = slot.name, 
                entryInit = paste(slot.value), 
                entryWidth = entryWidth, background = color)
        else if (size > 1) 
            t <- subListDialog(dlg, subtitle = slot.name, 
                               entryInit = slot.value, 
                               entryWidth = entryWidth, background = color)
        else t <- subTextDialog(dlg, subtitle = slot.name, 
                                entryInit = slot.value, 
            entryWidth = entryWidth, background = color)
        return(t)
    })
    init.class <- list()
    if (!is.null(classes)) {
        init.class <- which(class(object) == classes[, 2]) - 1
        if (length(init.class) == 0) 
            t.class <- subSelectDialog(dlg, class(object), init = NULL, 
                title = title, subtitle = "Class:", background = "grey")
        else {
            color <- "white"
            if ("class" %in% difficultSlots) 
                color <- "grey"
            if ("class" %in% fixedSlots) 
                color <- "DarkGrey"
            t.class <- subSelectDialog(dlg, classes[, 1], init = init.class, 
                title = title, subtitle = "Class:", background = color)
        }
    }
    OK.but <- tkbutton(dlg, text = "   OK   ", command = onOK)
    Cancel.but <- tkbutton(dlg, text = " Cancel ", command = onCancel)
    tkgrid(tklabel(dlg, text = "    "))
    if (TRUE || okReturn) 
        tkgrid(tklabel(dlg, text = "    "), OK.but, Cancel.but, 
            tklabel(dlg, text = "    "))
    else tkgrid(tklabel(dlg, text = "    "), Cancel.but, tklabel(dlg, 
        text = "    "))
    tkgrid(tklabel(dlg, text = "    "))
    tkfocus(dlg)
    tkbind(dlg, "<Destroy>", function() {
        tkgrab.release(dlg)
        if (!is.null(top)) 
            tkfocus(top)
    })
    tkwait.window(dlg)
    Result <- list()
    "isCancel" <- function(x) (is.null(x))
    old.slot.names <- slotNames(object)
    if (!is.null(classes) && !isCancel(r.class) && !(length(init.class) == 0) 
        && (init.class != r.class - 1)) {
        if ("class" %in% difficultSlots) 
            message(paste("Trying to change the difficult slot '", 
                "class", "' ; "))
        if ("class" %in% fixedSlots) 
            message(paste("Trying to change the fixed slot '", 
                "class", "' ; "))
        class(object) <- classes[r.class, 2]
        Result <- list(class = classes[r.class, 2])
    }
    if (!is.null(result) && !isCancel(result)) {
        new.slot.names <- slotNames(object)
        for (i in 1:length(result)) {
            slot.name <- old.slot.names[i]
            slot.value <- slot(object, slot.name)
            slot.class <- class(slot.value)
            r <- result[[i]]
            if (is.language(slot.value)) 
                size <- 0
            else size <- length(c(unlist(slot.value)))
            if (is.language(slot.value)) 
                r <- slot.value
            else if (size == 1) 
                class(r) <- class(slot.value)
            else if (slot.class == "matrix") {
                r <- as.numeric(r)
                dim(r) <- dim(slot.value)
            }
            else if (slot.class == "numeric") {
                r <- as.numeric(r)
            }
            else r <- r
            if (is.numeric(slot.value)) {
                if (any(is.na(r))) {
                  message("Invalid number;")
                  change <- FALSE
                }
                else {
                  diff <- r - slot.value
                  diff[r != 0] <- diff[r != 0]/r[r != 0]
                  change <- (sum(diff^2) > 0.001)
                }
            }
            else if (is.logical(slot.value)) {
                change <- !(all(r == slot.value))
                if (any(is.na(c(r, slot.value)))) 
                  change <- !(all(is.na(r) == is.na(slot.value)))
            }
            else change <- !(all(r == slot.value))
            if (change && is.element(slot.name, new.slot.names)) {
                if (slot.name %in% difficultSlots) 
                  message(paste("Trying to change the difficult slot '", 
                    slot.name, "' ; "))
                if (slot.name %in% fixedSlots) 
                  message(paste("Trying to change the fixed slot '", 
                    slot.name, "' ; "))
                slotname <- slot.name
                if (slot.name == "label.position") 
                  slotname <- "labelPosition"
                if (slot.name == "vertex.indices") 
                  slotname <- "nodeIndices"
                r <- slot(eval(call(paste(slotname, "<-", sep = ""), 
                  object, r)), slot.name)
                slot(object, slot.name) <- r
                lo <- slot(object, slot.name)
                item <- list(lo)
                names(item) <- slot.name
                Result <- append(Result, item)
            }
        }
    }
    return(list(object = object, values = Result))
}

".removeNull" <-
function (x) 
{
    result <- NULL
    for (i in seq(along = x)) if (!is.null(x[[i]])) 
        result <- append(result, list(x[[i]]))
    class(result) <- class(x)
    return(result)
}

".StrBlockTree" <-
function (tree, title = "blockTree") 
{
    message(title)
    "subStrBlockTree" <- function(tree, k, p, dept) {
        i <- abs(tree$block@index)
        s <- tree$block@stratum
        label <- tree$block@label
        d <- paste(tree$block@descendants, collapse = ",")
        a <- paste(tree$block@ancestors, collapse = ",")
        nc <- 10 - nchar(a)
        if (nc < 0) 
            nc <- 1
        cat(paste("|", paste(rep(".", dept), collapse = ""), 
            i, paste(rep(" ", 10 - dept), collapse = ""), "(", 
            k, "/", p, ")", ": ", label, "; Stratum: ", s, "; Ancestors: ", 
            a, paste(rep(" ", nc), collapse = ""), "; Descendants: ", 
            d, "\n", sep = ""))
        if (!is.null((tree$sub.blocks))) 
            for (j in 1:length(tree$sub.blocks)) {
                subStrBlockTree(tree$sub.blocks[[j]], j, i, dept + 
                  1)
            }
    }
    if (!.IsEmpty(tree) && (length(tree) > 0)) 
        subStrBlockTree(tree, 1, 0, 0)
    invisible()
}

".visit.envs" <-
function (X = ls(.DgRoot$env, all.names = all.names), all.names = TRUE) 
{
    print(X)
    for (x in X[X != "num.subenv"]) {
        env.x <- get(x, .DgRoot$env)
        if (class(env.x) == "dgenv") {
            print(x)
            print(ls(env.x$env, all.names = all.names))
            Y <- ls(env.x$env, all.names = all.names)
            for (y in Y[(Y != "num.subenv") & (Y != "parent")]) {
                env.y <- get(y, env.x$env)
                if (class(env.y) == "dgenv") {
                  print(y)
                  print(ls(env.y$env, all.names = all.names))
                  Z <- ls(env.y$env, all.names = all.names)
                  for (z in Z[(Z != "num.subenv") & (Z != "parent")]) {
                    env.z <- get(z, env.y$env)
                    if (class(env.z) == "dgenv") {
                      print(z)
                      print(ls(env.z$env, all.names = all.names))
                    }
                  }
                }
            }
        }
    }
}
