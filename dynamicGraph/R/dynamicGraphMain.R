"dynamicGraphMain" <-
function (vertexList = NULL, blockList = NULL, dg = NULL, object = NULL, 
    objectName = NULL, control = dg.control(...), ...) 
{
    if ((is.null(dg))) 
        dg <- .newDgGraphEdges(vertexList = vertexList, 
                               blockList = blockList, ...)
    args <- list(...)
    if (!is.null(args$title) && (control$label == "dynamicGraph")) {
        control$label <- args$title
    }
    blockTree <- NULL
    if (length(dg@visibleVertices) == 0) 
        dg@visibleVertices <- 1:length(vertexList)
    if (length(dg@visibleBlocks) == 0) 
        dg@visibleBlocks <- 1:length(blockList)
    if (length(blockList) == 0) 
        blockList <- NULL
    tkselectionForVisibleVertices <- FALSE
    permit.update.block.index <- TRUE
    font.vertex.label <- "8x16"
    font.edge.label <- "8x16"
    font.block <- "10x20"
    min.x <- -control$margin
    max.x <- control$width + control$margin
    d.x <- -min.x/(max.x - min.x)
    min.y <- -control$margin
    max.y <- control$height + control$margin
    d.y <- -min.y/(max.y - min.y)
    initial.set.popups <- FALSE
    colors <- c("DarkGreen", "navy", "NavyBlue", "DarkBlue", 
        "DarkRed", "MidnightBlue", "DarkSlateGray", "DarkSlateGrey", 
        "MediumBlue", "ForestGreen", "SaddleBrown", "DarkOliveGreen", 
        "firebrick", "brown", "blue", "green", "red", "DarkSlateBlue", 
        "SeaGreen", "DarkCyan", "DarkMagenta", "OliveDrab", "sienna", 
        "LimeGreen", "DimGray", "DimGrey", "maroon", "OrangeRed", 
        "DarkGoldenrod", "chocolate", "MediumSeaGreen", "DarkViolet", 
        "LawnGreen", "LightSeaGreen", "SteelBlue", "chartreuse", 
        "SpringGreen", "black", "SlateGray", "SlateGrey", "VioletRed", 
        "IndianRed", "DarkOrange", "RoyalBlue", "peru", "SlateBlue", 
        "BlueViolet", "DarkOrchid", "LightSlateGray", "LightSlateGrey", 
        "YellowGreen", "CadetBlue", "DarkTurquoise", "goldenrod", 
        "orange", "DeepPink", "tomato", "DodgerBlue", "purple", 
        "DeepSkyBlue", "coral", "gold", "DarkSeaGreen", "RosyBrown", 
        "GreenYellow", "MediumPurple", "PaleVioletRed", "DarkKhaki", 
        "MediumOrchid", "CornflowerBlue", "salmon", "LightCoral", 
        "turquoise", "LightSlateBlue", "SandyBrown", "DarkSalmon", 
        "DarkGray", "DarkGrey", "cyan", "magenta", "yellow", 
        "LightGreen", "tan", "LightSalmon", "HotPink", "burlywood", 
        "orchid", "PaleGreen", "gray", "grey", "SkyBlue", "LightGoldenrod", 
        "LightSkyBlue", "aquamarine", "LightSteelBlue", "plum", 
        "violet", "khaki", "LightBlue", "thistle", "LightPink", 
        "PowderBlue", "LightGray", "LightGrey", "PaleGoldenrod", 
        "wheat", "NavajoWhite", "pink", "PaleTurquoise", "PeachPuff", 
        "gainsboro", "moccasin", "bisque", "BlanchedAlmond", 
        "AntiqueWhite", "PapayaWhip", "MistyRose", "beige", "lavender", 
        "LemonChiffon", "linen", "cornsilk", "OldLace", "LightCyan", 
        "LightYellow", "honeydew", "WhiteSmoke", "seashell", 
        "LavenderBlush", "AliceBlue", "FloralWhite", "azure", 
        "ivory", "MintCream", "GhostWhite", "snow", "white")
    "myColor" <- function(i) colors[min(i%%137, length(colors))]
    "drawModel" <- function(frameModels = NULL, frameViews = NULL, 
        graphWindow = NULL, dg = NULL, object = NULL, 
        frameModelsEnv = .get.env.frameModels(frameModels = frameModels), 
        initialWindow = FALSE, returnNewMaster = FALSE, redraw = FALSE, 
        returnFrameModel = TRUE, control = NULL, ...) {
        if ((is.null(dg))) 
            dg <- .newDgGraphEdges(vertexList = vertexList, 
                                   blockList = blockList, ...)
        args <- list(...)
        if (!is.null(args$title)) {
            control$label <- args$title
        }
        "redrawView" <- function(frameModels = NULL, frameViews = NULL, 
            graphWindow = NULL, dg = NULL, initialWindow = FALSE, 
            returnNewMaster = FALSE, redraw = FALSE, returnFrameModel = TRUE, 
            setUpdateCountModelMain = FALSE, control = NULL, ...) {
            args <- list(...)
            if (!is.null(args$title)) {
                control$label <- args$title
            }
            "setSrcLabel" <- function(viewLabel) {
                "g" <- function(x, i, max) {
                  if (!is.null(zoomCenter)) 
                    x <- (x - c(zoomCenter[i]))/Scale + zoomCenter[i]
                  x <- c(x[1], mean(x), x[2])
                  x <- c(x, 100/max * x - 50)
                  return(round(x))
                }
                x <- (min.x + (max.x - min.x) * as.numeric(tkxview(canvas)))
                x <- g(x, 1, control$width)
                x.txt <- paste(format(x, digits = 3, trim = FALSE), 
                  collapse = ", ")
                y <- (min.y + (max.y - min.y) * as.numeric(tkyview(canvas)))
                y <- g(y, 2, control$height)
                y.txt <- paste(format(y, digits = 3, trim = FALSE), 
                  collapse = ", ")
                ViewType <- control$viewClasses[control$viewClasses[, 
                  2] == class(GraphWindow), 1]
                tkconfigure(viewLabel, text = paste(ViewType, 
                  " | ", "X:", x.txt, "Y:", y.txt))
            }
            "getLabel" <- function() {
                title <- paste("Model: ", dm.frameViews@index, 
                  sep = " ")
                if (exists("GraphWindow")) 
                  title <- paste(title, "; Graph: ", GraphWindow@index, 
                    GraphWindow@id, sep = " ")
                return(title)
            }
            "tkinsert.blockVertices" <- function(parent, stratum, 
                treeWidget, delete = FALSE) {
                for (j in seq(along = vertexList)) {
                  vertex <- vertexList[[j]]
                  vertex.stratum <- retStratum(j, vertex.type = "Vertex")
                  if (vertex.stratum == stratum) {
                    child <- name(vertex)
                    fill <- "DarkGreen"
                    if (!is.element(j, dg@visibleVertices)) {
                      child <- tdv(child)
                      fill <- "SpringGreen"
                    }
                    if (delete) 
                      tkdelete(treeWidget, child)
                    else tkinsert(treeWidget, "end", parent, 
                      child, text = child, fill = fill)
                  }
                }
            }
            "ubl" <- function(label = block@label, index = block@index, 
                block = NULL) paste(c("", abs(index), label), 
                collapse = ":")
            "iubl" <- function(s, n = nchar(s), m = 1 + (1:n)[substring(s, 
                1:n, 1:n) == ":"][2]) paste(substring(s, m:n, 
                m:n), collapse = "")
            "tdv" <- function(name) paste(c("(", name, ")"), collapse = "")
            "tkinsert.block" <- function(treeWidget, parent, 
                block, open = openTreeBlock[abs(block@index)], 
                delete = FALSE, fill = "Blue") {
                Child <- block@label
                child <- ubl(Child, block = block)
                if (delete) 
                  tkdelete(treeWidget, child)
                else tkinsert(treeWidget, "end", parent, child, 
                  open = open, text = Child, fill = fill)
            }
            "tkinsert.blockTree" <- function(treeWidget, tree, 
                delete = FALSE) {
                "subVisitBlockTree" <- function(tree) {
                  parent <- ubl(block = tree$block)
                  tkinsert.blockVertices(parent, tree$block@stratum, 
                    treeWidget, delete)
                  if (!is.null((tree$sub.blocks))) 
                    for (j in 1:length(tree$sub.blocks)) {
                      tkinsert.block(treeWidget, parent, 
                                     block = tree$sub.blocks[[j]]$block, 
                                     delete)
                      subVisitBlockTree(tree$sub.blocks[[j]])
                    }
                }
                if (!.IsEmpty(blockList) && !.IsEmpty(tree) && 
                  !(length(tree) == 0)) {
                  tkinsert.blockVertices(parent = "root", 0, 
                    treeWidget, delete)
                  tkinsert.block(treeWidget, parent = "root", 
                    block = tree$block, open = TRUE, delete = delete, 
                    fill = "DarkBlue")
                  subVisitBlockTree(tree)
                }
            }
            "oldtkinsert.blockList" <- function(treeWidget, blockList, 
                delete = FALSE) {
                if (!.IsEmpty(blockList)) 
                  for (j in 1:length(blockList)) {
                    parent <- "root"
                    ancestor <- max(ancestors(blockList[[j]]))
                    if (ancestor > 0) 
                      parent <- ubl(block = blockList[[ancestor]])
                    tkinsert.block(treeWidget, parent = parent, 
                      block = blockList[[j]], open = openTreeBlock[j], 
                      delete = delete)
                    parent <- ubl(block = blockList[[j]])
                    tkinsert.blockVertices(parent, blockList[[j]]@stratum, 
                      treeWidget, delete)
                  }
                tkinsert.blockVertices(parent = "root", 0, treeWidget, 
                  delete)
            }
            "tkinsert.blockList" <- function(treeWidget, blockList, 
                delete = FALSE) {
                tkinsert.blockVertices("root", stratum = 0, treeWidget, 
                  delete = delete)
                if (!.IsEmpty(blockList)) 
                  for (j in 1:length(blockList)) {
                    if (0 == parent(blockList[[j]])) {
                      tkinsert.subblockList(treeWidget, blockList, 
                        index = j, parent = 0, 
                        stratum = stratum(blockList[[j]]), 
                        delete = delete)
                    }
                  }
            }
            "tkinsert.subblockList" <- function(treeWidget, blockList, 
                index, parent = blockList[[index]]@parent, 
                stratum = blockList[[index]]@stratum, delete = FALSE) {
                parentName <- "root"
                if (parent > 0) 
                  parentName <- ubl(block = blockList[[parent]])
                tkinsert.block(treeWidget, parent = parentName, 
                  block = blockList[[index]], open = openTreeBlock[index], 
                  delete = delete)
                parentName <- ubl(block = blockList[[index]])
                tkinsert.blockVertices(parentName, stratum, treeWidget, 
                  delete)
                if (!.IsEmpty(blockList)) 
                  for (j in 1:length(blockList)) {
                    if (index == parent(blockList[[j]])) {
                      tkinsert.subblockList(treeWidget, blockList, 
                        index = j, parent = index, 
                        stratum = stratum(blockList[[j]]), 
                        delete = delete)
                    }
                  }
            }
            "grepVertices" <- function(name) {
                selectedVertices <- NULL
                for (j in 1:length(vertexList)) 
                  if (length(grep(namesVertices[j], name)) > 0) 
                  selectedVertices <- c(selectedVertices, j)
                return(selectedVertices)
            }
            "grepBlocks" <- function(name) {
                selectedBlocks <- NULL
                if (!.IsEmpty(blockList)) 
                  for (j in 1:length(blockList)) {
                    blockname <- ubl(block = blockList[[j]])
                    if (length(grep(blockname, name)) > 0) 
                      selectedBlocks <- c(selectedBlocks, j)
                  }
                return(selectedBlocks)
            }
            "popupSelectedInPanel" <- function() {
                selectedVertices <- NULL
                selectedBlocks <- NULL
                if (control$variableFrame) {
                  box <- GW.top$env$box
                  if ((get("type", box$env) == "variableList")) {
                    vv <- returnVisibleVertices()
                    sv <- as.numeric(tkcurselection(box)) + 1
                    if (tkselectionForVisibleVertices) 
                      selectedVertices <- c(setdiff(vv, sv), 
                        setdiff(sv, vv))
                    else selectedVertices <- sv
                  }
                  else {
                    sv <- tclvalue(tcl(box, "selection", "get"))
                    selectedVertices <- c(selectedVertices, grepVertices(sv))
                    selectedBlocks <- c(selectedBlocks, grepBlocks(sv))
                  }
                }
                for (i in selectedVertices) {
                  vertex.type <- "Vertex"
                  if (closedVertex[i]) 
                    vertex.type <- "Vertex"
                  callPopupInBox(vertex = vertexList[[i]], i = i, 
                    vertex.type = vertex.type, 
                    UserNodePopupItems = control$UserMenus)()
                }
                for (i in selectedBlocks) {
                  vertex.type <- "OpenBlock"
                  if (closedBlock[i]) 
                    vertex.type <- "ClosedBlock"
                  callPopupInBox(vertex = blockList[[i]], i = i, 
                    vertex.type = vertex.type, 
                    UserNodePopupItems = control$UserMenus)()
                }
            }
            "bindBox" <- function(box, label = "newGraph") {
                if (control$debug.position) 
                  print(label)
                tcl("bind", box, "<F11>", function(...) {
                  print("<<<F11>>>")
                })
                tkbind(box, "<F12>", function(...) {
                  print("<<<F12>>>")
                })
                tkbind(box, "<Button-1>", function() {
                  if (control$debug.position) 
                    print("<Button-1>")
                  vv <- returnVisibleVertices()
                  if (!(get("type", box$env) == "variableList")) {
                    sv <- tclvalue(tcl(box, "selection", "get"))
                  }
                  else {
                    sv <- as.numeric(tkcurselection(box)) + 1
                  }
                  if (control$debug.position) 
                    print(sv)
                  if (control$debug.position) 
                    print(vv)
                })
                if (control$debug.position) 
                  print("<Button-1>")
                tkbind(box, "<Double-Button-1>", function() {
                  if (control$debug.position) 
                    print("<Double-Button-1>")
                  if (tkselectionForVisibleVertices) {
                    vv <- returnVisibleVertices()
                    sv <- as.numeric(tkcurselection(box)) + 1
                    for (i in setdiff(vv, sv)) subDropVertex(i, 
                      slave = FALSE, upd = FALSE)
                    for (i in setdiff(sv, vv)) subAddVertex(i, 
                      slave = FALSE)
                  }
                  else {
                    for (i in (as.numeric(tkcurselection(box)) + 
                      1)) {
                      if (is.element(i, returnVisibleVertices())) 
                        subDropVertex(i, slave = FALSE, upd = FALSE)
                      else subAddVertex(i, slave = FALSE)
                    }
                  }
                })
                if (control$debug.position) 
                  print("<Double-Button-1>")
                tkbind(box, "<Button-2>", function() {
                  if (control$debug.position) 
                    print("<Button-2>")
                  if (!(get("type", box$env) == "variableList")) {
                    sv <- tclvalue(tcl(box, "selection", "get"))
                    if (control$debug.position) 
                      print(sv)
                  }
                  else {
                    vv <- returnVisibleVertices()
                    sv <- as.numeric(tkcurselection(box)) + 1
                    for (i in setdiff(vv, sv)) propertyNode(i, 
                      vertex.type = "Vertex")()
                    for (i in setdiff(sv, vv)) propertyNode(i, 
                      vertex.type = "Vertex")()
                  }
                })
                if (control$debug.position) 
                  print("<Button-2>")
                tkbind(box, "<Enter>", function() {
                  if (control$debug.position) 
                    print("<Enter>")
                  subUpdatePositions(label)
                  if ((get("type", box$env) == "variableList")) {
                    if (tkselectionForVisibleVertices) {
                      for (i in 1:length(vertexList)) {
                        tkselection.clear(box, i - 1)
                      }
                      for (i in returnVisibleVertices()) {
                        tkselection.set(box, i - 1)
                      }
                    }
                  }
                  else {
                    sv <- tclvalue(tcl(box, "selection", "get"))
                    if (control$debug.position) 
                      print(sv)
                  }
                })
                if (control$debug.position) 
                  print("<Enter>")
                tkbind(box, "<Leave>", function() {
                  if (control$debug.position) 
                    print("<Leave>")
                  subUpdatePositions(label)
                })
                if (control$debug.position) 
                  print("<Leave>")
                tkbind(box, "<ButtonRelease-1>", function() {
                  if (control$debug.position) 
                    print("<ButtonRelease-1>")
                  subUpdatePositions(label)
                  if ((get("type", box$env) == "variableList")) {
                  }
                  else {
                    sv <- tclvalue(tcl(box, "selection", "get"))
                    if (control$debug.position) 
                      print(sv)
                  }
                })
                if (control$debug.position) 
                  print("<ButtonRelease-1>")
                tkbind(box, "<Motion>", function() {
                })
                if (control$debug.position) 
                  print("<Motion>")
                tkbind(box, "<Button-3>", function() {
                  print("<Button-3....>")
                  if (control$debug.position) 
                    print("<Button-3>")
                  popupSelectedInPanel()
                })
                if (control$debug.position) 
                  print("<Button-3>")
            }

            "moveCanvas" <- function() {
                o.x <- NULL
                o.y <- NULL
                function(x, y) {
                  n.x <- as.numeric(x) * d.x / min.x
                  n.y <- as.numeric(y) * d.y / min.y
                  if (!is.null(o.y)) { 
                    d.x <- n.x - o.x
                    d.y <- n.y - o.y
                    if ((abs(d.x) < 0.1) && (abs(d.y) < 0.1)) {
                      r.x <- as.numeric(tkxview(GW.top$env$canvas))[1]
                      r.y <- as.numeric(tkyview(GW.top$env$canvas))[1]
                      tkxview.moveto(GW.top$env$canvas, r.x + d.x)
                      tkyview.moveto(GW.top$env$canvas, r.y + d.y)
                      setSrcLabel(GW.top$env$viewLabel)
                    }
                  }
                  o.x <<- n.x
                  o.y <<- n.y
                }
            }

            "newGraph" <- function(viewType, ldg, title = "Graph diddler", 
                index = -1, id = 0, close.enough = control$closeenough, 
                parent = "", background = "white", width = 400, 
                height = 400) {
                prototype <- "DynamicGraphView"
                x <- match(viewType, control$viewClasses[, 1])
                if (is.null(x) || all(is.na(x))) {
                  x <- match(viewType, control$viewClasses[, 
                    2])
                  viewType <- paste(control$viewClasses[, 1][x])
                }
                if (!is.null(x) && !all(is.na(x))) 
                  prototype <- paste(control$viewClasses[, 2][x])
                top <- tktoplevel()
                f0.box <- NULL
                viewLabel <- NULL
                tags <- NULL
                if (control$permitZoom) {
                  f <- tkframe(top)
                  tkpack(f, expand = "yes", side = "top", fill = "both")
                  if (control$variableFrame) {
                    w.pane <- FALSE
                    try(w.pane <- tcl("panedwindow", .Tk.subwin(f)), 
                      silent = TRUE)
                    has.paned <- (class(w.pane) != "logical")
                    if (!has.paned) {
                      control$variableFrame <<- has.paned
                      message("'panedwindow' not in your version of 'Tcl/tk'")
                    }
                  }
                  if (control$variableFrame) {
                    tkpack(w.pane, side = "top", expand = "yes", 
                      fill = "both", pady = 2, padx = "2m")
                    tcl("wm", "iconname", top, "Paned window")
                    "myTclRequire" <- function(package, warn = TRUE) {
                      a <- tclvalue(.Tcl(paste("package versions ", 
                        package)))
                      if (length(a) == 1 && nchar(a) == 0 && a == "") {
                        if (warn) 
                          warning(paste("Tcl package", package, 
                            "not found."))
                        return(FALSE)
                      }
                      else .Tcl(paste("package require ", package))
                    }
                    f0 <- tkframe(f)
                    if ((!.IsEmpty(blockList) || !.IsEmpty(blockTree)) && 
                      (class(myTclRequire("BWidget", warn = FALSE)) == 
                        "tclObj")) {
                      "selectcommand" <- function(...) {
                        sv <- NULL
                        if ((get("type", GW.top$env$box$env) == 
                          "variableList")) {
                          sv <- as.numeric(tkcurselection(box)) + 
                            1
                        }
                        else {
                          sv <- tclvalue(tcl(box, "selection", 
                            "get"))
                        }
                        print(paste(" Node: ", list(...)[[1]], 
                          "; ", list(...)[[2]], "; Cut: ", sv))
                      }
                      "dropcmd" <- function(...) {
                        moveVertex <- grepVertices(list(...)[[6]])
                        moveBlock <- grepBlocks(list(...)[[6]])
                        toVertex <- grepVertices(list(...)[[3]])
                        toBlock <- grepBlocks(list(...)[[3]])
                        toBlockIndex <- NULL
                        if (!is.null(toBlock)) {
                          if (length(toBlock) > 1) {
                            message(paste("Ignoring blocks: ", 
                              paste(toBlock[2:length(toBlock)], 
                                collapse = ",")))
                            toBlock <- toBlock[1]
                          }
                          if (!is.null(toVertex)) {
                            message(paste("Ignoring vertices: ", 
                              paste(toVertex, collapse = ",")))
                          }
                          toBlockIndex <- toBlock
                        }
                        else if (!is.null(toVertex)) {
                          if (length(toVertex) > 1) {
                            message(paste("Ignoring vertices: ", 
                              paste(toVertex[2:length(toVertex)], 
                                collapse = ",")))
                            toVertex <- toVertex[1]
                          }
                          toBlockIndex <- retBlockIndex(toVertex)
                        }
                        if (!is.null(moveVertex)) {
                          if (length(moveVertex) > 1) {
                            message(paste("Ignoring vertices: ", 
                              paste(moveVertex[2:length(moveVertex)], 
                                collapse = ",")))
                            moveVertex <- moveVertex[1]
                          }
                          if (!is.null(moveBlock)) 
                            message(paste("Ignoring blocks with same name: ", 
                              paste(moveBlock, collapse = ",")))
                          posFrom <- retVertexPos(moveVertex, 
                            vertex.type = "Vertex")
                          if (!is.null(toBlockIndex)) {
                            posTo <- retVertexPos(toBlockIndex, 
                              vertex.type = "ClosedBlock")
                            dxy = posTo - posFrom
                            if (setVertexPos(moveVertex, posTo, 
                              dxy, vertex.type = "Vertex")) {
                              if (toBlockIndex == retBlockIndex(moveVertex)) 
                                setUpdateBlockEdges("dropcmd")
                              if (closedVertex[moveVertex] && 
                                !closedBlock[toBlockIndex]) 
                                setCloseVertex(moveVertex, 
                                               !closedVertex[moveVertex], 
                                  vertex.type = "Vertex")
                              if (!closedVertex[moveVertex]) 
                                tkmove(canvas, vertexItem(moveVertex)$tag, 
                                  dxy[1], dxy[2])
                            }
                          }
                        }
                        else if (!is.null(moveBlock)) {
                          if (length(moveBlock) > 1) {
                            message(paste("Ignoring blocks: ", 
                              paste(moveBlock[2:length(moveBlock)], 
                                collapse = ",")))
                            moveBlock <- moveBlock[1]
                          }
                          old.parentName <- "root"
                          old.parent <- max(ancestors(blockList[[moveBlock]]))
                          if (old.parent > 0) 
                            old.parentName <- ubl(block = blockList[[old.parent]])
                          if (is.null(toBlockIndex)) 
                            toBlockIndex <- 0
                          moveBlockName <- ubl(block = blockList[[moveBlock]])
                          if (toBlockIndex == 0) 
                            toBlockName <- "root"
                          else toBlockName <- ubl(block = blockList[[toBlockIndex]])
                          if (moveBlock == toBlockIndex) 
                            message("No move!")
                          else if (!(toBlockIndex == 0) && is.element(toBlockIndex, 
                            blockList[[moveBlock]]@descendants)) 
                            message("Move to child!")
                          else {
                            setUpdatePositions("dropcmd")
                            tkinsert.subblockList(GW.top$env$box, 
                              blockList, index = moveBlock, delete = TRUE)
                            parentFromBlock <- blockList[[moveBlock]]@parent
                            if (parentFromBlock > 0) {
                              x <- blockList[[parentFromBlock]]@children
                              blockList[[parentFromBlock]]@children <<- x[x != 
                                moveBlock]
                            }
                            blockList[[moveBlock]]@parent <<- toBlockIndex
                            if (toBlockIndex > 0) {
                              x <- unique(c(blockList[[toBlockIndex]]@children, 
                                moveBlock))
                              blockList[[toBlockIndex]]@children <- x[x != 0]
                            }
                            x <- unique(c(moveBlock, blockList[[moveBlock]]@descendants))
                            x <- x[x != 0]
                            for (j in blockList[[moveBlock]]@ancestors) {
                              if (j > 0) {
                                blockList[[j]]@descendants <<- setdiff(blockList[[j]]@descendants, 
                                  x)
                              }
                            }
                            if (toBlockIndex > 0) {
                              for (j in blockList[[toBlockIndex]]@ancestors) {
                                if (j > 0) {
                                  z <- unique(c(blockList[[j]]@descendants, x))
                                  z <- z[z != 0]
                                  blockList[[j]]@descendants <<- z
                                }
                              }
                              z <- unique(c(blockList[[toBlockIndex]]@descendants, 
                                x))
                              z <- z[z != 0]
                              blockList[[toBlockIndex]]@descendants <<- z
                            }
                            if (toBlockIndex > -1) {
                              oldAnc <- blockList[[moveBlock]]@ancestors
                              if (toBlockIndex > 0) 
                                z <- unique(c(blockList[[toBlockIndex]]@ancestors, 
                                  toBlockIndex))
                              else z <- 0
                              z <- z[z != 0]
                              blockList[[moveBlock]]@ancestors <<- z
                              if (toBlockIndex > 0) 
                                newAnc <- unique(c(toBlockIndex, 
                                  blockList[[toBlockIndex]]@ancestors))
                              else newAnc <- 0
                              newAnc <- newAnc[newAnc != 0]
                              for (j in blockList[[moveBlock]]@descendants) {
                                if (j > 0) {
                                  z <- unique(c(newAnc, setdiff(blockList[[j]]@ancestors, 
                                    oldAnc)))
                                  z <- z[z != 0]
                                  blockList[[j]]@ancestors <<- z
                                }
                              }
                            }
                            blockList <<- checkBlockList(blockList)
                            tkinsert.subblockList(GW.top$env$box, 
                              blockList, index = moveBlock, delete = FALSE)
                            blockTree <<- list(NULL)
                          }
                        }
                      }
                      "opencmd" <- function(...) {
                        openTreeBlock[which(iubl(list(...)) == 
                          Names(blockList))] <<- TRUE
                      }
                      "closecmd" <- function(...) {
                        openTreeBlock[which(iubl(list(...)) == 
                          Names(blockList))] <<- FALSE
                      }
                      "dropovercmd" <- function(...) {
                        print(paste("DropOverCmd: ", paste(unlist(list(...)), 
                          collapse = ";")))
                      }
                      "draginitcmd" <- function(...) {
                        print(paste("DragInitCmd: ", paste(unlist(list(...)), 
                          collapse = ";")))
                      }
                      "dragendcmd" <- function(...) {
                        print(paste("DragEndCmd: ", paste(unlist(list(...)), 
                          collapse = ";")))
                      }
                      "f0.box" <- tkwidget(f0, "Tree", relief = "raised", 
                        background = background, highlightcolor = "Black", 
                        selectbackground = "gray", selectforeground = "white", 
                        dragenabled = TRUE, dragevent = 1, dropenabled = TRUE, 
                        dropcmd = function(...) dropcmd(...), 
                        opencmd = function(...) opencmd(...), 
                        closecmd = function(...) closecmd(...), 
                        selectcommand = function(...) selectcommand(...), 
                        yscrollcommand = function(...) tkset(f0.scr, 
                          ...))
                      if (FALSE && !.IsEmpty(blockTree)) {
                        assign("type", "blockTree", envir = f0.box$env)
                        tkinsert.blockTree(f0.box, blockTree)
                      }
                      else {
                        assign("type", "blockList", envir = f0.box$env)
                        tkinsert.blockList(f0.box, blockList, 
                          delete = FALSE)
                      }
                    }
                    else {
                      f0.box <- tklistbox(f0, relief = "raised", 
                        background = background, selectmode = "multiple", 
                        highlightcolor = "Black", foreground = "DarkGreen", 
                        yscrollcommand = function(...) tkset(f0.scr, 
                          ...))
                      assign("type", "variableList", envir = f0.box$env)
                      for (i in (1:length(namesVertices))) {
                        child <- namesVertices[i]
                        if (!is.element(i, ldg@visibleVertices)) 
                          child <- tdv(child)
                        tkinsert(f0.box, "end", child)
                      }
                      if (tkselectionForVisibleVertices) 
                        for (i in ldg@visibleVertices) {
                          tkselection.set(f0.box, i - 1)
                        }
                    }
                    tkpack(f0.box, expand = "yes", fill = "both")
                    f0.scr <- tkscrollbar(f0, repeatinterval = 5, 
                      orient = "vertical", command = function(...) {
                        tkyview(f0.box, ...)
                      })
                    tkgrid(f0.box, padx = 1, pady = 1, row = 0, 
                      column = 1, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(f0.scr, padx = 1, pady = 1, row = 0, 
                      column = 2, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid.columnconfigure(f0, 1, weight = 1, 
                      minsize = 0)
                    tkgrid.rowconfigure(f0, 0, weight = 1, minsize = 0)
                    bindBox(f0.box)
                    f1 <- tkframe(f)
                    viewLabel <- tklabel(f1, text = viewType, 
                      foreground = "DarkSlateBlue", background = "LightGrey")
                    xscr <- tkscrollbar(f1, repeatinterval = 5, 
                      orient = "horizontal", background = "white", 
                      command = function(...) {
                        tkxview(canvas, ...)
                        setSrcLabel(viewLabel)
                      })
                    yscr <- tkscrollbar(f1, repeatinterval = 5, 
                      orient = "vertical", command = function(...) {
                        tkyview(canvas, ...)
                        setSrcLabel(viewLabel)
                      })
                    canvas <- tkcanvas(f1, relief = "raised", 
                      background = background, closeenough = close.enough, 
                      borderwidth = 5, 
                      scrollregion = c(min.x, min.y, max.x, max.y), 
                      xscrollincrement = -4, 
                      yscrollincrement = -4, 
                      xscrollcommand = function(...) tkset(xscr, ...), 
                      yscrollcommand = function(...) tkset(yscr, ...), 
                      width = width, height = height)
                    tkpack(canvas, expand = "yes", fill = "both")
                    tkgrid(canvas, padx = 1, pady = 1, row = 0, 
                      column = 1, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(yscr, padx = 1, pady = 1, row = 0, 
                      column = 2, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(xscr, padx = 1, pady = 1, row = 1, 
                      column = 1, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(viewLabel, padx = 1, pady = 1, row = 2, 
                      column = 1, rowspan = 1, columnspan = 2, 
                      sticky = "news")
                    tkgrid.columnconfigure(f1, 1, weight = 1, 
                      minsize = 0)
                    tkgrid.rowconfigure(f1, 0, weight = 1, minsize = 0)
                    tkadd(w.pane, f0, f1)
                    tkfocus(canvas)
                    tkfocus(top)
                  }
                  else {
                    viewLabel <- tklabel(f, text = viewType, 
                      foreground = "DarkSlateBlue", background = "LightGrey")
                    xscr <- tkscrollbar(f, repeatinterval = 5, 
                      orient = "horizontal", background = "white", 
                      command = function(...) {
                        tkxview(canvas, ...)
                        setSrcLabel(viewLabel)
                      })
                    yscr <- tkscrollbar(f, repeatinterval = 5, 
                      orient = "vertical", command = function(...) {
                        tkyview(canvas, ...)
                        setSrcLabel(viewLabel)
                      })
                    f1 <- tkframe(f)
                    canvas <- tkcanvas(f1, relief = "raised", 
                      background = background, closeenough = close.enough, 
                      borderwidth = 5, highlightthickness = 2, 
                      scrollregion = c(min.x, min.y, max.x, max.y), 
                      xscrollincrement = -4, yscrollincrement = -4, 
                      xscrollcommand = function(...) tkset(xscr, ...), 
                      yscrollcommand = function(...) tkset(yscr, ...), 
                      width = width, height = height)
                    tkpack(canvas, expand = "yes", fill = "both")
                    tkgrid(f1, padx = 1, pady = 1, row = 0, column = 0, 
                      rowspan = 1, columnspan = 1, sticky = "news")
                    tkgrid(yscr, padx = 1, pady = 1, row = 0, 
                      column = 1, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(xscr, padx = 1, pady = 1, row = 1, 
                      column = 0, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid(viewLabel, padx = 1, pady = 1, row = 2, 
                      column = 0, rowspan = 1, columnspan = 1, 
                      sticky = "news")
                    tkgrid.columnconfigure(f, 0, weight = 1, 
                      minsize = 0)
                    tkgrid.rowconfigure(f, 0, weight = 1, minsize = 0)
                    tkpack(canvas, expand = "yes", fill = "both", 
                      padx = 1, pady = 1)
                    tkfocus(canvas)
                  }
                }
                else {
                  canvas <- tkcanvas(top, relief = "raised", 
                    background = background, closeenough = close.enough, 
                    width = width, height = height)
                  tkpack(canvas)
                }
                assign("box", f0.box, envir = top$env)
                assign("canvas", canvas, envir = top$env)
                assign("viewLabel", viewLabel, envir = top$env)
                assign("tags", tags, envir = top$env)
                GW.top <<- (top)
                GW.tags <<- list(NULL)
                GW.env <<- .Dg.toplevel(parent)
                prototype <- typeToPrototype(type = viewType, 
                  prototype = "DynamicGraphView", classes = control$viewClasses)
                result <- new(prototype, id.env = GW.env$ID, 
                  label = paste(title, viewType, sep = " / "), 
                  index = index, id = id, dg = ldg)
                tktitle(top) <- title
                return(result)
            }
            "Args" <- function(x.drawModel = drawModel, x.redrawView = redrawView, 
                x.frameModels = dgm.frameModels, x.frameViews = dm.frameViews, 
                x.graphWindow = GraphWindow, x.vertexList = vertexList, 
                x.blockList = blockList, x.dg = dg, x.visibleVertices = dg@visibleVertices, 
                x.visibleBlocks = dg@visibleBlocks, x.edgeList = currentEdges(edge.type = "VertexEdge"), 
                x.oriented = dg@oriented, x.blockEdgeList = currentEdges(edge.type = "BlockEdge"), 
                x.factorVertexList = dg@factorVertexList, x.factorEdgeList = currentEdges(edge.type = "FactorEdge"), 
                x.extraList = dg@extraList, x.extraEdgeList = currentEdges(edge.type = "ExtraEdge"), 
                x.object = object, x.objectName = objectName, 
                x.viewType = dg@viewType, x.top = GW.top, x.box = GW.top$env$box, 
                x.canvas = GW.top$env$canvas, x.viewLabel = GW.top$env$viewLabel, 
                x.tags = GW.tags, x.envir = GW.env$env, x.title = title, 
                x.selectedNodes = selectedNodes, x.selectedEdges = selectedEdges, 
                x.closedBlock = closedBlock, x.hiddenBlock = hiddenBlock, 
                x.control = control) return(list(drawModel = x.drawModel, 
                redrawView = x.redrawView, frameModels = x.frameModels, 
                frameViews = x.frameViews, graphWindow = x.graphWindow, 
                modelIndex = x.frameViews@index, viewIndex = x.graphWindow@index, 
                vertexList = x.vertexList, blockList = x.blockList, 
                dg = x.dg, visibleVertices = x.visibleVertices, 
                visibleBlocks = x.visibleBlocks, edgeList = x.edgeList, 
                oriented = x.oriented, blockEdgeList = x.blockEdgeList, 
                factorVertexList = x.factorVertexList, factorEdgeList = x.factorEdgeList, 
                extraList = x.extraList, extraEdgeList = x.extraEdgeList, 
                object = x.object, objectName = x.objectName, 
                viewType = x.viewType, top = x.top, box = x.box, 
                canvas = x.canvas, viewLabel = x.viewLabel, tags = x.tags, 
                envir = x.envir, title = x.title, selectedNodes = x.selectedNodes, 
                selectedEdges = x.selectedEdges, closedBlock = x.closedBlock, 
                hiddenBlock = x.hiddenBlock, control = x.control))
            "makeSlave" = function(sameModel = TRUE, local.Views = NULL, 
                Object = NULL, label = "Default", variableFrame = control$variableFrame) {
                control$variableFrame <- variableFrame
                sinkView(NULL, edges = FALSE, blocks = TRUE)
                edge.List <- currentEdges(edge.type = "VertexEdge")
                blockEdgeList <- currentEdges(edge.type = "BlockEdge")
                factorEdgeList <- currentEdges(edge.type = "FactorEdge")
                extraEdgeList <- currentEdges(edge.type = "ExtraEdge")
                Arguments <- Args()
                ldg <- dg
                ldg@edgeList <- edge.List
                if (any(slotNames(Object) == ".title")) {
                  label <- Object@.title
                }
                if (sameModel) 
                  redrawView(frameModels = dgm.frameModels, frameViews = local.Views, 
                    graphWindow = NULL, dg = ldg, control = control, 
                    Arguments = Arguments, title = label)
                else drawModel(frameModels = dgm.frameModels, 
                  frameViews = NULL, graphWindow = NULL, dg = ldg, 
                  object = Object, control = control, Arguments = Arguments, title = label)
            }
            "relativePositionsCanvas" <- function(positions) {
                if (!is.null(zoomPositions)) {
                  diff <- 100/(zoomPositions[, 2] - zoomPositions[, 
                    1])
                  p <- (diag(diff) %*% (.asRow(positions)) - 
                    0)
                }
                else p <- .asRow(positions)
                x <- t(diag(c(control$width, control$height, 
                  rep(100, local.N - 2))/100) %*% (p + 0))
                x <- (x - 0) * Scale + 0
                x <- round(x)
                return(x)
            }
            "inversCanvasRelativePosition" <- function(positions) {
                p <- .asRow(positions)
                p <- (p - 0)/Scale + 0
                if (is.null(dim(p)) && (length(p) < local.N)) 
                  if (length(p) == 1) 
                    p <- rep(p, local.N)
                  else p <- rep(0, local.N)
                p <- t(diag(c(100/control$width, 100/control$height, 
                  rep(1, local.N - 2))) %*% p - 0)
                if (!is.null(zoomPositions)) {
                  diff <- (zoomPositions[, 2] - zoomPositions[, 
                    1])/100
                  q <- t(diag(diff) %*% t(p + 0))
                  return(q)
                }
                else return(p)
            }
            "positionsCanvas" <- function(positions) {
                if (!is.null(zoomPositions)) {
                  a <- zoomPositions[, 1]
                  b <- zoomPositions[, 2]
                  A <- matrix(rep(a, ifelse(is.null(dim(positions)), 
                    1, nrow(positions))), byrow = TRUE, ncol = local.N)
                  diff <- 100/(b - a)
                  p <- (diag(diff) %*% (.asRow(positions) - t(A)) - 
                    50)
                }
                else p <- .asRow(positions)
                x <- t(diag(c(control$width, control$height, 
                  rep(100, local.N - 2))/100) %*% (p + 50))
                if (!is.null(zoomCenter)) 
                  x <- t((t(x) - c(zoomCenter)) * Scale + zoomCenter)
                x <- round(x)
                return(x)
            }
            "inversCanvasPosition" <- function(positions) {
                p <- .asRow(positions)
                if (!is.null(zoomCenter)) 
                  p <- t((t(p) - zoomCenter)/Scale + zoomCenter)
                if (is.null(dim(p)) && (length(p) < local.N)) 
                  if (length(p) == 1) 
                    p <- rep(p, local.N)
                  else p <- rep(0, local.N)
                p <- t(diag(c(100/control$width, 100/control$height, 
                  rep(1, local.N - 2))) %*% p - 50)
                if (!is.null(zoomPositions)) {
                  a <- zoomPositions[, 1]
                  b <- zoomPositions[, 2]
                  A <- matrix(rep(a, ifelse(is.null(dim(positions)), 
                    1, nrow(positions))), byrow = TRUE, ncol = local.N)
                  diff <- (b - a)/100
                  q <- t(diag(diff) %*% t(p + 50) + t(A))
                  return(q)
                }
                else return(p)
            }
            "replaceXY" <- function(x, y, position) {
                position[1] <- as.numeric(x)
                position[2] <- as.numeric(y)
                if (control$permitZoom) {
                  r.x <- as.numeric(tkxview(GW.top$env$canvas))
                  r.y <- as.numeric(tkyview(GW.top$env$canvas))
                  if (d.x == 0) 
                    m.x <- 0
                  else m.x <- -(r.x[1] - d.x)/d.x * min.x
                  if (d.y == 0) 
                    m.y <- 0
                  else m.y <- -(r.y[1] - d.y)/d.y * min.y
                  m <- c(m.x, m.y, rep(0, local.N - 2))
                  position <- position + m
                }
                return(position)
            }
            "callPopup" <- function(i, PopupMenu) {
                force(i)
                function(x, y) {
                  xCanvas <- as.integer(x) + as.integer(tkwinfo("rootx", 
                    canvas))
                  yCanvas <- as.integer(y) + as.integer(tkwinfo("rooty", 
                    canvas))
                  .Tcl(paste("tk_popup", .Tcl.args(PopupMenu, 
                    xCanvas, yCanvas)))
                }
            }
            "callPopupInBox" <- function(vertex, i, vertex.type, 
                UserNodePopupItems) {
                force(vertex)
                force(i)
                force(vertex.type)
                force(UserNodePopupItems)
                function(...) {
                  nodePopupMenu <- tkmenu(GW.top$env$box, tearoff = FALSE)
                  addNodePopups(vertex, i, vertex.type, nodePopupMenu, 
                    UserNodePopupItems, slave = FALSE)
                  xCanvas <- as.integer(tkwinfo("rootx", GW.top$env$box))
                  yCanvas <- as.integer(tkwinfo("rooty", GW.top$env$box))
                  .Tcl(paste("tk_popup", .Tcl.args(nodePopupMenu, 
                    xCanvas, yCanvas)))
                }
            }
            "callPopupNode" <- function(vertex, i, vertex.type, 
                UserNodePopupItems) {
                force(vertex)
                force(i)
                force(vertex.type)
                force(UserNodePopupItems)
                function(x, y) {
                  nodePopupMenu <- tkmenu(canvas, tearoff = FALSE)
                  addNodePopups(vertex, i, vertex.type, nodePopupMenu, 
                    UserNodePopupItems)
                  xCanvas <- as.integer(x) + as.integer(tkwinfo("rootx", 
                    canvas))
                  yCanvas <- as.integer(y) + as.integer(tkwinfo("rooty", 
                    canvas))
                  .Tcl(paste("tk_popup", .Tcl.args(nodePopupMenu, 
                    xCanvas, yCanvas)))
                }
            }
            "callPopupEdge" <- function(edge, i, f, t, edge.type, 
                U.Menus) {
                force(edge)
                force(i)
                force(f)
                force(t)
                force(edge.type)
                force(U.Menus)
                function(x, y) {
                  edgePopupMenu <- tkmenu(canvas, tearoff = FALSE)
                  addEdgePopups(canvas, edge, i, f, t, edgePopupMenu, 
                    U.Menus, edge.type)
                  xCanvas <- as.integer(x) + as.integer(tkwinfo("rootx", 
                    canvas))
                  yCanvas <- as.integer(y) + as.integer(tkwinfo("rooty", 
                    canvas))
                  .Tcl(paste("tk_popup", .Tcl.args(edgePopupMenu, 
                    xCanvas, yCanvas)))
                }
            }
            "getTag" <- function(text, number, setTag = TRUE) {
                tag <- paste(text, abs(number), GraphWindow@id, 
                  sep = ".")
                if (is.null(GW.tags[[1]])) 
                  GW.tags <<- list(tag)
                else if (!any(unlist(lapply(GW.tags, function(i) i == 
                  tag)))) 
                  GW.tags <<- append(list(tag), GW.tags)
                else if (setTag) 
                  message(paste("(( Duplicated tag: ", tag, " ))", 
                    sep = " "))
                if (control$saveTkReferences) 
                  assign("tags", GW.tags, envir = GW.env$env)
                return(tag)
            }
            "deleteTags" <- function(text = "deleteTags: ") {
                for (i in seq(along = GW.tags)) {
                  tkdelete(GW.top$env$canvas, GW.tags[[i]])
                }
                for (i in seq(along = GW.tags)) {
                  tag <- tkgettags(GW.top$env$canvas, i)
                  if (length(as.character(tag)) > 0) {
                    if (FALSE || control$debug.edges) {
                      print(paste(text, i, as.character(tag), 
                        sep = ": "))
                    }
                    tkdelete(GW.top$env$canvas, i)
                  }
                }
            }
            "destroyView" <- function(deleteTags = FALSE, txt = "") {
                function(...) {
                  if (deleteTags) 
                    deleteTags("destroyView")
                  updateWindow <<- FALSE
                  tkdestroy(GW.top)
                }
            }
            "sinkVertexList" <- function() {
                vertexList <<- dgm.frameModels@vertices
                for (i in seq(along = vertexList)) {
                  position <- positionsVertices[i, ]
                  position(vertexList[[i]]) <<- position
                  position <- positionsLabels[i, ]
                  labelPosition(vertexList[[i]]) <<- position
                  label(vertexList[[i]]) <<- Labels[i]
                  vertexList[[i]]@name <<- namesVertices[i]
                  color(vertexList[[i]]) <<- colorsVertices[i]
                  blockindex(vertexList[[i]]) <<- blocksVertices[i]
                  if (!.IsEmpty(blockList)) 
                    if (blocksVertices[i] == 0) 
                      stratum(vertexList[[i]]) <<- 0
                    else stratum(vertexList[[i]]) <<- strataBlocks[blocksVertices[i]]
                }
                dgm.frameModels@vertices <<- vertexList
                return(vertexList)
            }
            "sinkEdgeList" <- function() {
                dg@edgeList <<- GraphWindow@dg@edgeList
                for (i in seq(along = dg@edgeList)) {
                  position <- positionsEdgeLabels[i, ]
                  dg@edgeList[[i]]@label.position <<- position
                }
                GraphWindow@dg@edgeList <<- dg@edgeList
                return(dg@edgeList)
            }
            "sinkBlockTree" <- function(tree) {
                "subSinkBlockTree" <- function(tree) {
                  i <- abs(tree$block@index)
                  position <- positionsBlocks[i, , ]
                  position(tree$block) <<- position
                  position <- positionsBlockLabels[i, ]
                  labelPosition(tree$block) <<- position
                  tree$block@stratum <<- strataBlocks[i]
                  tree$block@label <<- blockLabels[i]
                  tree$block@closed <<- closedBlock[i]
                  tree$block@visible <<- !hiddenBlock[i]
                  if (!is.null((tree$sub.blocks))) 
                    for (j in 1:length(tree$sub.blocks)) {
                      tree$sub.blocks[[j]] <<- subSinkBlockTree(tree$sub.blocks[[j]])
                    }
                  return(tree)
                }
                if (!.IsEmpty(blockList) && !.IsEmpty(tree) && 
                  !(length(tree) == 0)) 
                  subSinkBlockTree(tree)
                return(tree)
            }
            "sinkBlockList" <- function(sinkTree = TRUE) {
                if (!.IsEmpty(blockList)) 
                  for (i in seq(along = blockList)) {
                    position <- positionsBlocks[i, , ]
                    position(blockList[[i]]) <<- position
                    position <- positionsBlockLabels[i, ]
                    labelPosition(blockList[[i]]) <<- position
                    blockList[[i]]@stratum <<- strataBlocks[i]
                    blockList[[i]]@label <<- blockLabels[i]
                    blockList[[i]]@closed <<- closedBlock[i]
                    blockList[[i]]@visible <<- !hiddenBlock[i]
                  }
                if (sinkTree && !is.null(blockTree)) 
                  sinkBlockTree(blockTree)
                if (is.null(blockList)) {
                  dgm.frameModels@blocks <<- new("dg.BlockList")
                  return(dgm.frameModels@blocks)
                }
                else {
                  dgm.frameModels@blocks <<- blockList
                  return(blockList)
                }
            }
            "sinkBlockEdges" <- function() {
                dg@blockEdgeList <<- GraphWindow@dg@blockEdgeList
                return(dg@blockEdgeList)
            }
            "sinkFactorVertexList" <- function() {
                dg@factorVertexList <<- GraphWindow@dg@factorVertexList
                for (i in seq(along = dg@factorVertexList)) {
                  position <- positionsFactorVertices[i, ]
                  position(dg@factorVertexList[[i]]) <<- position
                  position <- positionsFactorLabels[i, ]
                  labelPosition(dg@factorVertexList[[i]]) <<- position
                  label(dg@factorVertexList[[i]]) <<- factorLabels[i]
                  dg@factorVertexList[[i]]@name <<- namesFactorVertices[i]
                  color(dg@factorVertexList[[i]]) <<- colorsFactorVertices[i]
                  if (!.IsEmpty(blockList)) {
                    blockindex(dg@factorVertexList[[i]]) <<- blocksFactorVertices[i]
                    stratum(dg@factorVertexList[[i]]) <<- strataBlocks[blocksFactorVertices[i]]
                  }
                }
                GraphWindow@dg@factorVertexList <<- dg@factorVertexList
                return(dg@factorVertexList)
            }
            "sinkFactorEdgeList" <- function() {
                dg@factorEdgeList <<- GraphWindow@dg@factorEdgeList
                return(dg@factorEdgeList)
            }
            "sinkExtraVertexList" <- function() {
                dg@extraList <<- GraphWindow@dg@extraList
                for (i in seq(along = dg@extraList)) {
                  position <- positionsExtraVertices[i, ]
                  position(dg@extraList[[i]]) <<- position
                  position <- positionsExtraLabels[i, ]
                  labelPosition(dg@extraList[[i]]) <<- position
                  label(dg@extraList[[i]]) <<- extraLabels[i]
                  dg@extraList[[i]]@name <<- namesExtraVertices[i]
                  color(dg@extraList[[i]]) <<- colorsExtraVertices[i]
                  if (!.IsEmpty(blockList)) {
                    blockindex(dg@extraList[[i]]) <<- blocksExtraVertices[i]
                    stratum(dg@extraList[[i]]) <<- strataBlocks[blocksExtraVertices[i]]
                  }
                }
                GraphWindow@dg@extraList <<- dg@extraList
                return(dg@extraList)
            }
            "sinkExtraEdgeList" <- function() {
                dg@extraEdgeList <<- GraphWindow@dg@extraEdgeList
                return(dg@extraEdgeList)
            }
            "sinkView" <- function(menuItem, vertices = TRUE, 
                edges = TRUE, blocks = FALSE) {
                if (vertices && (is.null(menuItem$update.vertices) || 
                  menuItem$update.vertices)) {
                  V <- sinkVertexList()
                }
                if (edges && (is.null(menuItem$update.edges) || 
                  menuItem$update.edges)) {
                  E <- sinkEdgeList()
                  V <- sinkFactorVertexList()
                  E <- sinkFactorEdgeList()
                  V <- sinkExtraVertexList()
                  E <- sinkExtraEdgeList()
                }
                if (blocks && (is.null(menuItem$update.blocks) || 
                  menuItem$update.blocks)) 
                  B <- sinkBlockList()
                E <- sinkBlockEdges()
                if (!is.null(blockTree)) 
                  BT <- sinkBlockTree(blockTree)
                m <- dm.frameViews@index
                n <- GraphWindow@index
                dm.frameViews@graphs[[n]] <<- GraphWindow
                dgm.frameModels@models[[m]]@graphs[[n]] <<- GraphWindow
            }
            "update" <- function(type = "Arguments", ...) {
                if (type == "Arguments") 
                  sinkView(...)
                else subUpdatePositions(...)
            }
            "subSinkAllFrames" <- function(type = "Arguments", 
                ...) {
                V <- sinkVertexList()
                for (m in 1:length(dgm.frameModels@models)) {
                  title <- paste("Model:", m, sep = " ")
                  frame.view <- dgm.frameModels@models[[m]]
                  dgm.frameModels@models[[m]]@model <<- list(object)
                  for (n in 1:length(frame.view@graphs)) {
                    graph.window <- frame.view@graphs[[n]]
                    if (is.null(graph.window)) {
                      message(paste("Empty window: ", title, 
                        "; Graph:", n, sep = " "))
                    }
                    else {
                      gw.env <- .get.env.graphWindow(graphWindow = graph.window, 
                        frameViews = frame.view, frameModels = dgm.frameModels)$env
                      if (is.null(formals(gw.env$Update))) {
                        if (control$debug.update) {
                          txt <- paste("(( No update function: ", 
                            title, "; Graph:", n, graph.window@label, 
                            " ))", sep = " ")
                          message(txt)
                        }
                      }
                      else gw.env$Update(type, ...)
                    }
                  }
                }
            }
            "sinkAllFrames" <- function(type = "Arguments", txt) {
                force(type)
                force(txt)
                function(...) {
                  subSinkAllFrames(type, txt)
                }
            }
            "objectAssign" <- function(R) {
                if (!is.null(R) && !is.null(R$object)) {
                  if (!is.null(objectName)) 
                    assign(objectName, R$object, pos = 1)
                }
            }
            "setModel" <- function(R.object = NULL, dg = NULL, 
                txt = "", graphWindow = NULL, copyProperties = FALSE, 
                setUpdate = TRUE, RR = NULL) {
                if (is.null(dg)) 
                  dg <- .newDgGraphEdges(viewType = GraphWindow@dg@viewType, 
                    oriented = GraphWindow@dg@oriented, vertexList = vertexList, 
                    visibleVertices = GraphWindow@dg@visibleVertices, 
                    visibleBlocks = GraphWindow@dg@visibleBlocks, 
                    edgeList = copyCurrentEdges(edge.type = "VertexEdge", 
                      copyProperties = copyProperties), blockList = blockList, 
                    blockEdgeList = copyCurrentEdges(edge.type = "BlockEdge", 
                      copyProperties = copyProperties), factorVertexList = GraphWindow@dg@factorVertexList, 
                    factorEdgeList = copyCurrentEdges(edge.type = "FactorEdge", 
                      copyProperties = copyProperties), extraList = GraphWindow@dg@extraList, 
                    extraEdgeList = copyCurrentEdges(edge.type = "ExtraEdge", 
                      copyProperties = copyProperties))
                if (control$debug.update) {
                  print(paste("setModel:", txt))
                  print(c(updateCountModel, updateCountModelMain))
                }
                if (hasMethod("setGraphEdges", class(R.object))) {
                  # message("Using 'setGraphEdges' for your model class.")
                  object <<- setGraphEdges(R.object, dg = dg, 
                    ...)
                }
                else if (hasMethod("setGraphComponents", class(R.object))) {
                  message("Please implement 'setGraphEdges' for your model class.")
                  if (is.null(dg@edgeList)) 
                    dg@edgeList <- new("dg.VertexEdgeList")
                  object <<- setGraphComponents(R.object, viewType = dg@viewType, 
                    visibleVertices = dg@visibleVertices, visibleBlocks = dg@visibleBlocks, 
                    extraVertices = dg@extraList, vertexEdges = dg@edgeList, 
                    blockEdges = dg@blockEdgeList, factorVertices = dg@factorVertexList, 
                    factorEdges = dg@factorEdgeList, extraEdges = dg@extraEdgeList)
                }
                if (setUpdate) {
                  updateCountModelMain <<- updateCountModelMain + 
                    1
                  updateCountModel <<- updateCountModelMain
                }
            }
            "updateModel" <- function() {
                if (control$debug.update) 
                  print(paste("updateModel", getLabel()))
                if (hasMethod("returnGraphComponents", class(object)) || 
                  hasMethod("graphComponents", class(object)) || 
                  hasMethod("graphEdges", class(object))) {
                  tkconfigure(canvas, cursor = "watch")
                  tkconfigure(GW.top$env$viewLabel, text = paste(dg@viewType, 
                    " | Working !!!"))
                  tkfocus(GW.top)
                  Arguments <- Args()
                  if (hasMethod("graphEdges", class(object))) {
                    graphContent <- graphEdges(object, viewType = dg@viewType, 
                      Arguments = Arguments)
                  }
                  else {
                    message("Please implement 'graphEdges' for your model class.")
                    if (hasMethod("graphComponents", class(object))) 
                      graphContent <- graphComponents(object, 
                        viewType = dg@viewType, Arguments = Arguments)
                    else  {
                    warning("'returnGraphComponents' remove to aviod 'NOTE' by 'R CMD check'!")
                      # graphContent <- returnGraphComponents(object, 
                      #                    viewType = dg@viewType, Arguments = Arguments)
                    }
                  }
                  if ((class(graphContent) == "dg.graph") || 
                    (class(graphContent) == "dg.graphedges")) {
                    redrawView(frameModels = dgm.frameModels, 
                      frameViews = dm.frameViews, graphWindow = GraphWindow, 
                      dg = graphContent, control = control, Arguments = Arguments)
                  }
                  else {
                    message("Please return object of class 'dg.graphedges'.")
                    if (is.list(graphContent)) {
                      names <- names(graphContent)
                      checkClass <- function(name, class, z = "$", 
                        a = "graphContent", b = paste(a, z, name, 
                          sep = "")) {
                        text <- paste(c("if ((\"", name, "\" %in% names) ", 
                          "&& (class(", b, ") != \"", class, 
                          "\")) ", "{ message(paste(\"Invalid class of '", 
                          name, "' in list from 'graphComponents'; \")); ", 
                          b, " <<- new(\"", class, "\", .nullToList(", 
                          b, ")) }"), collapse = "")
                        eval(parse(text = text))
                      }
                      checkClass("extraVertices", "dg.VertexList")
                      checkClass("vertexEdges", "dg.VertexEdgeList")
                      checkClass("graphEdges", "dg.VertexEdgeList")
                      checkClass("blockEdges", "dg.BlockEdgeList")
                      checkClass("factorVertices", "dg.FactorVertexList")
                      checkClass("factorEdges", "dg.FactorEdgeList")
                      checkClass("extraEdges", "dg.ExtraEdgeList")
                    }
                    ldg <- .newDgGraphEdges(vertexList = vertexList, 
                      visibleVertices = graphContent$visibleVertices, 
                      visibleBlocks = graphContent$visibleBlocks, 
                      edgeList = graphContent$vertexEdges, blockList = blockList, 
                      blockEdgeList = graphContent$blockEdges, 
                      factorVertexList = graphContent$factorVertices, 
                      factorEdgeList = graphContent$factorEdges, 
                      extraList = graphContent$extraVertices, 
                      extraEdgeList = graphContent$extraEdges)
                    redrawView(frameModels = dgm.frameModels, 
                      frameViews = dm.frameViews, graphWindow = GraphWindow, 
                      dg = ldg, control = control, Arguments = Arguments)
                  }
                  tkconfigure(GW.top$env$viewLabel, text = dg@viewType)
                  tkconfigure(canvas, cursor = "arrow")
                }
            }
            "testUpdateModel" <- function() {
                if ((updateCountModel < updateCountModelMain)) {
                  updateModel()
                  updateCountModel <<- updateCountModelMain
                }
            }
            "extractEdgesResult" <- function(R, newEdges, from.R.edgeList = TRUE, 
                title) {
                if (!from.R.edgeList || is.null(R$edgeList)) 
                  if (is.null(R$newEdges$vertexEdges)) 
                    Edges <- newEdges$vertexEdges
                  else if (is.null(R$dg)) 
                    Edges <- R$newEdges$vertexEdges
                  else Edges <- R$dg@edgeList
                else Edges <- R$edgeList
                return(Edges)
            }
            "getEdges" <- function(edge.type = "VertexEdge") {
                if (edge.type == "VertexEdge") 
                  return(GraphWindow@dg@edgeList)
                else if (edge.type == "FactorEdge") 
                  return(GraphWindow@dg@factorEdgeList)
                else if (edge.type == "ExtraEdge") 
                  return(GraphWindow@dg@extraEdgeList)
                else if (edge.type == "BlockEdge") 
                  return(GraphWindow@dg@blockEdgeList)
                else return(NULL)
            }
            "edgesClass" <- function(edge.type = "VertexEdge") if (edge.type == 
                "VertexEdge") 
                return("dg.VertexEdgeList")
            else if (edge.type == "FactorEdge") 
                return("dg.FactorEdgeList")
            else if (edge.type == "ExtraEdge") 
                return("dg.ExtraEdgeList")
            else if (edge.type == "BlockEdge") 
                return("dg.BlockEdgeList")
            else return(NULL)
            "currentEdges" <- function(edge.type = "VertexEdge") {
                E <- getEdges(edge.type = edge.type)
                if (length(E) > 0) {
                  E <- lapply(E, function(egde) if (sum(abs(egde@vertex.indices)) > 
                    0) 
                    egde)
                  E <- .removeNull(E)
                }
                else E <- NULL
                if (!is.null(E)) {
                  class(E) <- edgesClass(edge.type = edge.type)
                }
                return(E)
            }
            "append.index.edge" <- function(e, edge.type = "VertexEdge", 
                edgeClass = NULL) {
                if (edge.type == "VertexEdge") 
                  new.edge <- returnEdgeList(list(e), vertexList, 
                    color = control$edgeColor, oriented = dg@oriented, 
                    types = edgeClass, N = local.N, edgeClasses = control$edgeClasses)
                else if (edge.type == "FactorEdge") 
                  new.edge <- returnFactorEdgeList(list(e), vertexList, 
                    color = control$factorEdgeColor, dg@factorVertexList)
                else if (edge.type == "ExtraEdge") 
                  new.edge <- returnExtraEdgeList(list(e), vertexList, 
                    color = control$extraEdgeColor, dg@extraList)
                else if (edge.type == "BlockEdge") 
                  new.edge <- new("dg.BlockEdgeList")
                if (!is.na(control$namesOnEdges) && !control$namesOnEdges) {
                  label(new.edge[[1]]) <- ""
                }
                E <- append(getEdges(edge.type = edge.type), 
                  new.edge)
                class(E) <- edgesClass(edge.type = edge.type)
                if (edge.type == "VertexEdge") 
                  GraphWindow@dg@edgeList <<- E
                else if (edge.type == "FactorEdge") 
                  GraphWindow@dg@factorEdgeList <<- E
                else if (edge.type == "ExtraEdge") 
                  GraphWindow@dg@extraEdgeList <<- E
                else if (edge.type == "BlockEdge") 
                  GraphWindow@dg@blockEdgeList <<- E
                return(E)
            }
            "append.edge" <- function(e, edge.type = "VertexEdge") {
                E <- append(getEdges(edge.type = edge.type), 
                  list(e))
                class(E) <- edgesClass(edge.type = edge.type)
                if (edge.type == "VertexEdge") 
                  GraphWindow@dg@edgeList <<- E
                else if (edge.type == "FactorEdge") 
                  GraphWindow@dg@factorEdgeList <<- E
                else if (edge.type == "ExtraEdge") 
                  GraphWindow@dg@extraEdgeList <<- E
                else if (edge.type == "BlockEdge") 
                  GraphWindow@dg@blockEdgeList <<- E
                return(E)
            }
            "selectCurrentEdges" <- function(omitEdges = FALSE, 
                edge.type = "VertexEdge") {
                E <- getEdges(edge.type = edge.type)
                if (length(E) > 0) {
                  j <- omitEdges | vertex.in.edge(0, edge.type = edge.type)
                  if (edge.type == "VertexEdge") 
                    j <- j | non.graph.edge(edge.type = edge.type)
                  E <- sapply(1:length(E), function(x) if (!j[x]) 
                    E[[x]])
                  E <- .removeNull(E)
                }
                else E <- NULL
                return(E)
            }
            "copyCurrentEdges" <- function(omitEdges = FALSE, 
                edge.type = "VertexEdge", copyProperties = FALSE) {
                E <- getEdges(edge.type = edge.type)
                if (length(E) > 0) {
                  j <- omitEdges | vertex.in.edge(0, edge.type = edge.type)
                  if (edge.type == "VertexEdge") 
                    j <- j | non.graph.edge(edge.type = edge.type)
                  edge.classes <- lapply(1:length(E), function(x) if (!j[x]) 
                    class(E[[x]]))
                  edge.classes <- .removeNull(edge.classes)
                  if (copyProperties) {
                    edge.widths <- lapply(1:length(E), function(x) if (!j[x]) 
                      width(E[[x]]))
                    edge.widths <- .removeNull(edge.widths)
                    edge.colors <- lapply(1:length(E), function(x) if (!j[x]) 
                      color(E[[x]]))
                    edge.colors <- .removeNull(edge.colors)
                    edge.dashs <- lapply(1:length(E), function(x) if (!j[x]) 
                      dash(E[[x]]))
                    edge.dashs <- .removeNull(edge.dashs)
                  }
                  if (edge.type == "BlockEdge") {
                    blockEdges <- lapply(1:length(E), function(x) if (!j[x]) 
                      E[[x]])
                    blockEdges <- .removeNull(blockEdges)
                    class(blockEdges) <- "dg.BlockEdgeList"
                  }
                  else {
                    E <- lapply(1:length(E), function(x) if (!j[x]) 
                      E[[x]]@vertex.indices)
                  }
                  E <- .removeNull(E)
                }
                else {
                  edge.classes <- NULL
                  if (copyProperties) {
                    edge.widths <- NULL
                    edge.colors <- NULL
                    edge.dashs <- NULL
                  }
                  E <- NULL
                  blockEdges <- NULL
                }
                if (edge.type == "VertexEdge") 
                  E <- returnEdgeList(E, vertexList, types = edge.classes, 
                    color = control$edgeColor, oriented = dg@oriented, 
                    N = local.N, edgeClasses = control$edgeClasses)
                else if (edge.type == "FactorEdge") 
                  E <- returnFactorEdgeList(E, vertexList, dg@factorVertexList, 
                    color = control$factorEdgeColor)
                else if (edge.type == "ExtraEdge") 
                  E <- returnExtraEdgeList(E, vertexList, dg@extraList, 
                    color = control$extraEdgeColor)
                else if (edge.type == "BlockEdge") 
                  E <- blockEdges
                if (copyProperties && !is.null(E)) {
                  Widths(E) <- unlist(edge.widths)
                  Colors(E) <- unlist(edge.colors)
                  Dashes(E) <- unlist(edge.dashs)
                }
                return(E)
            }
            "appendToCurrentEdges" <- function(omitEdges = FALSE, 
                new.edge = NULL, edge.type = "VertexEdge", edgeClass = NULL) {
                E <- getEdges(edge.type = edge.type)
                if (length(E) > 0) {
                  j <- omitEdges | vertex.in.edge(0, edge.type = edge.type)
                  if (edge.type == "VertexEdge") 
                    j <- j | non.graph.edge(edge.type = edge.type)
                  E <- lapply(1:length(E), function(x) if (!j[x]) 
                    E[[x]])
                  E <- .removeNull(E)
                  edge.list <- lapply(E, function(e) e@vertex.indices)
                  edge.classes <- lapply(E, function(e) class(e))
                  if (!is.null(new.edge)) {
                    edge.list <- append(edge.list, new.edge)
                    edge.classes <- append(edge.classes, edgeClass)
                  }
                }
                else {
                  edge.classes <- edgeClass
                  edge.list <- new.edge
                }
                if (edge.type == "VertexEdge") 
                  E <- returnEdgeList(edge.list, vertexList, 
                    types = edge.classes, color = control$edgeColor, 
                    oriented = dg@oriented, N = local.N, edgeClasses = control$edgeClasses)
                else if (edge.type == "FactorEdge") 
                  E <- NULL
                else if (edge.type == "ExtraEdge") 
                  E <- NULL
                else if (edge.type == "BlockEdge") 
                  E <- NULL
                return(E)
            }
            "which.unordered.edge" <- function(e, edge.type = "VertexEdge") {
                n <- length(e)
                unlist(lapply(getEdges(edge.type = edge.type), 
                  function(i) length(e[!is.na(match(e, i@vertex.indices))]) == 
                    n))
            }
            "which.edge" <- function(e, edge.type = "VertexEdge") unlist(lapply(getEdges(edge.type = edge.type), 
                function(i) all(i@vertex.indices == e)))
            "vertex.in.edge" <- function(e, edge.type = "VertexEdge") unlist(lapply(getEdges(edge.type = edge.type), 
                function(i) is.element(e, i@vertex.indices)))
            "non.graph.edge" <- function(edge.type = "VertexEdge") unlist(lapply(getEdges(edge.type = edge.type), 
                function(i) any(i@vertex.indices <= 0)))
            "edge.vertices" <- function(i, type.negative = "Factor", 
                edge.type = "VertexEdge") {
                E <- getEdges(edge.type = edge.type)
                edge <- E[[i]]@vertex.indices
                edge.vertices <- vector("list", length(edge))
                for (j in seq(along = edge)) if (edge[j] > 0) 
                  edge.vertices[[j]] <- vertexList[[edge[j]]]
                else if (type.negative == "Factor") 
                  edge.vertices[[j]] <- dg@factorVertexList[[-edge[j]]]
                else if (type.negative == "Extra") 
                  edge.vertices[[j]] <- dg@extraList[[-edge[j]]]
                else if (type.negative == "ClosedBlock") 
                  edge.vertices[[j]] <- blockList[[-edge[j]]]
                return(edge.vertices)
            }
            "edge.negative.type" <- function(edge.type = "VertexEdge") if (edge.type == 
                "VertexEdge") 
                return("Vertex")
            else if (edge.type == "FactorEdge") 
                return("Factor")
            else if (edge.type == "ExtraEdge") 
                return("Extra")
            else if (edge.type == "BlockEdge") 
                return("ClosedBlock")
            "edge.names" <- function(i, type.negative = edge.negative.type(edge.type), 
                edge.type = "VertexEdge") lapply(edge.vertices(i, 
                type.negative = type.negative, edge.type = edge.type), 
                function(v) retVertexName(v@index, vertex.type = ifelse(v@index > 
                  0, "Vertex", type.negative)))
            "edge.positions" <- function(i, type.negative = edge.negative.type(edge.type), 
                edge.type = "VertexEdge") {
                lapply(edge.vertices(i, type.negative = type.negative, 
                  edge.type = edge.type), function(v) retVertexPos(v@index, 
                  ifelse(v@index > 0, "Vertex", type.negative)))
            }
            "edge.strata" <- function(i, type.negative = edge.negative.type(edge.type), 
                edge.type = "VertexEdge") {
                lapply(edge.vertices(i, type.negative = type.negative, 
                  edge.type = edge.type), function(v) retStratum(v@index, 
                  vertex.type = ifelse(v@index > 0, "Vertex", 
                    type.negative)))
            }
            "clearEdge" <- function(i, edge.type = "VertexEdge") if (edge.type == 
                "VertexEdge") 
                GraphWindow@dg@edgeList[[i]]@vertex.indices <<- c(0, 
                  0)
            else if (edge.type == "FactorEdge") 
                GraphWindow@dg@factorEdgeList[[i]]@vertex.indices <<- c(0, 
                  0)
            else if (edge.type == "ExtraEdge") 
                GraphWindow@dg@extraEdgeList[[i]]@vertex.indices <<- c(0, 
                  0)
            else if (edge.type == "BlockEdge") 
                GraphWindow@dg@blockEdgeList[[i]]@vertex.indices <<- c(0, 
                  0)
            "from" <- function(i, edge.type = "VertexEdge") getEdges(edge.type = edge.type)[[i]]@vertex.indices[1]
            "to" <- function(i, edge.type = "VertexEdge") getEdges(edge.type = edge.type)[[i]]@vertex.indices[2]
            transformation <- control$transformation
            "setTransformation" <- function(value = NULL) {
                if (is.null(value) == (!is.null(transformation))) {
                  if (!.IsEmpty(blockList)) 
                    for (i in seq(along = blockList)) if ((closedBlock[i] || 
                      hiddenBlock[i])) {
                    }
                    else deleteBlock(i)
                  transformation <<- value
                  if (!.IsEmpty(blockList)) 
                    for (i in seq(along = blockList)) if ((closedBlock[i] || 
                      hiddenBlock[i])) {
                    }
                    else if (is.element(i, dg@visibleBlocks)) 
                      drawBlock(blockList[[i]], i, setTag = FALSE)
                }
                else transformation <<- value
                subUpdateGraphWindow("setTransformation", all.blockframes = TRUE)
            }
            "angle" <- function(value = NULL) if (!is.null(value)) 
                Angle <<- value
            else return(Angle)
            "project" <- function(position) if (!is.null(transformation)) 
                t(transformation %*% .asRow(position))
            else position
            "inversProject" <- function(position) if (!is.null(transformation)) 
                t(solve(transformation, .asRow(position)))
            else position
            "applyTransformation" <- function(trans, draw.box = FALSE, 
                redraw = TRUE) {
                if (!is.null(transformation)) {
                  transformation <<- transformation %*% trans
                  if (redraw) 
                    subUpdateGraphWindow("applyTransformation", 
                      all.blockframes = TRUE)
                }
            }
            "sphereRand" <- function(n) {
                nx2 <- 2
                while ((nx2 >= 1)) {
                  x <- 2 * runif(n) - 1
                  nx2 <- sum(x^2)
                }
                return(x/sqrt(nx2))
            }
            "makeRotation" <- function(x, y, alpha = 0, use.alpha = FALSE, 
                n = length(x)) {
                if (length(x) != length(y) || !is.null(dim(x)) || 
                  !is.null(dim(y))) 
                  stop("Invalid arguments")
                "dnrm2" <- function(x) sqrt(sum(x^2))
                rot <- diag(1, n)
                nx <- dnrm2(x)
                ny <- dnrm2(y)
                if ((nx == 0) || (ny == 0)) 
                  return(rot)
                x <- x * (1/nx)
                y <- y * (1/ny)
                xy <- t(x) %*% y
                c <- ifelse(use.alpha, cos(alpha), xy)
                cc <- 1 - c^2
                s <- ifelse(use.alpha, sin(alpha), ifelse(cc > 
                  0, sqrt(cc), 0))
                cm1 <- c - 1
                y <- y - xy * x
                ny <- dnrm2(y)
                if (ny == 0) 
                  return(rot)
                y <- y * (1/ny)
                a <- x * cm1 + y * s
                b <- -x * s + y * cm1
                rot <- rot + a %*% t(x) + b %*% t(y)
                return(rot)
            }
            "canvasToSphere" <- function(X) {
                rad <- 100
                pos <- inversCanvasRelativePosition(X)
                x <- pos[1]
                y <- pos[2]
                norm.2 <- x^2 + y^2
                rad.2 <- rad^2
                z <- sqrt(max(rad.2 - norm.2, 0))
                if (local.N > 2) 
                  res <- c(x, y, z, rep(0, local.N - 3))
                else res <- c(x, y)
                if (norm.2 < rad.2) 
                  return(res)
                else {
                  r <- sqrt(norm.2/rad.2)
                  return(res/r)
                }
            }
            "doHandRotate" <- function() {
                p <- NULL
                function(x, y) {
                  tkconfigure(canvas, cursor = "watch")
                  tkfocus(canvas)
                  X <- replaceXY(x, y, rep(50, local.N))
                  if (is.null(p)) 
                    p <<- canvasToSphere(X)
                  else {
                    oldp <- p
                    p <<- canvasToSphere(X)
                    applyTransformation(makeRotation(oldp, p), 
                      draw.box = FALSE, redraw = TRUE)
                    tkconfigure(canvas, cursor = "arrow")
                  }
                }
            }
            "rockPlot" <- function(k = 2) {
                function(x, y) {
                  tkconfigure(canvas, cursor = "watch")
                  tkfocus(GW.top)
                  print("rockPlot")
                  angle <- 10
                  p1 <- sphereRand(control$N) # 'N' ???
                  p2 <- sphereRand(control$N) # 'N' ???
                  for (i in 1:k) applyTransformation(makeRotation(p1, 
                    p2, alpha = angle, use.alpha = TRUE), draw.box = FALSE, 
                    redraw = TRUE)
                  for (i in 1:(2 * k)) applyTransformation(makeRotation(p1, 
                    p2, alpha = -angle, use.alpha = TRUE), draw.box = FALSE, 
                    redraw = TRUE)
                  for (i in 1:k) applyTransformation(makeRotation(p1, 
                    p2, alpha = angle, use.alpha = TRUE), draw.box = FALSE, 
                    redraw = TRUE)
                  print("Finished rocking!")
                  tkconfigure(canvas, cursor = "arrow")
                }
            }
            "keyRotate" <- function(v = 0, sign = 1) {
                force(v)
                force(sign)
                function(...) {
                  if ((v > 2) || (local.N > 2)) {
                    v1 <- ifelse(v == 1, 2, 1)
                    v2 <- ifelse(v == 3, 2, 3)
                    if (is.null(transformation)) 
                      X <- diag(local.N)
                    else X <- transformation
                    angle <- pi/16
                    applyTransformation(makeRotation((X[, v1]), 
                      (X[, v2]), alpha = ifelse(sign == 1, angle, 
                        -angle), use.alpha = TRUE), draw.box = FALSE, 
                      redraw = TRUE)
                  }
                }
            }
            "vertexItem" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (vertex.type == "ClosedBlock") 
                  return(itemsClosedBlocks[[i]])
                else if (vertex.type == "Vertex") 
                  return(itemsVertices[[i]])
                else if (vertex.type == "Factor") 
                  return(itemsFactors[[-i]])
                else if (vertex.type == "Extra") 
                  return(itemsExtras[[abs(i)]])
            }
            "setVertexItem" <- function(i, value, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (vertex.type == "ClosedBlock") 
                  itemsClosedBlocks[[i]] <<- value
                else if (vertex.type == "Vertex") 
                  itemsVertices[[i]] <<- value
                else if (vertex.type == "Factor") 
                  itemsFactors[[abs(i)]] <<- value
                else if (vertex.type == "Extra") 
                  itemsExtras[[i]] <<- value
            }
            "edgeItem" <- function(i, edge.type = "VertexEdge") {
                if (i > 0) 
                  return(itemsEdges[[i]])
                else if (edge.type == "BlockEdge") {
                  if (is.element(abs(i), dg@visibleBlocks)) 
                    return(itemsBlockEdges[[-i]])
                  else return(NULL)
                }
                else if (edge.type == "FactorEdge") 
                  return(itemsFactorEdges[[-i]])
                else return(itemsExtraEdges[[-i]])
            }
            "subSetEdgeItem" <- function(from, edge.type, i, 
                edgeNode) {
                if (from > 0) {
                  itemsEdges[[from]][[i]] <<- edgeNode
                }
                else if (edge.type == "BlockEdge") {
                  if (is.element(abs(i), dg@visibleBlocks)) 
                    itemsBlockEdges[[-from]][[i]] <<- edgeNode
                  else NULL
                }
                else if (edge.type == "FactorEdge") 
                  itemsFactorEdges[[-from]][[i]] <<- edgeNode
                else itemsExtraEdges[[-from]][[i]] <<- edgeNode
            }
            "reinsertEdgeItem" <- function(edgeNode, from, to, 
                nr, edge.type = "VertexEdge") {
                edges <- edgeItem(from, edge.type = edge.type)
                if (length(edges) > 0) 
                  for (i in seq(along = edges)) {
                    e <- edges[[i]]
                    if (!(is.null(e))) 
                      if ((e$nr == nr) && (e$type == edge.type)) 
                        if (e$to == to) {
                          subSetEdgeItem(from, edge.type, i, 
                            edgeNode)
                        }
                  }
            }
            "setEdgeItem" <- function(i, edge.type = "VertexEdge", 
                edges = NULL) {
                if (i > 0) 
                  itemsEdges[[i]] <<- edges
                else if (edge.type == "BlockEdge") {
                  itemsBlockEdges[[-i]] <<- edges
                }
                else if (edge.type == "FactorEdge") 
                  itemsFactorEdges[[-i]] <<- edges
                else itemsExtraEdges[[-i]] <<- edges
            }
            "openBlockItem" <- function(i) return(itemsOpenBlocks[[i]])
            "setOpenBlockItem" <- function(i, blocks) itemsOpenBlocks[[i]] <<- blocks
            "closedBlockItem" <- function(i) {
                return(itemsClosedBlocks[[i]])
            }
            "setCloseVertex" <- function(i, value, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) if (vertex.type == "Vertex") {
                closedVertex[i] <<- value
                if (value) {
                  tkdelete(canvas, vertexItem(i)$tag)
                  updateBlockEdges()
                  updateCountBlockEdges <<- updateCountBlockEdgesMain
                }
                else {
                  vertexColor <- retVertexColor(i, vertex.type)
                  drawVertex(i, w = control$w, vertexcolor = control$vertexColor, 
                    vertex.type = "Vertex")
                  setVertexColor(i, color = vertexColor, vertex.type = vertex.type)
                }
            }
            "setClosedBlock" <- function(i, value, update = TRUE) {
                if (i > 0) {
                  closedBlock[i] <<- value
                  if (all(is.na(positionsClosedBlocks[i, ]))) 
                    positionsClosedBlocks[i, ] <<- apply(positionsBlocks[i, 
                      , ], 1, mean)
                  if (value) 
                    tkdelete(canvas, openBlockItem(i)$tag)
                  else tkdelete(canvas, closedBlockItem(i)$tag)
                  if (update) 
                    if ((updateCountBlockEdges < updateCountBlockEdgesMain)) {
                      updateBlockEdges()
                      updateCountBlockEdges <<- updateCountBlockEdgesMain
                    }
                }
            }
            "isInClosedBlock" <- function(i) {
                a <- blockList[[i]]@ancestors
                result <- FALSE
                if (length(a) > 1) {
                  a <- a[a != 0]
                  result <- any(closedBlock[a])
                }
                return(result)
            }
            "setHiddenBlock" <- function(i, value, update = TRUE) {
                if (i > 0) {
                  hiddenBlock[i] <<- value
                  if (value) {
                    if (all(is.na(positionsClosedBlocks[i, ]))) 
                      positionsClosedBlocks[i, ] <<- apply(positionsBlocks[i, 
                        , ], 1, mean)
                    if (closedBlock[i]) 
                      tkdelete(canvas, closedBlockItem(i)$tag)
                    else tkdelete(canvas, openBlockItem(i)$tag)
                  }
                }
            }
            "retStratum" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                "strata" <- function(i) if (i > 0) 
                  strataBlocks[i]
                else i
                if (vertex.type == "ClosedBlock") 
                  strataBlocks[abs(i)]
                else if (vertex.type == "Vertex") {
                  if (.IsEmpty(blockList)) 
                    strataVertices[i]
                  else if (blocksVertices[i] == 0) 
                    return(0)
                  else strata(blocksVertices[i])
                }
                else if (vertex.type == "Factor") {
                  if (.IsEmpty(blockList)) 
                    strataFactorVertices[-i]
                  else strata(blocksFactorVertices[-i])
                }
                else if (vertex.type == "Extra") {
                  if (.IsEmpty(blockList)) 
                    strataExtraVertices[abs(i)]
                  else strata(blocksExtraVertices[i])
                }
            }
            "retBlockIndex" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (vertex.type == "ClosedBlock") 
                  abs(i)
                else if (vertex.type == "Vertex") 
                  blocksVertices[i]
                else if (vertex.type == "Factor") 
                  blocksFactorVertices[-i]
                else if (vertex.type == "Extra") 
                  blocksExtraVertices[i]
            }
            "setBlockIndex" <- function(i, value, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                update <- FALSE
                if (permit.update.block.index) {
                  if (vertex.type == "Vertex") {
                    blocksVertices[i] <<- value
                    b <- closedBlock[value] || hiddenBlock[value]
                    if ((value > 0) && (b != closedVertex[i])) {
                      if (!(control$constrained || constrainedVertices[i])) {
                        if (b %in% dg@visibleBlocks) 
                          setCloseVertex(i, !closedVertex[i], 
                            vertex.type)
                        if (!closedVertex[i]) {
                          pos <- retVertexPos(i, vertex.type)
                          moveEdgesToVertex(pos, i, edge.type = "VertexEdge")
                        }
                      }
                      update <- TRUE
                    }
                  }
                  else if (vertex.type == "Factor") 
                    blocksFactorVertices[-i] <<- value
                  else if (vertex.type == "Extra") 
                    blocksExtraVertices[i] <<- value
                }
                return(update)
            }
            "updateVertexInBlock" <- function(i, k, visibleBefore, 
                visibleAfter) {
                if (TRUE) {
                  child <- namesVertices[i]
                  if (!visibleBefore) 
                    child <- tdv(child)
                  if (k == 0) 
                    parent = "root"
                  else {
                    parent <- blockLabels[k]
                    parent <- ubl(label = parent, index = k)
                  }
                  tkdelete(GW.top$env$box, child)
                  m <- 0
                  if (control$debug.strata) 
                    print(c(i, k))
                  for (j in seq(along = vertexList)) {
                    vertex <- vertexList[[j]]
                    stratum <- retStratum(j, vertex.type = "Vertex")
                    if (control$debug.strata && FALSE) {
                      STRATUM <- stratum(vertex)
                      a <- ""
                      if (stratum != STRATUM) 
                        a <- "%"
                      b <- ""
                      if (j != index(vertex)) 
                        b <- "#"
                      print(paste(c(a, b, name(vertex), index(vertex), 
                        as.numeric(stratum), names(STRATUM), 
                        as.numeric(STRATUM)), collapse = ", "))
                    }
                    if ((stratum == k) && (index(vertex) < i)) {
                      m <- m + 1
                      if (control$debug.strata) 
                        print(m)
                    }
                  }
                  child <- namesVertices[i]
                  fill <- "ForestGreen"
                  if (!visibleAfter) {
                    child <- tdv(child)
                    fill <- "LimeGreen"
                  }
                  tkinsert(GW.top$env$box, m, parent, child, 
                    text = child, fill = fill)
                }
            }
            "updateVertexBlockIndex" <- function(position, i) {
                currentIndex <- retBlockIndex(i, vertex.type = "Vertex")
                update <- FALSE
                if (permit.update.block.index) {
                  if (!.IsEmpty(blockList)) {
                    k <- 0
                    for (j in seq(along = blockList)) if (inBlock(position, 
                      j)) {
                      k <- j
                    }
                    change <- setBlockIndex(i, k, vertex.type = "Vertex")
                    update <- update || change
                  }
                }
                update <- update || (currentIndex != retBlockIndex(i, 
                  vertex.type = "Vertex"))
                if (control$variableFrame && update) {
                  if ((get("type", GW.top$env$box$env) == "variableList")) {
                  }
                  else {
                    v <- is.element(i, dg@visibleVertices)
                    updateVertexInBlock(i, k, visibleBefore = v, 
                      visibleAfter = v)
                  }
                }
                return(update)
            }
            "updateAllBlockIndices" <- function() {
                updateEdges <- FALSE
                if (permit.update.block.index) {
                  if (!is.null(vertexList)) 
                    for (i in seq(along = vertexList)) {
                      update <- updateVertexBlockIndex(positionsVertices[i, 
                        ], i)
                      updateEdges <- updateEdges || update
                    }
                  if (updateEdges) {
                    setUpdateBlockEdges("updateAllBlockIndices")
                  }
                }
                return(updateEdges)
            }
            "findMove" <- function(position, dxy = rep(0, local.N)) {
                return(inversProject(inversCanvasPosition(positionsCanvas(project(position)) + 
                  dxy)))
            }
            "findDifference" <- function(p1, p2) {
                return(relativePositionsCanvas(project(inversProject(inversCanvasPosition(p1)) - 
                  inversProject(inversCanvasPosition(p2)))))
            }
            "retVertexPos" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (vertex.type == "ClosedBlock") 
                  position <- positionsClosedBlocks[abs(i), ]
                else if (vertex.type == "Vertex") {
                  if (closedVertex[i]) 
                    position <- positionsClosedBlocks[blockReferences[retBlockIndex(i, 
                      vertex.type)], ]
                  else position <- positionsVertices[i, ]
                }
                else if (vertex.type == "Factor") 
                  position <- positionsFactorVertices[-i, ]
                else if (vertex.type == "Extra") 
                  position <- positionsExtraVertices[abs(i), 
                    ]
                return(positionsCanvas(project(position)))
            }
            "setVertexPos" <- function(i, xy, dxy = rep(0, local.N), 
                vertex.type = ifelse(i > 0, "Vertex", "Factor")) {
                ok <- TRUE
                position <- inversProject(inversCanvasPosition(xy))
                if (vertex.type == "ClosedBlock") 
                  positionsClosedBlocks[i, ] <<- position
                else if (vertex.type == "Vertex") {
                  old.position <- positionsVertices[i, ]
                  old.labelposition <- positionsLabels[i, ]
                  positionsVertices[i, ] <<- position
                  positionsLabels[i, ] <<- findMove(positionsLabels[i, 
                    ], dxy)
                  if (permit.update.block.index) {
                    if (updateVertexBlockIndex(position, i)) {
                      if (control$constrained || constrainedVertices[i]) {
                        positionsVertices[i, ] <<- old.position
                        positionsLabels[i, ] <<- old.labelposition
                        updateVertexBlockIndex(old.position, 
                          i)
                        ok <- FALSE
                      }
                      else setUpdateBlockEdges("setVertexPos")
                    }
                  }
                }
                else if (vertex.type == "Factor") {
                  positionsFactorVertices[-i, ] <<- position
                  positionsFactorLabels[-i, ] <<- findMove(positionsFactorLabels[-i, 
                    ], dxy)
                }
                else if (vertex.type == "Extra") {
                  positionsExtraVertices[abs(i), ] <<- position
                  positionsExtraLabels[abs(i), ] <<- findMove(positionsExtraLabels[abs(i), 
                    ], dxy)
                }
                return(ok)
            }
            "changeVertexPos" <- function(i, dxy = rep(0, local.N), 
                vertex.type = ifelse(i > 0, "Vertex", "Factor")) {
                ok <- TRUE
                if (vertex.type == "Vertex") {
                  old.position <- positionsVertices[i, ]
                  old.labelposition <- positionsLabels[i, ]
                  position <- findMove(positionsVertices[i, ], 
                    dxy)
                  positionsVertices[i, ] <<- position
                  positionsLabels[i, ] <<- findMove(positionsLabels[i, 
                    ], dxy)
                  if (updateVertexBlockIndex(position, i)) {
                    if (control$constrained || constrainedVertices[i]) {
                      positionsVertices[i, ] <<- old.position
                      positionsLabels[i, ] <<- old.labelposition
                      updateVertexBlockIndex(old.position, i)
                      ok <- FALSE
                    }
                  }
                }
                else if (vertex.type == "Factor") {
                  positionsFactorVertices[-i, ] <<- findMove(positionsFactorVertices[-i, 
                    ], dxy)
                  positionsFactorLabels[-i, ] <<- findMove(positionsFactorLabels[-i, 
                    ], dxy)
                }
                else if (vertex.type == "Extra") {
                  positionsExtraVertices[i, ] <<- findMove(positionsExtraVertices[i, 
                    ], dxy)
                  positionsExtraLabels[i, ] <<- findMove(positionsExtraLabels[i, 
                    ], dxy)
                }
                return(ok)
            }
            "retVertexName" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) if (vertex.type == "OpenBlock") 
                blockLabels[i]
            else if (vertex.type == "ClosedBlock") 
                blockLabels[abs(i)]
            else if (vertex.type == "Vertex") 
                namesVertices[i]
            else if (vertex.type == "Factor") 
                namesFactorVertices[-i]
            else if (vertex.type == "Extra") 
                extraLabels[abs(i)]
            "selectedNodesMatrix" <- function() {
                if (length(selectedNodes) > 0) {
                  r <- data.frame(index = unlist(lapply(selectedNodes, 
                    function(k) k$index)), hit = unlist(lapply(selectedNodes, 
                    function(k) k$hit.type)), type = unlist(lapply(selectedNodes, 
                    function(k) k$node.type)))
                  r
                }
            }
            "selectedEdgesMatrix" <- function() {
                if (length(selectedEdges) > 0) {
                  data.frame(index = unlist(lapply(selectedEdges, 
                    function(k) k$index)), from = unlist(lapply(selectedEdges, 
                    function(k) k$from)), to = unlist(lapply(selectedEdges, 
                    function(k) k$to)), hit = unlist(lapply(selectedEdges, 
                    function(k) k$hit.type)), type = unlist(lapply(selectedEdges, 
                    function(k) k$edge.type)))
                }
            }
            "retVertexColor" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (length(selectedNodes) > 0) {
                  x <- lapply(selectedNodes, function(k) ((i == 
                    k$index) && ("none" != k$hit.type) && (vertex.type == 
                    k$node.type)))
                  if (any(unlist(x))) 
                    return("YellowGreen")
                }
                if (vertex.type == "ClosedBlock") 
                  color(blockList[[i]])
                else if (vertex.type == "Vertex") 
                  colorsVertices[i]
                else if (vertex.type == "Factor") 
                  colorsFactorVertices[-i]
                else if (vertex.type == "Extra") 
                  colorsExtraVertices[abs(i)]
            }
            "setEdgeColor" <- function(i = 0, edge.type = "VertexEdge", 
                color = NULL) {
                activefill <- "LimeGreen"
                if (is.null(color)) 
                  activefill <- "DarkSlateGray"
                "f" <- function(x) if (is.null(x)) 
                  "NULL"
                else x
                E <- getEdges(edge.type = edge.type)[[i]]
                if (is.null(color)) 
                  color <- E@color
                "setEdge" <- function(k, edges) {
                  if (k != 0) {
                    edges <- edgeItem(k, edge.type = edge.type)
                    if (length(edges) > 0) 
                      for (e in edges) if (!(is.null(e))) 
                        if ((e$nr == i) && (e$type == edge.type)) 
                          if (e$to > k) {
                            for (l in 1:length(e$edges)) tkitemconfigure(canvas, 
                              e$edges[[l]], fill = color, activefill = activefill)
                            for (l in 1:length(e$tags)) tkitemconfigure(canvas, 
                              e$tags[[l]], fill = color, activefill = activefill)
                          }
                  }
                }
                for (j in E@vertex.indices) setEdge(j, edgeItem(j, 
                  edge.type = edge.type))
            }
            "setEdgeDash" <- function(i = 0, edge.type = "VertexEdge", 
                dash = NULL) {
                "f" <- function(x) if (is.null(x)) 
                  "NULL"
                else x
                E <- getEdges(edge.type = edge.type)[[i]]
                if (is.null(color)) 
                  color <- E@color
                "setEdge" <- function(k, edges) {
                  if (k != 0) {
                    edges <- edgeItem(k, edge.type = edge.type)
                    if (length(edges) > 0) 
                      for (e in edges) if (!(is.null(e))) 
                        if ((e$nr == i) && (e$type == edge.type)) 
                          if (e$to > k) {
                            for (l in 1:length(e$edges)) tkitemconfigure(canvas, 
                              e$edges[[l]], dash = dash)
                            for (l in 1:length(e$tags)) tkitemconfigure(canvas, 
                              e$tags[[l]], dash = dash)
                          }
                  }
                }
                for (j in E@vertex.indices) setEdge(j, edgeItem(j, 
                  edge.type = edge.type))
            }
            "setVertexColor" <- function(i, color = retVertexColor(i, 
                vertex.type), vertex.type = ifelse(i > 0, "Vertex", 
                "Factor"), permanent = FALSE) {
                if (!(length(color) == 1)) {
                  print(color)
                  print(i)
                }
                if (color == "cyan") 
                  activefill <- "DarkCyan"
                else activefill <- "LimeGreen"
                if (color == retVertexColor(i, vertex.type)) 
                  activefill <- "IndianRed"
                items <- vertexItem(i, vertex.type)$dot$dynamic
                if (!is.null(items)) 
                  if (length(items) > 0) 
                    for (k in seq(length(items))) tkitemconfigure(canvas, 
                      items[[k]], fill = color[[1]], activefill = activefill)
                if (permanent) {
                  if (vertex.type == "ClosedBlock") 
                    color(blockList[[i]]) <<- color
                  else if (vertex.type == "Vertex") 
                    colorsVertices[i] <<- color
                  else if (vertex.type == "Factor") 
                    colorsFactorVertices[abs(i)] <<- color
                  else if (vertex.type == "Extra") 
                    colorsExtraVertices[abs(i)] <<- color
                }
            }
            "retVertexLabel" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) if (control$useNamesForLabels) {
                if ((vertex.type == "OpenBlock") || (vertex.type == 
                  "ClosedBlock")) 
                  blockLabels[i]
                else if (vertex.type == "Vertex") 
                  namesVertices[i]
                else if (vertex.type == "Factor") 
                  namesFactorVertices[-i]
                else if (vertex.type == "Extra") 
                  namesExtraVertices[i]
            }
            else {
                if ((vertex.type == "OpenBlock") || (vertex.type == 
                  "ClosedBlock")) 
                  blockLabels[i]
                else if (vertex.type == "Vertex") 
                  Labels[i]
                else if (vertex.type == "Factor") 
                  factorLabels[-i]
                else if (vertex.type == "Extra") 
                  extraLabels[i]
            }
            "setVertexLabel" <- function(i, label, vertex.type) {
                if (vertex.type == "ClosedBlock") {
                  blockLabels[i] <<- label
                  tkitemconfigure(canvas, itemsClosedBlocks[[i]]$l, 
                    text = label)
                }
                else if (vertex.type == "Vertex") {
                  Labels[i] <<- label
                  tkitemconfigure(canvas, itemsVertices[[i]]$l, 
                    text = label)
                }
                else if (vertex.type == "Factor") {
                  factorLabels[-i] <<- label
                  tkitemconfigure(canvas, itemsFactors[[-i]]$l, 
                    text = label)
                }
                else if (vertex.type == "Extra") {
                  extraLabels[i] <<- label
                  tkitemconfigure(canvas, itemsExtras[[i]]$l, 
                    text = label)
                }
            }
            "retLabelPos" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                if (vertex.type == "ClosedBlock") 
                  position <- positionsBlockLabels[i, ]
                else if (vertex.type == "Vertex") 
                  position <- positionsLabels[i, ]
                else if (vertex.type == "Factor") 
                  position <- positionsFactorLabels[-i, ]
                else if (vertex.type == "Extra") 
                  position <- positionsExtraLabels[i, ]
                positionsCanvas(project(position))
            }
            "setLabelPos" <- function(i, xy, dxy = rep(0, local.N), 
                vertex.type = ifelse(i > 0, "Vertex", "Factor")) if (vertex.type == 
                "ClosedBlock") {
                positionsBlockLabels[i, ] <<- positionsBlockLabels[i, 
                  ] + inversCanvasRelativePosition(dxy)
            }
            else if (vertex.type == "Vertex") 
                positionsLabels[i, ] <<- findMove(positionsLabels[i, 
                  ], dxy)
            else if (vertex.type == "Factor") 
                positionsFactorLabels[-i, ] <<- findMove(positionsFactorLabels[-i, 
                  ], dxy)
            else if (vertex.type == "Extra") 
                positionsExtraLabels[i, ] <<- findMove(positionsExtraLabels[i, 
                  ], dxy)
            "retEdgeLabelPos" <- function(label.number, f = 0, 
                t = 0) relativePositionsCanvas(positionsEdgeLabels[label.number, 
                ])
            "setEdgeLabelPos" <- function(edgeNode, label.number, 
                xy, dxy = rep(0, local.N), f = 0, t = edgeNode$to, 
                edge.type = edgeNode$type) {
                positionsEdgeLabels[label.number, ] <<- positionsEdgeLabels[label.number, 
                  ] + inversCanvasRelativePosition(dxy)
                E <- getEdges(edge.type = edge.type)
                labelPosition(E[[edgeNode$nr]]) <- positionsEdgeLabels[label.number, 
                  ]
            }
            "setEdgeLabel" <- function(edgeNode, label = "", 
                i = edgeNode$label.number, f = 0, t = edgeNode$to, 
                edge.type = edgeNode$type, permanent = TRUE) {
                if (control$debug.strata && (label != " ")) {
                  text <- paste(paste(i, paste(f, t, sep = "-"), 
                    sep = "<"), label, sep = ">")
                  tkitemconfigure(canvas, edgeNode$label, text = text)
                  tkitemconfigure(canvas, edgeNode$label, fill = myColor(i + 
                    2))
                }
                else tkitemconfigure(canvas, edgeNode$label, 
                  text = label)
                if (permanent) {
                  if (edgeNode$type == "VertexEdge") 
                    label(GraphWindow@dg@edgeList[[edgeNode$nr]]) <<- label
                  else if (edgeNode$type == "FactorEdge") 
                    label(GraphWindow@dg@factorEdgeList[[edgeNode$nr]]) <<- label
                  else if (edgeNode$type == "ExtraEdge") 
                    label(GraphWindow@dg@extraEdgeList[[edgeNode$nr]]) <<- label
                  else if (edgeNode$type == "BlockEdge") 
                    label(GraphWindow@dg@blockEdgeList[[edgeNode$nr]]) <<- label
                }
            }
            "retEdgeLabel" <- function(edgeNode, i, f, t, edge.type = "VertexEdge") {
                if (is.na(control$namesOnEdges)) 
                  label <- ""
                else if (edgeNode$type == "VertexEdge") 
                  label(GraphWindow@dg@edgeList[[edgeNode$nr]])
                else if (edgeNode$type == "FactorEdge") 
                  label(GraphWindow@dg@factorEdgeList[[edgeNode$nr]])
                else if (edgeNode$type == "ExtraEdge") 
                  label(GraphWindow@dg@extraEdgeList[[edgeNode$nr]])
                else if (edgeNode$type == "BlockEdge") 
                  label(GraphWindow@dg@blockEdgeList[[edgeNode$nr]])
            }
            "setEdgeWidth" <- function(edgeNode, width = 1, i = edgeNode$label.number, 
                f = 0, t = edgeNode$to, edge.type = edgeNode$type) {
                for (l in 1:length(edgeNode$edges)) tkitemconfigure(canvas, 
                  edgeNode$edges[[l]], width = width)
                if (edgeNode$type == "VertexEdge") 
                  width(GraphWindow@dg@edgeList[[edgeNode$nr]]) <<- width
                else if (edgeNode$type == "FactorEdge") 
                  width(GraphWindow@dg@factorEdgeList[[edgeNode$nr]]) <<- width
                else if (edgeNode$type == "ExtraEdge") 
                  width(GraphWindow@dg@extraEdgeList[[edgeNode$nr]]) <<- width
                else if (edgeNode$type == "BlockEdge") 
                  width(GraphWindow@dg@blockEdgeList[[edgeNode$nr]]) <<- width
            }
            "vertexTypeOfEdge" <- function(index, edge.type = "VertexEdge", 
                edgeObject = NULL) if (edge.type == "factorBlockEdge") 
                ifelse(index > 0, "Factor", "ClosedBlock")
            else ifelse(index > 0, "Vertex", ifelse(edge.type == 
                "FactorEdge", "Factor", ifelse(edge.type == "ExtraEdge", 
                "Extra", "ClosedBlock")))
            "displayNode" <- function(i, type) {
                display <- TRUE
                if (type == "Factor") {
                }
                else if (type == "Extra") {
                }
                else if (type == "ClosedBlock") {
                  if (!closedBlock[abs(i)]) 
                    display <- FALSE
                  if (hiddenBlock[abs(i)]) 
                    display <- FALSE
                  if (!(abs(i) %in% dg@visibleBlocks)) 
                    display <- FALSE
                }
                else if (closedVertex[abs(i)]) 
                  display <- FALSE
                return(display)
            }
            "displayEdge" <- function(f, t, from.type = vertexTypeOfEdge(f, 
                edgeNode$type), to.type = vertexTypeOfEdge(t, 
                edgeNode$type), edgeNode = NULL) {
                display <- displayNode(f, from.type)
                if (display) 
                  display <- displayNode(t, to.type)
                return(display)
            }
            "setEdgeCoords" <- function(edgeNode, edge.type = "Edge", 
                f = edgeNode$from, t = edgeNode$to, posFrom = retVertexPos(f, 
                  from.type), posTo = retVertexPos(t, to.type), 
                from.type = vertexTypeOfEdge(f, edgeNode$type), 
                to.type = vertexTypeOfEdge(t, edgeNode$type), 
                raise = TRUE, setEdgeLabel = TRUE, width = control$w) {
                edgeNode <- subSetEdgeCoords(edgeNode, f, t, 
                  posFrom, posTo, from.type, to.type)
                l <- sqrt(sum((posTo - posFrom)^2))
                doDrawEdge <- FALSE
                if (is.null(edgeNode$edges)) 
                  doDrawEdge <- TRUE
                else {
                  if (l < 1) 
                    setEdgeLabel(edgeNode, " ", edgeNode$label.number, 
                      f = f, edge.type = edge.type, permanent = FALSE)
                  else {
                    if (raise) 
                      for (l in 1:length(edgeNode$edges)) tkitemraise(canvas, 
                        edgeNode$edges[[l]])
                    if (length(edgeNode$tags) > 0) 
                      for (l in 1:length(edgeNode$tags)) tkitemraise(canvas, 
                        edgeNode$tags[[l]])
                    tkitemraise(canvas, edgeNode$label)
                  }
                  if (setEdgeLabel) {
                    display <- displayEdge(f = f, t, edgeNode = edgeNode)
                    if (display) {
                      label <- retEdgeLabel(edgeNode, edgeNode$label.number, 
                        f = f)
                      setEdgeLabel(edgeNode, label, edgeNode$label.number, 
                        f = f, edge.type = edge.type, permanent = FALSE)
                    }
                  }
                  posLabel <- (posFrom + posTo)/2 + retEdgeLabelPos(edgeNode$label.number, 
                    f, t)
                  tkcoords(canvas, edgeNode$label, posLabel[1], 
                    posLabel[2])
                }
            }
            "subSetEdgeCoords" <- function(edgeNode, f, t, posFrom = retVertexPos(f, 
                from.type), posTo = retVertexPos(t, to.type), 
                from.type = vertexTypeOfEdge(f, edgeNode$type), 
                to.type = vertexTypeOfEdge(t, edgeNode$type), 
                width = control$w) {
                stratumFrom <- retStratum(f, from.type)
                stratumTo <- retStratum(t, to.type)
                reverse <- edgeNode$reverse
                display <- displayEdge(f, t, from.type, to.type, 
                  edgeNode = edgeNode)
                drawnEdge <- TRUE
                if (display) {
                  diff <- posTo - posFrom
                  l <- sqrt(sum(diff^2))
                  posTo <- posTo - diff * min(2 * control$w, 
                    l)/l
                  posFrom <- posFrom + diff * min(2 * control$w, 
                    l)/l
                  if (is.null(edgeNode$edges)) {
                    edges <- getEdges(edge.type = edgeNode$type)
                    drawEdge(edges[[edgeNode$nr]], edgeNode$nr, 
                      lower = TRUE, edge.type = edgeNode$type, 
                      reinsert = TRUE)
                    edgeNodes <- edgeItem(f, edge.type = edgeNode$type) # 'edge.type' = edgeNode$type ???
                    if (length(edgeNodes) > 0) 
                      for (e in edgeNodes) {
                        if (!(is.null(e))) 
                          if ((e$nr == edgeNode$nr) && (e$type == 
                            edgeNode$type)) 
                            if (e$to == edgeNode$to) 
                              edgeNode <- e
                      }
                  }
                  else {
                    label <- retEdgeLabel(edgeNode, edgeNode$label.number, 
                      f, t, edge.type = edgeNode$type)
                    tkitemconfigure(canvas, edgeNode$label, text = label)
                  }
                }
                else {
                  posTo <- c(0, 0)
                  posFrom <- c(0, 0)
                  if (is.null(edgeNode$edges)) 
                    drawnEdge <- FALSE
                  else {
                    tkitemconfigure(canvas, edgeNode$label, text = "")
                  }
                }
                if (drawnEdge) {
                  if ((stratumFrom != 0) || (stratumTo != 0) || 
                    !Oriented) {
                    if (is.na(edgeNode$oriented)) 
                      none <- stratumFrom == stratumTo
                    else none <- !edgeNode$oriented
                    for (l in 1:length(edgeNode$edges)) if (none) 
                      tkitemconfigure(canvas, edgeNode$edges[[l]], 
                        arrow = "none")
                    else tkitemconfigure(canvas, edgeNode$edges[[l]], 
                      arrow = "last")
                    edge.oriented <- FALSE
                    if (!is.na(edgeNode$oriented)) 
                      edge.oriented <- edgeNode$oriented
                    if (edge.oriented) 
                      reverse <- edgeNode$reverse
                    else reverse <- (stratumFrom > stratumTo)
                  }
                  "g" <- function(pos, ll) {
                    dxy <- tkcoords(canvas, edgeNode$tags[[ll]])
                    dxy <- apply(matrix(as.numeric(dxy), ncol = 2, 
                      byrow = 2), 2, mean)
                    dxy <- pos[1:2] - dxy
                    tkmove(canvas, edgeNode$tags[[ll]], dxy[1], 
                      dxy[2])
                  }
                  pos <- (posFrom + posTo)/2
                  if (length(edgeNode$tags) > 0) 
                    for (ll in 1:length(edgeNode$tags)) g(pos, 
                      ll)
                  "ff" <- function(posFrom, posTo, ll) {
                    if (reverse) 
                      tkcoords(canvas, edgeNode$edges[[ll]], 
                        posTo[1], posTo[2], posFrom[1], posFrom[2])
                    else tkcoords(canvas, edgeNode$edges[[ll]], 
                      posFrom[1], posFrom[2], posTo[1], posTo[2])
                  }
                  if (length(edgeNode$edges) == 1) 
                    ff(posFrom, posTo, 1)
                  else {
                    d <- posFrom - posTo
                    ld <- sqrt(sum(d[1:2]^2))
                    e <- width * d[1:2]/ld/2 * (length(edgeNode$edges) - 
                      1)/4
                    d <- width * c(-d[2], d[1])/ld/2 * (length(edgeNode$edges) - 
                      1)
                    if (length(edgeNode$edges) == 2) {
                      ff(posFrom[1:2] + d, posTo[1:2] + d, 1)
                      ff(posFrom[1:2] - d, posTo[1:2] - d, 2)
                    }
                    else {
                      for (lll in 1:length(edgeNode$edges)) {
                        kk <- lll - (length(edgeNode$edges) + 
                          1)/2
                        ff(posFrom[1:2] + kk * d + abs(kk) * 
                          e, posTo[1:2] + kk * d - abs(kk) * 
                          e, lll)
                      }
                    }
                  }
                }
                return(edgeNode)
            }
            "retBlockPos" <- function(i, j) {
                position <- positionsBlocks[i, , j]
                positionsCanvas(project(position))
            }
            "changeBlockCornerPos" <- function(i, A, dxy) {
                db <- toBlockPoints(A, dxy)
                positionsBlocks[i, , 1] <<- findMove(positionsBlocks[i, 
                  , 1], db[1, ])
                positionsBlocks[i, , 2] <<- findMove(positionsBlocks[i, 
                  , 2], db[2, ])
            }
            "changeBlockPos" <- function(i, A, dxy) {
                positionsBlocks[i, , 1] <<- findMove(positionsBlocks[i, 
                  , 1], dxy)
                positionsBlocks[i, , 2] <<- findMove(positionsBlocks[i, 
                  , 2], dxy)
            }
            "inBlock" <- function(position, block) {
                if (is.null(blockList[[block]])) 
                  return(FALSE)
                else {
                  block.position <- t(positionsBlocks[block, 
                    , ])
                  if (!all((block.position[1, ] < block.position[2, 
                    ]))) 
                    warning("Invalid block positions")
                  return(all((block.position[1, ] < position) & 
                    (position < block.position[2, ])))
                }
            }
            "retBlockPoints" <- function(i, header = FALSE, box = FALSE, 
                n) {
                A <- positionsBlocks[i, , 1]
                if (header) {
                  if (box) {
                    A <- A + c(1, 1, rep(0, local.N - 2))
                    B <- A + c(3 * n, 5, rep(0, local.N - 2))
                  }
                  else {
                    B <- positionsBlocks[i, , 2]
                    A <- A + c(1, 5, rep(0, local.N - 2))
                    B[1] <- B[1] - 1
                    B[2] <- A[2] + 1
                  }
                }
                else B <- positionsBlocks[i, , 2]
                delta <- c(0, 0, 0)
                position <- matrix(c(c(A[1], A[2], A[3]), c(B[1], 
                  B[2], B[3]) + delta, c(A[1], B[2], A[3]), c(B[1], 
                  A[2], A[3]), c(A[1], A[2], B[3]) + delta, c(A[1], 
                  B[2], B[3]) + delta, c(B[1], A[2], B[3]) + 
                  delta, c(B[1], B[2], A[3])), ncol = 3, byrow = TRUE)
                if (local.N < 3) 
                  position <- position[, 1:local.N]
                else if (local.N > 3) 
                  for (i in 4:local.N) position <- cbind(position, 
                    rep(0, 8))
                positionsCanvas(project(position))
            }
            "toBlockPoints" <- function(n, p) {
                result <- switch(EXPR = paste(n - 1), "0" = c(c(p[1], 
                  p[2], p[3]), c(0, 0, 0)), "1" = c(c(0, 0, 0), 
                  c(p[1], p[2], p[3])), "2" = c(c(p[1], 0, p[3]), 
                  c(0, p[2], 0)), "3" = c(c(0, p[2], p[3]), c(p[1], 
                  0, 0)), "4" = c(c(p[1], p[2], 0), c(0, 0, p[3])), 
                  "5" = c(c(p[1], 0, 0), c(0, p[2], p[3])), "6" = c(c(0, 
                    p[2], 0), c(p[1], 0, p[3])), "7" = c(c(0, 
                    0, p[3]), c(p[1], p[2], 0)))
                result <- matrix(result, ncol = 3, byrow = TRUE)
                if (local.N < 3) 
                  result <- result[, 1:local.N]
                if (local.N > 3) 
                  for (i in 4:local.N) result <- cbind(result, 
                    rep(0, 8))
                return(result)
            }
            "retBlockLabelPos" <- function(i) relativePositionsCanvas(positionsBlockLabels[i, 
                ])
            "addEdgePopups" <- function(canvas, edge, i, f, t, 
                edgePopupMenu, U.Menus, edge.type = "VertexEdge") {
                tkadd(edgePopupMenu, "command", label = paste("Edge from", 
                  retVertexLabel(f), "to", retVertexLabel(t), 
                  "(echo indices)"), command = function() {
                  print("Hej from edge")
                  print(c(f, t))
                })
                tkadd(edgePopupMenu, "command", label = paste("Delete edge (Here: Slave view!)"), 
                  accelerator = "[ double click edge ]", command = function() subDropEdge(i, 
                    f, t, edge.type = edge.type, slave = TRUE))
                tkadd(edgePopupMenu, "command", label = paste("Delete all edges to/from blocks"), 
                  command = function() subDropEdge(i, f, t, from.all = TRUE, 
                    to.all = TRUE, edge.type = edge.type, slave = FALSE))
                propEdgeMenu <- tkmenu(edgePopupMenu, tearoff = FALSE)
                tkadd(propEdgeMenu, "command", label = paste("Open dialog box for slot values"), 
                  command = function() propertyEdge(i, f, t, 
                    edge.type = edge.type)())
                tkadd(propEdgeMenu, "command", label = paste("/ Change edge class"), 
                  command = function() changeEdgeClass(i, f, 
                    t, edge.type = edge.type)())
                tkadd(propEdgeMenu, "command", label = paste("/ Set edge label"), 
                  command = function() {
                    activateEdge(i, from = f, to = t, edge.type = edge.type)()
                    changeEdgeLabel(i, f, t, edge.type = edge.type)()
                  })
                tkadd(propEdgeMenu, "command", label = paste("/ Compute edge label"), 
                  accelerator = "[ click label ]", command = function() {
                    activateEdge(i, from = f, to = t, edge.type = edge.type)()
                    computeEdgeLabel(i, f, t, FALSE, edge.type = edge.type)()
                  })
                tkadd(propEdgeMenu, "command", label = paste("/ Force computation of edge label"), 
                  accelerator = "[ double click label ]", command = function() {
                    activateEdge(i, from = f, to = t, edge.type = edge.type)()
                    computeEdgeLabel(i, f, t, TRUE, edge.type = edge.type)()
                  })
                tkadd(propEdgeMenu, "command", label = paste("/ Delete label of edge"), 
                  accelerator = "[ triple click label ]", command = function() deleteEdgeLabel(i, 
                    f, t, edge.type = edge.type)())
                tkadd(edgePopupMenu, "cascade", label = "Properties", 
                  menu = propEdgeMenu)
                helpEdgeMenu <- tkmenu(edgePopupMenu, tearoff = FALSE)
                tkadd(helpEdgeMenu, "command", label = paste(" - Add edge: Left click the vertices of the edge to add"), 
                  command = function() message("Left click the vertices of the edge to add"))
                tkadd(helpEdgeMenu, "command", label = paste(" - Drag edge: Move edge with two vertices"), 
                  command = function() message("Left click edge and drag edge"))
                tkadd(helpEdgeMenu, "command", label = paste(" - Drag label: Move label of edge"), 
                  command = function() message("Left click edge label and drag label"))
                tkadd(edgePopupMenu, "cascade", label = "Help on edges", 
                  menu = helpEdgeMenu)
                methEdgeMenu <- tkmenu(edgePopupMenu, tearoff = FALSE)
                if (hasMethod("addToPopups", class(edge))) 
                  addToPopups(edge, edge.type, methEdgeMenu, 
                    i, sinkView, Args)
                tkadd(edgePopupMenu, "cascade", label = "Items by method 'addToPopups'", 
                  menu = methEdgeMenu)
                userEdgeMenu <- tkmenu(edgePopupMenu, tearoff = FALSE)
                "UserEdgePopup" <- function(item) {
                  force(item)
                  force(f)
                  force(t)
                  force(edge.type)
                  force(edge)
                  function(...) {
                    sinkView(U.Menus[[item]])
                    j <- which.unordered.edge(c(t, f), edge.type = edge.type)
                    from.type <- vertexTypeOfEdge(f, edge.type, 
                      edge)
                    to.type <- vertexTypeOfEdge(t, edge.type, 
                      edge)
                    U.Menus[[item]]$command(object, retVertexName(f, 
                      from.type), retVertexName(t, to.type), 
                      from = f, to = t, from.type = from.type, 
                      to.type = to.type, edge.index = i, which.edge = j, 
                      edge.type = edge.type, Arguments = Args())
                  }
                }
                if (length(U.Menus) > 0) 
                  for (item in seq(along = U.Menus)) if (names(U.Menus[item]) == 
                    "Edge") 
                    tkadd(userEdgeMenu, "command", label = U.Menus[[item]]$label, 
                      command = UserEdgePopup(item))
                tkadd(edgePopupMenu, "cascade", label = "User defined items", 
                  menu = userEdgeMenu)
            }
            "setBindEdge" <- function(canvas, edge, line, label, 
                i, f, t, U.Menus, edge.type = "VertexEdge") {
                if (initial.set.popups) {
                  edgePopupMenu <- tkmenu(canvas, tearoff = FALSE)
                  addEdgePopups(canvas, edge, i, f, t, edgePopupMenu, 
                    U.Menus, edge.type)
                }
                tkitembind(canvas, label, "<Leave>", function() tkconfigure(canvas, 
                  cursor = "arrow"))
                tkitembind(canvas, label, "<Enter>", function() tkconfigure(canvas, 
                  cursor = "hand1"))
                tkitembind(canvas, line, "<Leave>", function() tkconfigure(canvas, 
                  cursor = "arrow"))
                tkitembind(canvas, line, "<Enter>", function() tkconfigure(canvas, 
                  cursor = "tcross"))
                tkitembind(canvas, label, "<Button-1>", activateEdge(i, 
                  from = f, to = t, edge.type = edge.type))
                tkitembind(canvas, label, "<B1-Motion>", moveEdgeLabel(i, 
                  f, t, edge.type = edge.type))
                tkitembind(canvas, label, "<ButtonRelease-1>", 
                  computeEdgeLabel(i, f, t, FALSE, edge.type = edge.type))
                tkitembind(canvas, label, "<Double-Button-1>", 
                  computeEdgeLabel(i, f, t, TRUE, edge.type = edge.type))
                tkitembind(canvas, label, "<Triple-Button-1>", 
                  deleteEdgeLabel(i, f, t, edge.type = edge.type))
                tkitembind(canvas, label, "<Shift-1>", changeEdgeClass(i, 
                  f, t, edge.type = edge.type))
                if (initial.set.popups) 
                  tkitembind(canvas, label, "<Button-3>", callPopup(i, 
                    edgePopupMenu))
                else tkitembind(canvas, label, "<Button-3>", 
                  callPopupEdge(edge, i, f, t, edge.type, U.Menus))
                tkitembind(canvas, line, "<Option-1>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "option-1", 
                  color = "DarkGreen"))
                tkitembind(canvas, line, "<Shift-1>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "shift-1", 
                  color = "DarkGreen"))
                tkitembind(canvas, line, "<Control-1>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "control-1", 
                  color = "SeaGreen"))
                tkitembind(canvas, line, "<Shift-Control-1>", 
                  activateEdge(i, from = f, to = t, edge.type, 
                    hit.type = "shift-control-1", color = "LightSeaGreen"))
                tkitembind(canvas, line, "<Option-3>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "option-3", 
                  color = "LightGreen"))
                tkitembind(canvas, line, "<Shift-3>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "shift-3", 
                  color = "LightGreen"))
                tkitembind(canvas, line, "<Control-3>", activateEdge(i, 
                  from = f, to = t, edge.type, hit.type = "control-3", 
                  color = "SpringGreen"))
                tkitembind(canvas, line, "<Shift-Control-3>", 
                  activateEdge(i, from = f, to = t, edge.type, 
                    hit.type = "shift-control-3", color = "LimeGreen"))
                tkitembind(canvas, line, "<Button-1>", activateEdge(i, 
                  from = f, to = t, edge.type = edge.type))
                tkitembind(canvas, line, "<Double-Button-1>", 
                  deleteEdge(i, f, t, edge.type = edge.type))
                tkitembind(canvas, line, "<B1-Motion>", moveEdge(i, 
                  f, t, edge.type = edge.type))
                if (initial.set.popups) 
                  tkitembind(canvas, line, "<Button-3>", callPopup(i, 
                    edgePopupMenu))
                else tkitembind(canvas, line, "<Button-3>", callPopupEdge(edge, 
                  i, f, t, edge.type, U.Menus))
            }
            "subDrawEdge" <- function(edge, i, edgecolor = "black", 
                lower = FALSE, edge.type = "VertexEdge", newE = FALSE) {
                "emptyEdge" <- function() list(list(lines = NULL, 
                  tags = NULL, from = f, to = t, label = NULL, 
                  label.position = NULL))
                type.negative <- ifelse(edge.type == "BlockEdge", 
                  "ClosedBlock", ifelse(edge.type == "FactorEdge", 
                    "Factor", "Extra"))
                useMethod <- FALSE
                result <- FALSE
                if (!is.null(edge)) 
                  useMethod <- hasMethod("draw", class(edge))
                if (useMethod) {
                  position <- edge.positions(i, type.negative = type.negative, 
                    edge.type = edge.type)
                  if (length(position) == 2) {
                    f <- from(i, edge.type = edge.type)
                    t <- to(i, edge.type = edge.type)
                    from.type <- vertexTypeOfEdge(f, edge.type, 
                      edge)
                    to.type <- vertexTypeOfEdge(t, edge.type, 
                      edge)
                    display <- displayEdge(f, t, from.type, to.type)
                  }
                  else {
                    display <- TRUE
                  }
                  if (display) {
                    strata <- edge.strata(i, type.negative = type.negative, 
                      edge.type = edge.type)
                    x <- lapply(position, function(e) e[1])
                    y <- lapply(position, function(e) e[2])
                    result <- draw(edge, canvas, position, x, 
                      y, stratum = strata, w = edge@width * Scale, 
                      color = edge@color, font.edge.label = font.edge.label, 
                      background = control$background)
                  }
                  else result <- emptyEdge()
                }
                else {
                  f <- from(i, edge.type = edge.type)
                  t <- to(i, edge.type = edge.type)
                  from.type <- vertexTypeOfEdge(f, edge.type, 
                    edge)
                  to.type <- vertexTypeOfEdge(t, edge.type, edge)
                  display <- displayEdge(f, t, from.type, to.type)
                  if (!display) {
                    print("Drawing edge to undisplay")
                    result <- emptyEdge
                  }
                  else {
                    posFrom <- retVertexPos(f, from.type)
                    posTo <- retVertexPos(t, to.type)
                    stratumFrom <- retStratum(f, from.type)
                    stratumTo <- retStratum(t, to.type)
                    if (stratumFrom == stratumTo) 
                      arrowhead = "none"
                    else if (stratumFrom < stratumTo) 
                      arrowhead = "last"
                    else arrowhead = "first"
                    E <- getEdges(edge.type = edge.type)[[i]]
                    line <- tkcreate(canvas, "line", posFrom[1], 
                      posFrom[2], posTo[1], posTo[2], arrow = arrowhead, 
                      width = E@width, fill = E@color)
                    label.position <- (posFrom + posTo)/2
                    pos <- label.position + rep(0, local.N)
                    txt <- E@label
                    label <- tkcreate(canvas, "text", pos[1], 
                      pos[2], text = txt, anchor = "nw", font = font.edge.label, 
                      activefill = "DarkSlateGray")
                    result <- list(list(lines = list(line), from = f, 
                      to = t, label = label, label.position = label.position))
                  }
                }
                return(result)
            }
            "insertEdgeItems" <- function(R, edge, i, edge.type = "VertexEdge", 
                newE = FALSE, reinsert = FALSE) {
                if (!is.null(R)) {
                  for (k in 1:length(R)) {
                    if (!reinsert) 
                      positionsEdgeLabels <<- rbind(positionsEdgeLabels, 
                        rep(0, local.N))
                    f <- R[[k]]$from
                    t <- R[[k]]$to
                    edge.oriented <- NA
                    if (is.element("oriented", slotNames(edge))) 
                      edge.oriented <- edge@oriented
                    if (is.null(R[[k]]$lines)) {
                      edgeNode <- list(nr = i, type = edge.type, 
                        to = t, tag = NULL, oriented = edge.oriented, 
                        reverse = FALSE, edges = NULL, tags = NULL, 
                        label = NULL, label.number = nrow(positionsEdgeLabels))
                      setEdgeItem(f, edge.type = edge.type, c(edgeItem(f, 
                        edge.type = edge.type), list(edgeNode)))
                      edgeNode <- list(nr = i, type = edge.type, 
                        to = f, tag = NULL, oriented = edge.oriented, 
                        reverse = TRUE, edges = NULL, tags = NULL, 
                        label = NULL, label.number = nrow(positionsEdgeLabels))
                      setEdgeItem(t, edge.type = edge.type, c(edgeItem(t, 
                        edge.type = edge.type), list(edgeNode)))
                    }
                    else {
                      if (length(R) > 1) 
                        tag <- getTag(edge.type, round(i + (k - 
                          1)/length(R), digits = 2))
                      else tag <- getTag(edge.type, i)
                      tkaddtag(canvas, tag, "withtag", R[[k]]$label)
                      for (l in 1:length(R[[k]]$lines)) tkaddtag(canvas, 
                        tag, "withtag", R[[k]]$lines[[l]])
                      if (length(R[[k]]$tags) > 0) 
                        for (l in 1:length(R[[k]]$tags)) tkaddtag(canvas, 
                          tag, "withtag", R[[k]]$tags[[l]])
                      edgeNode <- list(nr = i, type = edge.type, 
                        to = t, tag = tag, oriented = edge.oriented, 
                        reverse = FALSE, edges = R[[k]]$lines, 
                        tags = R[[k]]$tags, label = R[[k]]$label, 
                        label.number = nrow(positionsEdgeLabels))
                      if (reinsert) {
                        reinsertEdgeItem(edgeNode, from = f, 
                          to = t, nr = i, edge.type = edge.type)
                      }
                      else setEdgeItem(f, edge.type = edge.type, 
                        c(edgeItem(f, edge.type = edge.type), 
                          list(edgeNode)))
                      subSetEdgeCoords(edgeNode, f, t, width = control$w)
                      if (newE) 
                        tkitemconfigure(canvas, R[[k]]$label, 
                          text = edge@label)
                      edgeNode <- list(nr = i, type = edge.type, 
                        to = f, tag = tag, oriented = edge.oriented, 
                        reverse = TRUE, edges = R[[k]]$lines, 
                        tags = R[[k]]$tags, label = R[[k]]$label, 
                        label.number = nrow(positionsEdgeLabels))
                      if (reinsert) {
                        reinsertEdgeItem(edgeNode, from = t, 
                          to = f, nr = i, edge.type = edge.type)
                      }
                      else setEdgeItem(t, edge.type = edge.type, 
                        c(edgeItem(t, edge.type = edge.type), 
                          list(edgeNode)))
                      for (l in 1:length(R[[k]]$lines)) setBindEdge(canvas, 
                        edge, R[[k]]$lines[[l]], R[[k]]$label, 
                        i, f, t, control$UserMenus, edge.type = edge.type)
                      if (length(R[[k]]$tags) > 0) 
                        for (l in 1:length(R[[k]]$tags)) setBindEdge(canvas, 
                          edge, R[[k]]$tags[[l]], R[[k]]$label, 
                          i, f, t, control$UserMenus, edge.type = edge.type)
                    }
                  }
                }
            }
            "drawEdge" <- function(edge, i, edgecolor = "black", 
                lower = FALSE, edge.type = "VertexEdge", newE = FALSE, 
                reinsert = FALSE) {
                if (control$debug.edges) 
                  print(paste("drawEdge", i, edge.type))
                if (!any(nodeIndices(edge) == 0)) {
                  result <- subDrawEdge(edge, i, edgecolor, lower, 
                    edge.type, newE)
                  insertEdgeItems(result, edge, i, edge.type = edge.type, 
                    newE = newE, reinsert = reinsert)
                }
            }
            "tkcoordsBlock" <- function(i, color = "black", lower = FALSE) {
                "tkcoordsRectangleLine" <- function(line, i, A, 
                  B, positions, color = "black", width = 1) {
                  posA <- positions[A, ]
                  posB <- positions[B, ]
                  tkcoords(canvas, line, posA[1], posA[2], posB[1], 
                    posB[2])
                }
                "tkcoordsCornerLine" <- function(line, i, A, posA, 
                  posB, color = "black", width = 1) {
                  diff <- posB - posA
                  l <- sqrt(sum(diff^2))
                  posB <- posA + diff * min(30, l)/l
                  posA <- posA - diff * min(control$w/2, l)/l
                  tkcoords(canvas, line, posA[1], posA[2], posB[1], 
                    posB[2])
                }
                "tkcoordsRectangleCorner" <- function(line, i, 
                  A, B, C, D, positions, color = "black", width = 2) {
                  posA <- positions[A, ]
                  tkcoordsCornerLine(line[[1]], i, A, posA, positions[B, 
                    ], color, width)
                  tkcoordsCornerLine(line[[2]], i, A, posA, positions[C, 
                    ], color, width)
                  if (!is.null(transformation)) 
                    tkcoordsCornerLine(line[[3]], i, A, posA, 
                      positions[D, ], color, width)
                }
                "tkcoordsRectangle" <- function(rectangle, i, 
                  positions, color = "black", width = 1) {
                  line <- rectangle$Lines
                  tkcoordsRectangleLine(line[[1]], i, 1, 3, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[2]], i, 4, 8, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[3]], i, 1, 4, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[4]], i, 3, 8, positions, 
                    color, width)
                  if (!is.null(transformation)) {
                    tkcoordsRectangleLine(line[[5]], i, 5, 6, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[6]], i, 7, 2, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[7]], i, 5, 7, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[8]], i, 6, 2, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[9]], i, 1, 5, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[10]], i, 3, 6, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[11]], i, 4, 7, 
                      positions, color, width)
                    tkcoordsRectangleLine(line[[12]], i, 8, 2, 
                      positions, color, width)
                  }
                  corner <- rectangle$Corners
                  tkcoordsRectangleCorner(corner[[1]], i, 1, 
                    3, 4, 5, positions, color, width + 2)
                  tkcoordsRectangleCorner(corner[[2]], i, 8, 
                    4, 3, 2, positions, color, width + 2)
                  tkcoordsRectangleCorner(corner[[3]], i, 4, 
                    1, 8, 7, positions, color, width + 2)
                  tkcoordsRectangleCorner(corner[[4]], i, 3, 
                    8, 1, 6, positions, color, width + 2)
                  if (!is.null(transformation)) {
                    tkcoordsRectangleCorner(corner[[5]], i, 5, 
                      6, 7, 1, positions, color, width + 2)
                    tkcoordsRectangleCorner(corner[[6]], i, 2, 
                      7, 6, 8, positions, color, width + 2)
                    tkcoordsRectangleCorner(corner[[7]], i, 7, 
                      5, 2, 4, positions, color, width + 2)
                    tkcoordsRectangleCorner(corner[[8]], i, 6, 
                      2, 5, 3, positions, color, width + 2)
                  }
                }
                "tkcoordsBar" <- function(line, i, positions, 
                  color = "black", width = 1) {
                  tkcoordsRectangleLine(line[[1]], i, 1, 3, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[2]], i, 4, 8, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[3]], i, 1, 4, positions, 
                    color, width)
                  tkcoordsRectangleLine(line[[4]], i, 3, 8, positions, 
                    color, width)
                }
                positions <- retBlockPoints(i)
                if (!is.null(openBlockItem(i)$canvas)) 
                  tkcoords(canvas, openBlockItem(i)$canvas, positions[1, 
                    1], positions[1, 2], positions[8, 1], positions[8, 
                    2])
                if (!is.null(openBlockItem(i)$rectangle)) 
                  tkcoordsRectangle(openBlockItem(i)$rectangle, 
                    i, positions, color = color, width = 1)
                txt <- blockList[[i]]@label
                positions <- retBlockPoints(i, header = TRUE, 
                  n = nchar(txt))
                if (!is.null(openBlockItem(i)$bar)) 
                  tkcoordsBar(openBlockItem(i)$bar, i, positions, 
                    color = color, width = 1)
                pos <- retBlockPos(i, 1) + c(8, 4, rep(0, local.N - 
                  2))
                if (!is.null(openBlockItem(i)$label)) 
                  tkcoords(canvas, openBlockItem(i)$label, pos[1], 
                    pos[2])
            }
            "drawBlock" <- function(block, i, color = "Grey", 
                box = FALSE, lower = TRUE, setTag = TRUE) {
                "drawRectangleLine" <- function(i, A, B, positions, 
                  tag, color = "black", width = 1) {
                  posA <- positions[A, ]
                  posB <- positions[B, ]
                  line <- tkcreate(canvas, "line", posA[1], posA[2], 
                    posB[1], posB[2], width = width, fill = color)
                  tkaddtag(canvas, tag, "withtag", line)
                  tkitembind(canvas, line, "<B1-Motion>", moveBlockLine(i, 
                    A, B))
                  tkitembind(canvas, line, "<Leave>", function() tkconfigure(canvas, 
                    cursor = "arrow"))
                  if ((A == 1) && (B == 3)) 
                    cursor <- "left_side"
                  else if ((A == 4) && (B == 8)) 
                    cursor <- "right_side"
                  else if ((A == 1) && (B == 4)) 
                    cursor <- "top_side"
                  else if ((A == 3) && (B == 8)) 
                    cursor <- "bottom_side"
                  else if ((A == 5) && (B == 6)) 
                    cursor <- "left_side"
                  else if ((A == 7) && (B == 2)) 
                    cursor <- "right_side"
                  else if ((A == 5) && (B == 7)) 
                    cursor <- "top_side"
                  else if ((A == 6) && (B == 2)) 
                    cursor <- "bottom_side"
                  else if ((A == 1) && (B == 5)) 
                    cursor <- "top_side"
                  else if ((A == 3) && (B == 6)) 
                    cursor <- "bottom_side"
                  else if ((A == 4) && (B == 7)) 
                    cursor <- "top_side"
                  else if ((A == 8) && (B == 2)) 
                    cursor <- "bottom_side"
                  tkitembind(canvas, line, "<Enter>", function() tkconfigure(canvas, 
                    cursor = cursor))
                  return(line)
                }
                "drawCornerLine" <- function(i, A, posA, posB, 
                  tag, color = "black", width = 1) {
                  diff <- posB - posA
                  l <- sqrt(sum(diff^2))
                  posB <- posA + diff * min(25, l)/l
                  posA <- posA - diff * min(control$w/2, l)/l
                  line <- tkcreate(canvas, "line", posA[1], posA[2], 
                    posB[1], posB[2], width = width, fill = color)
                  tkaddtag(canvas, tag, "withtag", line)
                  tkitembind(canvas, line, "<B1-Motion>", moveBlockPoint(i, 
                    A))
                  tkitembind(canvas, line, "<Leave>", function() tkconfigure(canvas, 
                    cursor = "arrow"))
                  if ((A == 1) || (A == 5)) 
                    cursor <- "top_left_corner"
                  else if ((A == 8) || (A == 2)) 
                    cursor <- "bottom_right_corner"
                  else if ((A == 4) || (A == 7)) 
                    cursor <- "top_right_corner"
                  else if ((A == 3) || (A == 6)) 
                    cursor <- "bottom_left_corner"
                  tkitembind(canvas, line, "<Enter>", function() tkconfigure(canvas, 
                    cursor = cursor))
                  return(line)
                }
                "drawRectangleCorner" <- function(i, A, B, C, 
                  D, positions, tag, color = "black", width = 2) {
                  posA <- positions[A, ]
                  l <- vector("list", 3)
                  l[[1]] <- drawCornerLine(i, A, posA, positions[B, 
                    ], tag, color, width)
                  l[[2]] <- drawCornerLine(i, A, posA, positions[C, 
                    ], tag, color, width)
                  if (!is.null(transformation)) 
                    l[[3]] <- drawCornerLine(i, A, posA, positions[D, 
                      ], tag, color, width)
                  return(l)
                }
                "drawRectangle" <- function(i, positions, tag, 
                  color = "black", width = 1) {
                  l <- vector("list", 12)
                  l[[1]] <- drawRectangleLine(i, 1, 3, positions, 
                    tag, color, width)
                  l[[2]] <- drawRectangleLine(i, 4, 8, positions, 
                    tag, color, width)
                  l[[3]] <- drawRectangleLine(i, 1, 4, positions, 
                    tag, color, width)
                  l[[4]] <- drawRectangleLine(i, 3, 8, positions, 
                    tag, color, width)
                  if (!is.null(transformation)) {
                    l[[5]] <- drawRectangleLine(i, 5, 6, positions, 
                      tag, color, width)
                    l[[6]] <- drawRectangleLine(i, 7, 2, positions, 
                      tag, color, width)
                    l[[7]] <- drawRectangleLine(i, 5, 7, positions, 
                      tag, color, width)
                    l[[8]] <- drawRectangleLine(i, 6, 2, positions, 
                      tag, color, width)
                    l[[9]] <- drawRectangleLine(i, 1, 5, positions, 
                      tag, color, width)
                    l[[10]] <- drawRectangleLine(i, 3, 6, positions, 
                      tag, color, width)
                    l[[11]] <- drawRectangleLine(i, 4, 7, positions, 
                      tag, color, width)
                    l[[12]] <- drawRectangleLine(i, 8, 2, positions, 
                      tag, color, width)
                  }
                  c <- vector("list", 8)
                  c[[1]] <- drawRectangleCorner(i, 1, 3, 4, 5, 
                    positions, tag, color, width + 2)
                  c[[2]] <- drawRectangleCorner(i, 8, 4, 3, 2, 
                    positions, tag, color, width + 2)
                  c[[3]] <- drawRectangleCorner(i, 4, 1, 8, 7, 
                    positions, tag, color, width + 2)
                  c[[4]] <- drawRectangleCorner(i, 3, 8, 1, 6, 
                    positions, tag, color, width + 2)
                  if (!is.null(transformation)) {
                    c[[5]] <- drawRectangleCorner(i, 5, 6, 7, 
                      1, positions, tag, color, width + 2)
                    c[[6]] <- drawRectangleCorner(i, 2, 7, 6, 
                      8, positions, tag, color, width + 2)
                    c[[7]] <- drawRectangleCorner(i, 7, 5, 2, 
                      4, positions, tag, color, width + 2)
                    c[[8]] <- drawRectangleCorner(i, 6, 2, 5, 
                      3, positions, tag, color, width + 2)
                  }
                  return(list(Lines = l, Corners = c))
                }
                "drawBar" <- function(i, positions, tag, box = FALSE, 
                  color = "black", width = 1) {
                  l <- vector("list", 4)
                  l[[1]] <- drawRectangleLine(i, 1, 3, positions, 
                    tag, color, width)
                  l[[2]] <- drawRectangleLine(i, 4, 8, positions, 
                    tag, color, width)
                  l[[3]] <- drawRectangleLine(i, 1, 4, positions, 
                    tag, color, width)
                  l[[4]] <- drawRectangleLine(i, 3, 8, positions, 
                    tag, color, width)
                  return(l)
                }
                tag <- getTag("block", i, setTag = setTag)
                positions <- retBlockPoints(i)
                posA <- positions[1, ]
                posB <- positions[8, ]
                popupitems <- NULL
                blockcanvas <- NULL
                if (control$drawBlockBackground) 
                  if (is.null(transformation)) {
                    blockcanvas <- tkcreate(canvas, "rectangle", 
                      posA[1], posA[2], posB[1], posB[2], fill = color(block))
                    tkaddtag(canvas, tag, "withtag", blockcanvas)
                    popupitems <- append(popupitems, list(blockcanvas))
                  }
                if (control$drawBlockFrame) 
                  Rectangle <- drawRectangle(i, positions, tag, 
                    color = color, width = 1)
                else Rectangle <- NULL
                txt <- blockLabels[i]
                positions <- retBlockPoints(i, header = TRUE, 
                  box = box, n = nchar(txt))
                if (control$drawBlockFrame) {
                  Bar <- drawBar(i, positions, tag, box = box, 
                    color = color, width = 2)
                  popupitems <- append(popupitems, Bar)
                }
                else Bar <- NULL
                posA <- retBlockPos(i, 1)
                pos <- posA + c(8, 4, rep(0, local.N - 2))
                label <- tkcreate(canvas, "text", pos[1], pos[2], 
                  text = txt, anchor = "nw", font = font.block, 
                  activefill = "DarkSlateGray")
                setOpenBlockItem(i, list(tag = tag, rectangle = Rectangle, 
                  canvas = blockcanvas, bar = Bar, label = label, 
                  block = i))
                setBindNode(canvas, blockList[[i]], tag, popupitems, 
                  label, i, "OpenBlock", control$UserMenus)
            }
            "addNodePopups" <- function(vertex, i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor"), nodePopupMenu, UserNodePopupItems, 
                slave = TRUE) {
                label <- retVertexLabel(i, vertex.type)
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(nodePopupMenu, "command", label = paste("Vertex", 
                    label, "(echo index)"), command = function() print(paste("Hej from vertex", 
                    label, "with index", i)))
                else if ((vertex.type == "OpenBlock") || (vertex.type == 
                  "ClosedBlock")) 
                  tkadd(nodePopupMenu, "command", label = paste("Block", 
                    label, "(echo index)"), command = function() print(paste("Hej from block", 
                    label, "with index", i)))
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(nodePopupMenu, "command", label = paste("Highlight for adding edge"), 
                    accelerator = "[ Click vertex ]", command = function() {
                      subActivateVertex(i, color = "green", vertex.type = vertex.type)
                      message("Click the other vertex")
                    })
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(nodePopupMenu, "command", label = paste("Highlight block for adding edges"), 
                    accelerator = "[ Click block ]", command = function() {
                      subActivateVertex(i, color = "green", vertex.type = "ClosedBlock")
                      message("Click vertex or block")
                    })
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(nodePopupMenu, "command", 
                        label = paste("Add edge after highlighting with selecting class"), 
                    command = newEdge(i, vertex.type = "Vertex", 
                      slave = FALSE, selectClass = TRUE))
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(nodePopupMenu, "command", label = paste("Add edge after highlighting", 
                    if (slave) 
                      "(Here: Slave view!)"
                    else "", collapse = ""), accelerator = "[ Click vertex ]", 
                    command = newEdge(i, vertex.type = "Vertex", 
                      slave = slave))
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(nodePopupMenu, "command", 
                        label = paste("Adds edges from/to block after highlight", 
                    if (slave) 
                      "(Here: Slave view!)"
                    else "", collapse = ""), accelerator = "[ Click block ]", 
                    command = newEdge(i, vertex.type = "ClosedBlock", 
                      slave = slave))
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) {
                }
                else if ((vertex.type == "OpenBlock")) {
                  tkadd(nodePopupMenu, "command", label = paste("Mark vertices of block"), 
                    command = function() {
                      markVerticesOfBlock(i, descendants = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste(" --- and descendants"), 
                    command = function() {
                      markVerticesOfBlock(i, slave = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste("Undisplay frames of block"), 
                    command = function() {
                      undisplayBlock(i, descendants = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste(" --- and descendants"), 
                    command = function() {
                      undisplayBlock(i, slave = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste("'Delete' (close) block"), 
                    command = function() {
                      removeBlock(i, descendants = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste(" --- and descendants"), 
                    command = function() {
                      removeBlock(i, slave = FALSE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste("New sub block"), 
                    accelerator = "(~ <F4>)", command = function() {
                      position <- position(blockList[[i]])
                      position <- apply(position, 1, mean)
                      position <- matrix(c(position - 10, position + 
                        10), nrow = 2, byrow = TRUE)
                      new.Block(position, get.name = TRUE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste("Minimize (shade) block"), 
                    accelerator = "[ Double click block head ]", 
                    command = function() {
                      closeBlock(i)()
                    })
                  tkadd(nodePopupMenu, "command", label = paste("Maximize block"), 
                    command = function() {
                      zoomPositions <- positionsBlocks[i, , ]
                      zoomPositions[, 1] <- zoomPositions[, 1] - 
                        2
                      zoomPositions[, 2] <- zoomPositions[, 2] + 
                        2
                      zoomPositions <<- zoomPositions
                      subUpdateGraphWindow("Maximize", redrawVertices = TRUE, 
                        all.blockframes = TRUE)
                    })
                  tkadd(nodePopupMenu, "command", label = paste("Redraw full graph"), 
                    command = function() {
                      zoomPositions <<- NULL
                      subUpdateGraphWindow("Redraw", redrawVertices = TRUE, 
                        all.blockframes = TRUE)
                    })
                }
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(nodePopupMenu, "command", label = paste("Open block"), 
                    accelerator = "[ Double click minimized block ]", 
                    command = function() {
                      openBlock(i)()
                    })
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) {
                  if (!(is.element(i, returnVisibleVertices()))) 
                    tkadd(nodePopupMenu, "command", label = paste("Add vertex", 
                      if (slave) 
                        "(Here: Slave view!)"
                      else "", collapse = ""), command = function() subAddVertex(i, 
                      vertex.type = vertex.type, slave = slave))
                  else tkadd(nodePopupMenu, "command", label = paste("Delete vertex", 
                    if (slave) 
                      "(Here: Slave view!)"
                    else "", collapse = ""), command = function() subDropVertex(i, 
                    vertex.type = vertex.type, slave = slave))
                }
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) {
                }
                propNodeMenu <- tkmenu(nodePopupMenu, tearoff = FALSE)
                tkadd(propNodeMenu, "command", label = paste("Open dialog box for slot values"), 
                  command = function() propertyNode(i, vertex.type = vertex.type)())
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Change label"), 
                    accelerator = "[ Double click label ]", command = changeVertexLabel(i, 
                      vertex.type = vertex.type))
                else if ((vertex.type == "OpenBlock")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Change label of block"), 
                    command = changeVertexLabel(i, vertex.type = "ClosedBlock"))
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Change label"), 
                    accelerator = "[ Double click label of minimized block ]", 
                    command = changeVertexLabel(i, vertex.type = "ClosedBlock"))
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Delete vertex label"), 
                    command = deleteVertexLabel(i, vertex.type = "Vertex"))
                else if ((vertex.type == "OpenBlock")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Delete label of block"), 
                    command = deleteVertexLabel(i, vertex.type = "ClosedBlock"))
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(propNodeMenu, "command", label = paste("/ Delete label of block"), 
                    command = deleteVertexLabel(i, vertex.type = "ClosedBlock"))
                tkadd(nodePopupMenu, "cascade", label = "Properties", 
                  menu = propNodeMenu)
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) {
                }
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) {
                }
                helpNodeMenu <- tkmenu(nodePopupMenu, tearoff = FALSE)
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(helpNodeMenu, "command", label = paste(" - Delete vertex"), 
                    accelerator = "[ Double click vertex ]", 
                    command = function() {
                    })
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) {
                }
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(helpNodeMenu, "command", label = paste(" - Drag vertex: Move vertex"), 
                    command = function() {
                    })
                else if ((vertex.type == "OpenBlock")) {
                  tkadd(helpNodeMenu, "command", label = paste(" - Drag block head:    Move block"), 
                    command = function() {
                    })
                  tkadd(helpNodeMenu, "command", label = paste(" - Drag block corner:  Resize block"), 
                    command = function() {
                    })
                }
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(helpNodeMenu, "command", label = paste(" - Drag block:  Move minimized block"), 
                    command = function() {
                    })
                if ((vertex.type == "Vertex") || (vertex.type == 
                  "Factor")) 
                  tkadd(helpNodeMenu, "command", label = paste(" - Drag label:    Move vertex label"), 
                    command = function() {
                    })
                else if ((vertex.type == "OpenBlock")) {
                }
                else if ((vertex.type == "ClosedBlock")) 
                  tkadd(helpNodeMenu, "command", 
                        label = paste(" - Drag label:   Move label of minimized block"), 
                    command = function() {
                    })
                tkadd(nodePopupMenu, "cascade", label = "Help on node", 
                  menu = helpNodeMenu)
                methNodeMenu <- tkmenu(nodePopupMenu, tearoff = FALSE)
                if (hasMethod("addToPopups", class(vertex))) 
                  addToPopups(vertex, vertex.type, methNodeMenu, 
                    i, sinkView, Args)
                tkadd(nodePopupMenu, "cascade", label = "Items by method 'addToPopups'", 
                  menu = methNodeMenu)
                userNodeMenu <- tkmenu(nodePopupMenu, tearoff = FALSE)
                "NodePopup" <- function(UserNodePopupItems, item, 
                  vertex.type) {
                  force(UserNodePopupItems)
                  force(item)
                  force(vertex.type)
                  function(...) {
                    sinkView(UserNodePopupItems[[item]], blocks = TRUE)
                    UserNodePopupItems[[item]]$command(object, 
                      retVertexName(i, vertex.type), type = vertex.type, 
                      index = i, Arguments = Args())
                  }
                }
                if (length(UserNodePopupItems) > 0) 
                  for (item in seq(along = UserNodePopupItems)) 
                    if ((names(UserNodePopupItems[item]) == vertex.type)) 
                    tkadd(userNodeMenu, "command", label = UserNodePopupItems[[item]]$label, 
                      command = NodePopup(UserNodePopupItems, 
                        item, vertex.type))
                tkadd(nodePopupMenu, "cascade", label = "User defined items", 
                  menu = userNodeMenu)
            }
            "setBindNode" <- function(canvas, vertex, tag, result, 
                label, i, vertex.type, UserNodePopupItems) {
                "f" <- function(item, i, label = FALSE) {
                  tkitembind(canvas, item, "<Leave>", function() tkconfigure(canvas, 
                    cursor = "arrow"))
                  if (label) {
                    if (vertex.type == "ClosedBlock") 
                      tkitembind(canvas, item, "<Enter>", function() tkconfigure(canvas, 
                        cursor = "diamond_cross"))
                    else tkitembind(canvas, item, "<Enter>", 
                      function() tkconfigure(canvas, cursor = "hand2"))
                  }
                  else if (vertex.type == "ClosedBlock") 
                    tkitembind(canvas, item, "<Enter>", function() tkconfigure(canvas, 
                      cursor = "cross"))
                  else tkitembind(canvas, item, "<Enter>", function() tkconfigure(canvas, 
                    cursor = "crosshair"))
                  if (label) 
                    tkitembind(canvas, item, "<Button-1>", newEdge(i, 
                      vertex.type, slave = FALSE))
                  else tkitembind(canvas, item, "<Button-1>", 
                    newEdge(i, vertex.type, slave = FALSE))
                  # tkitembind(canvas, item, "<Button-2>", function(...) {
                  #   print("--2--")
                  # })
                  tkitembind(canvas, item, "<Up>", function(...) {
                    function(...) print("--UP%--")
                  })
                  tkitembind(canvas, item, "<Down>", function(...) {
                    print("--Down%--")
                  })
                  tkitembind(canvas, item, "<Left>", function(...) {
                    print("--Left%--")
                  })
                  tkitembind(canvas, item, "<Right>", function(...) {
                    print("--Right%--")
                  })
                  tkitembind(canvas, item, "<Home>", function(...) {
                    print("--PgUp%--")
                  })
                  tkitembind(canvas, item, "<End>", function(...) {
                    print("--PgDn%--")
                  })
                  tkitembind(canvas, item, "<Delete>", function(...) {
                    print("--Delete%--")
                  })
                  tkitembind(canvas, item, "<F1>", function(...) {
                    print("--F1%--")
                  })
                  tkitembind(canvas, item, "<Alt-1>", function(...) {
                    print("--A1%--")
                  })
                  tkitembind(canvas, item, "<Alt_L>", function(...) {
                    print("--A%--")
                  })
                  tkitembind(canvas, item, "<Option-1>", activateVertex(i, 
                    vertex.type, hit.type = "option-1", color = "DarkGreen"))
                  tkitembind(canvas, item, "<Shift-1>", activateVertex(i, 
                    vertex.type, hit.type = "shift-1", color = "DarkGreen"))
                  tkitembind(canvas, item, "<Control-1>", activateVertex(i, 
                    vertex.type, hit.type = "control-1", color = "SeaGreen"))
                  tkitembind(canvas, item, "<Shift-Control-1>", 
                    activateVertex(i, vertex.type, hit.type = "shift-control-1", 
                      color = "LightSeaGreen"))
                  tkitembind(canvas, item, "<Option-3>", activateVertex(i, 
                    vertex.type, hit.type = "option-3", color = "LightGreen"))
                  tkitembind(canvas, item, "<Shift-3>", activateVertex(i, 
                    vertex.type, hit.type = "shift-3", color = "LightGreen"))
                  tkitembind(canvas, item, "<Control-3>", activateVertex(i, 
                    vertex.type, hit.type = "control-3", color = "SpringGreen"))
                  tkitembind(canvas, item, "<Shift-Control-3>", 
                    activateVertex(i, vertex.type, hit.type = "shift-control-3", 
                      color = "LimeGreen"))
                  if (label) 
                    tkitembind(canvas, item, "<B1-Motion>", moveVertexLabel(i, 
                      vertex.type))
                  else tkitembind(canvas, item, "<B1-Motion>", 
                    moveVertex(i, vertex.type))
                  if (label) 
                    tkitembind(canvas, item, "<Double-Button-1>", 
                      changeVertexLabel(i, vertex.type))
                  else if (vertex.type == "ClosedBlock") 
                    tkitembind(canvas, item, "<Double-Button-1>", 
                      openBlock(i))
                  else tkitembind(canvas, item, "<Double-Button-1>", 
                    undisplayVertex(i, vertex.type, slave = FALSE))
                  if (label) 
                    tkitembind(canvas, item, "<Triple-Button-1>", 
                      deleteVertexLabel(i, vertex.type))
                  if (initial.set.popups) 
                    tkitembind(canvas, item, "<Button-3>", callPopup(i, 
                      nodePopupMenu))
                  else tkitembind(canvas, item, "<Button-3>", 
                    callPopupNode(vertex, i, vertex.type, UserNodePopupItems))
                  tkaddtag(canvas, tag, "withtag", item)
                }
                "blockitembind" <- function(item, label = FALSE) {
                  tkitembind(canvas, item, "<Leave>", function() tkconfigure(canvas, 
                    cursor = "arrow"))
                  if (label) 
                    tkitembind(canvas, item, "<Enter>", function() tkconfigure(canvas, 
                      cursor = "left_ptr"))
                  else tkitembind(canvas, item, "<Enter>", function() tkconfigure(canvas, 
                    cursor = "right_ptr"))
                  tkitembind(canvas, item, "<B1-Motion>", moveBlock(i, 
                    1))
                  tkitembind(canvas, item, "<Double-Button-1>", 
                    closeBlock(i))
                  if (initial.set.popups) 
                    tkitembind(canvas, item, "<Button-3>", callPopup(i, 
                      nodePopupMenu))
                  else tkitembind(canvas, item, "<Button-3>", 
                    callPopupNode(vertex, i, vertex.type, UserNodePopupItems))
                  tkaddtag(canvas, tag, "withtag", item)
                }
                if (initial.set.popups) {
                  nodePopupMenu <- tkmenu(canvas, tearoff = FALSE)
                  addNodePopups(vertex, i, vertex.type, nodePopupMenu, 
                    UserNodePopupItems)
                }
                if (vertex.type == "OpenBlock") {
                  blockitembind(label, label = TRUE)
                  if (!is.null(result)) 
                    if (length(result) > 0) 
                      for (k in seq(length(result))) blockitembind(result[[k]])
                }
                else {
                  f(label, i, TRUE)
                  if (!is.null(result$dynamic)) 
                    if (length(result$dynamic) > 0) 
                      for (k in seq(length(result$dynamic))) f(result$dynamic[[k]], 
                        i, FALSE)
                  if (!is.null(result$fixed)) 
                    if (length(result$fixed) > 0) 
                      for (k in seq(length(result$fixed))) f(result$fixed[[k]], 
                        i, FALSE)
                }
            }
            "subDrawVertex" <- function(vertex, i, w = control$w, 
                vertexcolor = vertexcolor, vertex.type = ifelse(i > 
                  0, "Vertex", "Factor"), setTag = TRUE) {
                tag <- getTag(vertex.type, i, setTag = setTag)
                pos <- retVertexPos(i, vertex.type)
                if (hasMethod("draw", class(vertex))) 
                  dot <- draw(vertex, canvas, pos, x = pos[1], 
                    y = pos[2], stratum = retStratum(i, vertex.type = vertex.type), 
                    w = w * Scale, color = vertexcolor, background = control$background)
                else {
                  s <- w * sqrt(4/pi) * Scale
                  p <- tkcreate(canvas, "oval", pos[1] - s, pos[2] - 
                    s, pos[1] + s, pos[2] + s, fill = vertexcolor, 
                    activefill = "OrangeRed")
                  dot <- list(dynamic = list(p), fixed = NULL)
                }
                label <- tkcreate(canvas, "text", pos[1] + w, 
                  pos[2], text = retVertexLabel(i, vertex.type), 
                  anchor = "nw", font = font.vertex.label, activefill = "DarkSlateGray")
                if (control$debug.strata && (vertex.type != "Factor") && 
                  (vertex.type != "Extra")) {
                  strata <- retStratum(i, vertex.type)
                  block <- retBlockIndex(i, vertex.type)
                  color <- myColor(strata)
                  numbers <- tkcreate(canvas, "text", pos[1] - 
                    4 * w, pos[2] - 4 * w, text = paste(i, strata, 
                    block, sep = "."), fill = color, anchor = "nw", 
                    font = "12x30", activefill = "DarkSlateGray")
                  tkaddtag(canvas, tag, "withtag", numbers)
                }
                else numbers <- NULL
                if (vertex.type != "OpenBlock") 
                  setBindNode(canvas, vertex, tag, dot, label, 
                    i, vertex.type, control$UserMenus)
                else setBindNode(canvas, vertex, tag, dot, label, 
                  i, vertex.type, control$UserMenus)
                return(list(tag = tag, dot = dot, label = label, 
                  numbers = numbers))
            }
            "drawVertex" <- function(i, w = control$w, vertexcolor = vertexcolor, 
                vertex.type = ifelse(i > 0, "Vertex", "Factor"), 
                setTag = TRUE) {
                if (vertex.type == "ClosedBlock") 
                  itemsClosedBlocks[[i]] <<- subDrawVertex(blockList[[i]], 
                    i, w = w, vertexcolor = vertexcolor, vertex.type = vertex.type, 
                    setTag = setTag)
                else if (vertex.type == "Vertex") 
                  itemsVertices[[i]] <<- subDrawVertex(vertexList[[i]], 
                    i, w = w, vertexcolor = vertexcolor, vertex.type = vertex.type, 
                    setTag = setTag)
                else if (vertex.type == "Factor") 
                  itemsFactors[[-i]] <<- subDrawVertex(dg@factorVertexList[[-i]], 
                    i, w = w, vertexcolor = vertexcolor, vertex.type = vertex.type, 
                    setTag = setTag)
                else if (vertex.type == "Extra") 
                  itemsExtras[[i]] <<- subDrawVertex(dg@extraList[[i]], 
                    i, w = w, vertexcolor = vertexcolor, vertex.type = vertex.type, 
                    setTag = setTag)
            }
            "clearSelectedVertices" <- function() {
                if (length(selectedNodes) > 0) {
                  lapply(selectedNodes, function(k) setVertexColor(k$index, 
                    color = "red", vertex.type = k$node.type))
                  selectedNodes <<- list()
                }
            }
            "setActivatedVertex" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor"), hit.type = "none") {
                if (hit.type == "none") 
                  activatedNode <<- list(number = i, vertex.type = vertex.type)
                if (!((i == 0) && (hit.type == "none"))) {
                  a <- list(index = i, node.type = vertex.type, 
                    hit.type = hit.type)
                  if (!any(unlist(lapply(selectedNodes, function(i) all(unlist(i) == 
                    unlist(a)))))) 
                    selectedNodes <<- append(list(a), selectedNodes)
                }
            }
            "retActivatedVertex" <- function() return(activatedNode[[1]])
            "retActivatedVertexVertex.Type" <- function() return(activatedNode[[2]])
            "deActivateVertex" <- function(i, color = retVertexColor(i, 
                vertex.type), vertex.type = ifelse(i > 0, "Vertex", 
                "Factor"), new.edge = FALSE) {
                if (length(selectedNodes) > 0) {
                  same <- (retActivatedVertex() == i)
                  x <- lapply(selectedNodes, function(k) if ((k$index == 
                    i) && (k$node.type == vertex.type) && ((new.edge && 
                    (k$hit.type == "none")) || (!new.edge) || 
                    same)) {
                    setVertexColor(i, color = color, vertex.type = vertex.type)
                    return(TRUE)
                  }
                  else return(FALSE))
                  selectedNodes <<- selectedNodes[!unlist(x)]
                }
                if ((retActivatedVertex() == i) && (retActivatedVertexVertex.Type() == 
                  vertex.type)) {
                  setActivatedVertex(0, "Null")
                  setVertexColor(i, color = color, vertex.type = vertex.type)
                  return(TRUE)
                }
                else return(FALSE)
            }
            "subActivateVertex" <- function(i, color = "green", 
                vertex.type = ifelse(i > 0, "Vertex", "Factor"), 
                hit.type = "none", new.edge = FALSE) {
                if (!(hit.type == "none")) {
                  setActivatedVertex(i, vertex.type, hit.type = hit.type)
                  setVertexColor(i, color = color, vertex.type = vertex.type)
                  return(TRUE)
                }
                else if (!deActivateVertex(i, "cyan", vertex.type, 
                  new.edge = new.edge)) 
                  if ((retActivatedVertex() == 0) && (vertex.type != 
                    "Extra")) {
                    setActivatedVertex(i, vertex.type)
                    setVertexColor(i, color = color, vertex.type = vertex.type)
                    return(TRUE)
                  }
                  else return(FALSE)
                else return(TRUE)
            }
            "activateVertex" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor"), color = "green", hit.type = "none") {
                force(i)
                force(vertex.type)
                force(hit.type)
                force(color)
                function(...) subActivateVertex(i, color = color, 
                  vertex.type = vertex.type, hit.type = hit.type)
            }
            "clearSelectedEdges" <- function() {
                if (length(selectedEdges) > 0) {
                  lapply(selectedEdges, function(k) setEdgeColor(k$index, 
                    color = NULL, edge.type = k$edge.type))
                  selectedEdges <<- list()
                }
            }
            "setActivatedEdge" <- function(i, from = -1, to = -1, 
                edge.type = ifelse(i > 0, "Edge", "Factor"), 
                hit.type = "none") {
                if (hit.type == "none") 
                  activatedEdge <<- list(number = i, edge.type = edge.type)
                if (!((i == 0) && (hit.type == "none"))) {
                  a <- list(index = i, from = from, to = to, 
                    edge.type = edge.type, hit.type = hit.type)
                  if (!any(unlist(lapply(selectedEdges, function(i) all(unlist(i) == 
                    unlist(a)))))) 
                    selectedEdges <<- append(list(a), selectedEdges)
                }
            }
            "subActivateEdge" <- function(i, from = -1, to = -1, 
                edge.type = "VertexEdge", color = "green", hit.type = "none") {
                if ((hit.type != "none")) {
                  setActivatedEdge(i, from = from, to = to, edge.type, 
                    hit.type = hit.type)
                  if (i > 0) 
                    setEdgeColor(i, edge.type = edge.type, color)
                }
                else {
                  if (activatedEdge$number > 0) {
                    setEdgeColor(activatedEdge$number, edge.type = activatedEdge$edge.type, 
                      color = NULL)
                    if ((i != 0) && (length(selectedEdges) > 
                      0)) {
                      same <- (retActivatedEdge() == i)
                      new.edge <- TRUE
                      x <- lapply(selectedEdges, function(k) if ((k$hit.type == 
                        "none")) {
                        if (k$index != 0) 
                          setEdgeColor(k$index, color = NULL, 
                            edge.type = k$edge.type)
                        return(TRUE)
                      }
                      else return(FALSE))
                      assign("x", x, pos = 1)
                      selectedEdges <<- selectedEdges[!unlist(x)]
                    }
                  }
                  setActivatedEdge(i, from = from, to = to, edge.type, 
                    hit.type = hit.type)
                  if (i > 0) 
                    setEdgeColor(i, edge.type = edge.type, "Green")
                }
            }
            "activateEdge" <- function(i, from = -1, to = -1, edge.type = "VertexEdge", 
                color = "green", hit.type = "none") {
                force(i)
                force(from)
                force(to)
                force(edge.type)
                force(hit.type)
                force(color)
                function(...) subActivateEdge(i, from = from, 
                  to = to, color = color, edge.type = edge.type, 
                  hit.type = hit.type)
            }
            "retActivatedEdge" <- function(edge.type = "VertexEdge") return(activatedEdge$number)
            "updateSelectedEdges" <- function(R, vertexEdges, edge.type = "VertexEdge", 
                setVertex = FALSE, updateBE = FALSE) {
                if (updateBE) {
                  updateBlockEdges()
                  sinkBlockEdges()
                }
                if (is.element("dg", names(R))) 
                  ldg <- R$dg
                else ldg <- .newDgGraphEdges(vertexList = vertexList, 
                  visibleVertices = R$VisibleVertices, visibleBlocks = R$VisibleBlocks, 
                  edgeList = vertexEdges, blockList = blockList, 
                  factorVertexList = R$FactorVertices, factorEdgeList = R$FactorEdges, 
                  extraList = R$ExtraVertices, extraEdgeList = R$ExtraEdges)
                setModel(R$object, dg = ldg, txt = "updateSelectedEdges", 
                  RR = R)
                activateEdge(0, edge.type = NULL)()
                if (setVertex) 
                  setActivatedVertex(0, "Vertex")
                clearSelectedVertices()
                clearSelectedEdges()
                clearFactorEdges()
                clearExtraEdges()
                update.edge.labels()
                if (!is.null(ldg@factorVertexList) && !is.null(ldg@factorEdgeList)) 
                  drawFactors(ldg@factorEdgeList, ldg@factorVertexList)
            }
            "drawResult" <- function(newEdges, R, slave, txt, Arguments = Args()) {
                Edges <- extractEdgesResult(R, newEdges, TRUE, 
                  txt)
                if (is.null(Edges)) 
                  Edges <- new("dg.VertexEdgeList")
                if (is.element("dg", names(R))) {
                  ldg <- R$dg
                }
                else {
                  ldg <- .newDgGraphEdges(vertexList = vertexList, 
                    visibleVertices = R$VisibleVertices, visibleBlocks = R$VisibleBlocks, 
                    edgeList = Edges, blockList = blockList, 
                    blockEdgeList = R$BlockEdges, factorVertexList = R$FactorVertices, 
                    factorEdgeList = R$FactorEdges, extraList = R$ExtraVertices, 
                    extraEdgeList = R$ExtraEdges)
                }
                if (any(slotNames(R$object) == ".title"))
                  newtitle <- R$object@.title
                else
                  newtitle <- dgm.frameModels@label
                if (slave) {
                  drawModel(frameModels = dgm.frameModels,
                  # frameViews = dm.frameViews, 
                    graphWindow = NULL, dg = ldg, object = R$object, 
                    title = newtitle, control = control, Arguments = Arguments)
                }
                else {
                  tktitle(GW.top) <- newtitle
                  setModel(R$object, dg = ldg, txt = txt, RR = R)
                  redrawView(frameModels = dgm.frameModels, frameViews = dm.frameViews, 
                    graphWindow = GraphWindow, dg = ldg,
                    control = control, Arguments = Arguments)
                }
            }
            "updateNode" <- function(i, vertex.type = "Vertex", 
                items = itemsVertices, ii = i, pos = retVertexPos(ii, 
                  vertex.type), redrawVertices = TRUE) {
                update <- FALSE
                update.edges <- FALSE
                if (redrawVertices) {
                  tkdelete(canvas, vertexItem(ii, vertex.type = vertex.type)$tag)
                  drawVertex(ii, w = control$w, vertexcolor = control$vertexColor, 
                    vertex.type = vertex.type, setTag = FALSE)
                  vertexColor <- retVertexColor(ii, vertex.type = vertex.type)
                  setVertexColor(ii, color = vertexColor, vertex.type = vertex.type)
                }
                else {
                  xy <- tkcoords(canvas, items[[i]]$tag)
                  xy <- apply(matrix(as.numeric(xy), ncol = 2, 
                    byrow = 2), 2, mean)
                  if (!(all(is.na(xy)))) {
                    update <- TRUE
                    if (!any(is.nan(xy)) && (length(xy) > 0)) {
                      dxy <- findDifference(pos, c(xy, rep(0, 
                        local.N - 2)))
                      ll <- sum(dxy[1:2]^2)
                      if ((vertex.type == "Vertex") || (vertex.type == 
                        "ClosedBlock")) 
                        if (is.numeric(ll)) {
                          if ((length(ll) > 0) && (ll > 0)) 
                            update.edges <- TRUE
                        }
                        else message(paste("Invalid length: ", 
                          ll))
                      tkmove(canvas, items[[i]]$tag, dxy[1], 
                        dxy[2])
                    }
                    if (vertex.type == "ClosedBlock") {
                      posLabel <- pos + retBlockLabelPos(ii)
                      tkcoords(canvas, items[[i]]$label, posLabel[1], 
                        posLabel[2])
                      tkitemconfigure(canvas, items[[i]]$label, 
                        text = retVertexLabel(i, vertex.type = vertex.type))
                    }
                    else {
                      posLabel <- retLabelPos(ii, vertex.type = vertex.type)
                      xyl <- as.numeric(tkcoords(canvas, items[[i]]$l))
                      tkitemconfigure(canvas, items[[i]]$l, text = retVertexLabel(ii, 
                        vertex.type = vertex.type))
                      vertexColor <- retVertexColor(ii, vertex.type = vertex.type)
                      setVertexColor(ii, color = vertexColor, 
                        vertex.type = vertex.type)
                      if (!any(is.nan(xyl)) && (length(xyl) > 
                        0)) {
                        dxy <- findDifference(posLabel, c(xyl, 
                          rep(0, local.N - 2)))
                        tkmove(canvas, items[[i]]$l, dxy[1], 
                          dxy[2])
                      }
                    }
                  }
                }
                return(list(update = update, update.edges = update.edges))
            }
            "subUpdatePanel" <- function() {
                if (control$variableFrame) {
                  if ((get("type", GW.top$env$box$env) == "variableList")) {
                  }
                  else {
                    tkinsert.blockList(GW.top$env$box, blockList, 
                      delete = TRUE)
                    tkinsert.blockList(GW.top$env$box, blockList, 
                      delete = FALSE)
                  }
                }
            }
            "subUpdateGraphWindow" <- function(txt = "", redrawVertices = FALSE, 
                raiseEdges = FALSE, updateEdges = FALSE, all.blockframes = FALSE, 
                blockframes = NULL) {
                if (control$debug.update) 
                  print(paste("subUpdateGraphWindow:", txt, " (0)", 
                    getLabel()))
                if (control$variableFrame) 
                  subUpdatePanel()
                pos <- NULL
                update.edges <- TRUE
                if (!.IsEmpty(blockList)) 
                  for (i in seq(along = blockList)) {
                    if ((closedBlock[i] || hiddenBlock[i]) && 
                      (i %in% dg@visibleBlocks)) {
                      if (all(positionsClosedBlocks[i, ] < rep(-100, 
                        local.N))) {
                        deleteBlock(i)
                        visible.Blocks <- returnVisibleBlocks()
                        visible.Blocks <- visible.Blocks[visible.Blocks != 
                          i]
                        setVisibleBlocks(visible.Blocks)
                        openBlock(i, update = FALSE)()
                      }
                      else {
                        pos <- retVertexPos(i, "ClosedBlock")
                        update.edges <- updateEdges || is.element(i, 
                          blockframes)
                        if (!hiddenBlock[i] && !is.null(itemsClosedBlocks[[i]])) {
                          R <- updateNode(i, vertex.type = "ClosedBlock", 
                            items = itemsClosedBlocks)
                          update.edges <- update.edges | R$update.edges
                        }
                      }
                      if (!is.null(pos) && !is.null(dg@blockEdgeList)) 
                        if ((update.edges || raiseEdges) && (length(itemsBlockEdges[[i]]) > 
                          0)) 
                          for (e in itemsBlockEdges[[i]]) if (!(is.null(e))) 
                            if (TRUE) 
                              setEdgeCoords(e, edge.type = "BlockEdge", 
                                f = i, posFrom = pos, from.type = vertexTypeOfEdge(-i, 
                                  e$type), raise = FALSE, setEdgeLabel = FALSE)
                    }
                    else if (all.blockframes || (is.element(i, 
                      blockframes))) {
                      if (all(positionsClosedBlocks[i, ] < rep(-100, 
                        local.N))) {
                        deleteBlock(i)
                      }
                      else {
                        tkitemconfigure(canvas, itemsOpenBlocks[[i]]$label, 
                          text = retVertexLabel(i, vertex.type = "ClosedBlock"))
                        tkcoordsBlock(i, lower = FALSE)
                      }
                    }
                  }
                if (!is.null(itemsExtras)) 
                  for (i in seq(along = itemsExtras)) if (!is.null(itemsExtras[[i]]) && 
                    !is.null(itemsExtras[[i]][[1]])) {
                    updateNode(i, vertex.type = "Extra", items = itemsExtras, 
                      redrawVertices = redrawVertices)
                    if (!is.null(pos) && (update.edges || raiseEdges) && 
                      (length(itemsExtraEdges[[i]]) > 0)) 
                      for (e in itemsExtraEdges[[i]]) if (!(is.null(e))) 
                        if (TRUE) 
                          setEdgeCoords(e, edge.type = "ExtraEdge", 
                            f = -i, posFrom = pos)
                  }
                if (!is.null(itemsFactors)) 
                  for (i in seq(along = itemsFactors)) if (!is.null(itemsFactors[[i]]) && 
                    !is.null(itemsFactors[[i]][[1]])) {
                    vertex.indices <- dg@factorVertexList[[i]]@vertex.indices
                    setFactorVertexPosition(i, vertex.indices)
                    updateNode(i, vertex.type = "Factor", items = itemsFactors, 
                      ii = -i, redrawVertices = redrawVertices)
                  }
                for (i in seq(along = itemsVertices)) if (!is.null(itemsVertices[[i]]) && 
                  !is.null(itemsVertices[[i]][[1]])) {
                  pos <- retVertexPos(i, "Vertex")
                  update.edges <- updateEdges || redrawVertices
                  if (!closedVertex[i]) {
                    R <- updateNode(i, pos = pos, redrawVertices = redrawVertices)
                    update.edges <- update.edges || R$update.edges
                    if (R$update && !redrawVertices) {
                      if (control$debug.strata) {
                        strata <- retStratum(i, vertex.type = "Vertex")
                        block <- retBlockIndex(i, vertex.type = "Vertex")
                        color <- myColor(strata)
                        tkitemconfigure(canvas, itemsVertices[[i]]$numbers, 
                          text = paste(i, strata, block, sep = "."))
                        tkitemconfigure(canvas, itemsVertices[[i]]$numbers, 
                          fill = color)
                      }
                    }
                  }
                  else update.edges <- TRUE
                  if (!is.null(pos) && (update.edges || raiseEdges) && 
                    (length(itemsEdges[[i]]) > 0)) 
                    for (e in itemsEdges[[i]]) if (!(is.null(e))) 
                      if (TRUE) 
                        setEdgeCoords(e, edge.type = "vertexEdge", 
                          f = i, posFrom = pos)
                }
            }
            "setUpdateVertices" <- function(txt = "") {
                updateCountVerticesMain <<- updateCountVerticesMain + 
                  1
                updateCountVertices <<- updateCountVerticesMain
            }
            "setUpdatePositions" <- function(txt = "") {
                updateCountPositionsMain <<- updateCountPositionsMain + 
                  1
                updateCountPositions <<- updateCountPositionsMain
            }
            "setUpdateBlocks" <- function(txt = "") {
                updateCountBlocksMain <<- updateCountBlocksMain + 
                  1
                updateCountBlocks <<- updateCountBlocksMain
            }
            "setUpdateBlockEdges" <- function(txt = "", local = TRUE) {
                updateCountBlockEdgesMain <<- updateCountBlockEdgesMain + 
                  1
                if (!local) 
                  updateCountBlockEdges <<- updateCountBlockEdgesMain
            }
            "setUpdateAll" <- function(txt = "") {
                updateCountVerticesMain <<- updateCountVerticesMain + 
                  1
                updateCountPositionsMain <<- updateCountPositionsMain + 
                  1
                updateCountBlocksMain <<- updateCountBlocksMain + 
                  1
                updateCountBlockEdgesMain <<- updateCountBlockEdgesMain + 
                  1
                updateCountVertices <<- updateCountVerticesMain
                updateCountPositions <<- updateCountPositionsMain
                updateCountBlocks <<- updateCountBlocksMain
                updateCountBlockEdges <<- updateCountBlockEdgesMain
            }
            "subUpdatePositions" <- function(txt = "") {
                if (updateWindow) {
                  testUpdateModel()
                  n <- length(vertexList)
                  m <- length(itemsVertices)
                  if (n > m) 
                    for (i in seq(n - m)) {
                      closedVertex <<- c(closedVertex, FALSE)
                      itemsVertices <<- append(itemsVertices, 
                        list(NULL))
                      itemsEdges <<- append(itemsEdges, list(NULL))
                    }
                  n <- length(blockList)
                  m <- length(itemsClosedBlocks)
                  if (n > m) {
                    for (i in seq(n - m)) {
                      itemsBlockEdges <<- append(itemsBlockEdges, 
                        list(NULL))
                      itemsClosedBlocks <<- append(itemsClosedBlocks, 
                        list(NULL))
                      itemsOpenBlocks <<- append(itemsOpenBlocks, 
                        list(NULL))
                      closedBlock <<- c(closedBlock, FALSE)
                      hiddenBlock <<- c(hiddenBlock, FALSE)
                      openTreeBlock <<- c(openTreeBlock, FALSE)
                      dg@visibleBlocks <- unique(sort(c(dg@visibleBlocks, 
                        m + i)))
                      dg@visibleBlocks <<- dg@visibleBlocks[dg@visibleBlocks != 
                        0]
                    }
                    Arguments <- Args()
                    redrawView(frameModels = dgm.frameModels, 
                      frameViews = dm.frameViews, graphWindow = GraphWindow, 
                      dg = dg, control = control, Arguments = Arguments)
                  }
                  all.blockframes <- updateCountBlocks < updateCountBlocksMain
                  updateEdges <- updateCountBlockEdges < updateCountBlockEdgesMain
                  updateVertices <- updateCountVertices < updateCountVerticesMain
                  updatePositions <- updateCountPositions < updateCountPositionsMain
                  if (updateCountBlockEdges < updateCountBlockEdgesMain) 
                    updateAllBlockIndices()
                  if (updatePositions || updateVertices || updateEdges || 
                    all.blockframes) 
                    subUpdateGraphWindow(txt, updateEdges = updateEdges, 
                      redrawVertices = updateVertices, all.blockframes = all.blockframes)
                  if (updateEdges) 
                    updateBlockEdges()
                  updateCountVertices <<- updateCountVerticesMain
                  updateCountPositions <<- updateCountPositionsMain
                  updateCountBlocks <<- updateCountBlocksMain
                  updateCountBlockEdges <<- updateCountBlockEdgesMain
                }
            }
            "updatePositions" <- function(txt = "") {
                force(txt)
                function(...) {
                  subUpdatePositions(txt)
                }
            }
            "tkdeleteRectangleCorner" <- function(line) {
                tkdelete(canvas, line[[1]])
                tkdelete(canvas, line[[2]])
                if (!is.null(transformation)) 
                  tkdelete(canvas, line[[3]])
            }
            "tkdeleteRectangle" <- function(rectangle) {
                line <- rectangle$Lines
                tkdelete(canvas, line[[1]])
                tkdelete(canvas, line[[2]])
                tkdelete(canvas, line[[3]])
                tkdelete(canvas, line[[4]])
                if (!is.null(transformation)) {
                  tkdelete(canvas, line[[5]])
                  tkdelete(canvas, line[[6]])
                  tkdelete(canvas, line[[7]])
                  tkdelete(canvas, line[[8]])
                  tkdelete(canvas, line[[9]])
                  tkdelete(canvas, line[[10]])
                  tkdelete(canvas, line[[11]])
                  tkdelete(canvas, line[[12]])
                }
                corner <- rectangle$Corners
                tkdeleteRectangleCorner(corner[[1]])
                tkdeleteRectangleCorner(corner[[2]])
                tkdeleteRectangleCorner(corner[[3]])
                tkdeleteRectangleCorner(corner[[4]])
                if (!is.null(transformation)) {
                  tkdeleteRectangleCorner(corner[[5]])
                  tkdeleteRectangleCorner(corner[[6]])
                  tkdeleteRectangleCorner(corner[[7]])
                  tkdeleteRectangleCorner(corner[[8]])
                }
            }
            "tkdeleteBar" <- function(line) {
                tkdelete(canvas, line[[1]])
                tkdelete(canvas, line[[2]])
                tkdelete(canvas, line[[3]])
                tkdelete(canvas, line[[4]])
            }
            "deleteBlock" <- function(i) {
                tkdeleteRectangle(openBlockItem(i)$rectangle)
                tkdeleteBar(openBlockItem(i)$bar)
                tkdelete(canvas, openBlockItem(i)$label)
                if (!is.null(openBlockItem(i)$canvas)) 
                  tkdelete(canvas, openBlockItem(i)$canvas)
            }
            "subSubDeleteEdge" <- function(i, f, t, edge.type = "VertexEdge") {
                E <- getEdges(edge.type = edge.type)[[i]]
                "delete" <- function(k, edges) {
                  edges <- edgeItem(k, edge.type = edge.type)
                  if (length(edges) > 0) 
                    for (e in edges) if (!(is.null(e))) 
                      if ((e$nr == i) && (e$type == edge.type)) 
                        if (e$to > k) {
                          if (control$debug.edges) {
                            cat("Tkdelete, subSubDeleteEdge: ")
                            for (l in 1:length(e$edges)) cat(paste(e$edges[[l]], 
                              " "))
                            for (l in 1:length(e$tags)) cat(paste(e$tags[[l]], 
                              " "))
                            cat(paste(e$label))
                            cat("\n")
                          }
                          for (l in 1:length(e$edges)) tkdelete(canvas, 
                            e$edges[[l]])
                          for (l in 1:length(e$tags)) tkdelete(canvas, 
                            e$tags[[l]])
                          tkdelete(canvas, e$label)
                        }
                }
                for (j in E@vertex.indices) delete(j, edgeItem(j, 
                  edge.type = edge.type))
                "remove" <- function(k, edges) if (length(edges) > 
                  0) {
                  result <- NULL
                  for (e in edges) if (!(is.null(e))) 
                    if (!((e$nr == i) && (e$type == edge.type))) 
                      result <- c(result, list(e))
                  if (is.null(result)) 
                    setEdgeItem(k, edge.type = edge.type, list(NULL))
                  else setEdgeItem(k, edge.type = edge.type, 
                    result)
                }
                for (j in E@vertex.indices) remove(j, edgeItem(j, 
                  edge.type = edge.type))
            }
            "subSubUndisplayFactorVertex" <- function(i, edge.type = "FactorEdge") {
                edges <- edgeItem(i, edge.type = edge.type)
                if (length(edges) > 0) 
                  for (e in edges) if (!(is.null(e))) {
                    subSubDeleteEdge(e$nr, i, e$to, edge.type = edge.type)
                    clearEdge(e$nr, edge.type = edge.type)
                  }
                visible.Vertices <- returnVisibleVertices()
                visible.Vertices <- visible.Vertices[visible.Vertices != 
                  i]
                setVisibleVertices(visible.Vertices)
                tkdelete(canvas, vertexItem(i)$tag)
                setVertexItem(i, list(NULL))
            }
            "clearFactorEdges" <- function() {
                for (f in seq(along = itemsFactors)) subSubUndisplayFactorVertex(-f)
            }
            "clearExtraEdges" <- function() {
            }
            "setVisibleVertices" <- function(i) {
                dg@visibleVertices <<- i
                GraphWindow@dg@visibleVertices <<- dg@visibleVertices
            }
            "returnVisibleVertices" <- function() {
                return(dg@visibleVertices)
            }
            "setVisibleBlocks" <- function(i) {
                dg@visibleBlocks <<- i
                GraphWindow@dg@visibleBlocks <<- dg@visibleBlocks
            }
            "returnVisibleBlocks" <- function() {
                return(dg@visibleBlocks)
            }
            "subSubUndisplayVertex" <- function(i, edge.type = "VertexEdge") {
                if (control$debug.position) 
                  print(paste("subSubUndisplayVertex", i, edge.type))
                edges <- edgeItem(i, edge.type = edge.type)
                if (length(edges) > 0) 
                  for (e in edges) if (!(is.null(e))) 
                    if ((e$type == edge.type)) {
                      subSubDeleteEdge(e$nr, i, e$to, edge.type = edge.type)
                      clearEdge(e$nr, edge.type = edge.type)
                    }
                visible.Vertices <- returnVisibleVertices()
                visible.Vertices <- visible.Vertices[visible.Vertices != 
                  i]
                setVisibleVertices(visible.Vertices)
                tkdelete(canvas, vertexItem(i)$tag)
                setVertexItem(i, list(NULL))
                if (control$variableFrame) {
                }
            }
            "update.edge.labels" <- function() {
                "subUpdateEdgeLabels" <- function(itemsNodes, edge.type = "VertexEdge") 
                  for (f in seq(along = itemsNodes)) 
                    if (!is.null(itemsNodes[[f]]) && !is.null(itemsNodes[[f]][[1]])) {
                  if (edge.type != "VertexEdge") 
                    f <- -f
                  edges <- edgeItem(f, edge.type = edge.type)
                  if (length(edges) > 0) 
                    for (e in edges) if (!(is.null(e))) 
                      if (e$to < f) {
                        display <- displayEdge(f, e$to, edgeNode = e)
                        if (!is.na(control$namesOnEdges) && control$namesOnEdges && 
                          display) {
                          vertexnames <- c(retVertexName(f, vertexTypeOfEdge(f, 
                            e$type)), retVertexName(e$to, vertexTypeOfEdge(e$to, 
                            e$type)))
                          if (e$reverse) 
                            vertexnames <- rev(vertexnames)
                          label <- paste(vertexnames, collapse = "~")
                        }
                        else label <- ""
                        setEdgeLabel(e, label, e$label.number, 
                          f = f, permanent = TRUE)
                        setEdgeWidth(e, 2, e$label.number, f = f)
                      }
                }
                if (control$updateEdgeLabels) {
                  subUpdateEdgeLabels(itemsVertices, edge.type = "VertexEdge")
                  subUpdateEdgeLabels(itemsFactors, edge.type = "FactorEdge")
                  subUpdateEdgeLabels(itemsExtras, edge.type = "ExtraEdge")
                  subUpdateEdgeLabels(itemsClosedBlocks, edge.type = "BlockEdge")
                }
                if (!is.na(control$namesOnEdges) && !control$namesOnEdges) {
                  for (i in seq(along = GraphWindow@dg@edgeList)) 
                    label(GraphWindow@dg@edgeList[[i]]) <<- ""
                  for (i in seq(along = GraphWindow@dg@factorEdgeList)) 
                    label(GraphWindow@dg@factorEdgeList[[i]]) <<- ""
                  for (i in seq(along = GraphWindow@dg@extraEdgeList)) 
                    label(GraphWindow@dg@extraEdgeList[[i]]) <<- ""
                  for (i in seq(along = GraphWindow@dg@blockEdgeList)) 
                    label(GraphWindow@dg@blockEdgeList[[i]]) <<- ""
                }
            }
            "updateBlockEdges" <- function() {
                Edges <- selectCurrentEdges(omitEdges = FALSE, 
                  edge.type = "VertexEdge")
                edge.list <- lapply(Edges, function(i) i@vertex.indices)
                sinkVertexList()
                NewBlockEdges <- returnBlockEdgeList(edge.list, 
                  vertexList, blockList, color = control$blockEdgeColor, 
                  visibleBlocks = dg@visibleBlocks, oriented = dg@oriented)
                new.list <- lapply(NewBlockEdges, function(i) {
                  x <- i@vertex.indices
                  names(x) <- NULL
                  x
                })
                old.list <- lapply(getEdges(edge.type = "BlockEdge"), 
                  function(i) {
                    x <- i@vertex.indices
                    names(x) <- NULL
                    x
                  })
                match.old.new <- match(old.list, new.list)
                for (i in seq(along = match.old.new)) if (is.na(match.old.new[i])) 
                  if (all(abs(old.list[[i]]) > 0)) {
                    subSubDeleteEdge(i, old.list[[i]][1], old.list[[i]][2], 
                      edge.type = "BlockEdge")
                    clearEdge(i, edge.type = "BlockEdge")
                  }
                match.new.old <- match(new.list, old.list)
                for (i in seq(along = match.new.old)) if (is.na(match.new.old[i])) {
                  E <- append.edge(NewBlockEdges[[i]], edge.type = "BlockEdge")
                  drawEdge(E[[length(E)]], length(E), lower = TRUE, 
                    edge.type = "BlockEdge")
                }
            }
            "deleteAllEdgeLabels" <- function(permanent = TRUE) {
                function(...) {
                  "g" <- function(items, edge.type = "VertexEdge") 
                    for (f in seq(along = items)) 
                      if (TRUE || !is.null(items[[f]]) && !is.null(items[[f]][[1]])) {
                    ff <- ifelse(edge.type == "VertexEdge", f, 
                      -f)
                    if (control$debug.edges) 
                      cat(paste(edge.type, f, ff, ": "))
                    edges <- edgeItem(ff, edge.type = edge.type)
                    if (length(edges) > 0) 
                      for (e in edges) if (!(is.null(e))) {
                        if (control$debug.edges) 
                          cat(paste(" ", e$to))
                        if ((e$to < f) || ((ff < 0))) {
                          setEdgeLabel(e, "", e$label.number, 
                            f = ff, permanent = permanent)
                          setEdgeWidth(e, 2, e$label.number, 
                            f = ff)
                        }
                      }
                    if (control$debug.edges) 
                      cat("\n")
                  }
                  g(itemsVertices, edge.type = "VertexEdge")
                  g(itemsClosedBlocks, edge.type = "BlockEdge")
                  g(itemsFactors, edge.type = "FactorEdge")
                  g(itemsExtras, edge.type = "ExtraEdge")
                }
            }
            "moveEdgesToVertex" <- function(X, v, edge.type = "VertexEdge") {
                edges <- edgeItem(v, edge.type = edge.type)
                if (length(edges) > 0) 
                  for (e in edges) if (!(is.null(e))) 
                    setEdgeCoords(e, edge.type = "vertexEdge", 
                      f = v, posFrom = X, raise = FALSE, setEdgeLabel = FALSE)
            }
            "moveVerticesInBlock" <- function(i, dxy, move.vertices = TRUE) {
                for (v in seq(along = vertexList)) 
                  if (is.element(v, dg@visibleVertices)) {
                  blockIndex <- retBlockIndex(v, vertex.type = "Vertex")
                  if ((blockIndex > 0) && (blockReferences[blockIndex] == i)) {
                    if (move.vertices) {
                      if (changeVertexPos(v, dxy)) 
                        if (!closedVertex[v]) 
                          tkmove(canvas, vertexItem(v)$tag, dxy[1], 
                            dxy[2])
                    }
                    pos <- retVertexPos(v, "Vertex")
                    moveEdgesToVertex(pos, v, edge.type = "VertexEdge")
                  }
                }
            }
            "subMoveVertex" <- function(i, vertex.type, posFrom, posTo) {
                dxy <- findDifference(posTo, posFrom)
                tag <- vertexItem(i, vertex.type)$tag
                if (setVertexPos(i, posTo, dxy, vertex.type)) {
                  tkmove(canvas, tag, dxy[1], dxy[2])
                  tkitemraise(canvas, tag)
                  if (control$debug.strata && (vertex.type != 
                    "Factor") && (vertex.type != "Extra") && 
                    (vertex.type != "ClosedBlock")) {
                    strata <- retStratum(i, vertex.type)
                    block <- retBlockIndex(i, vertex.type)
                    color <- myColor(strata)
                    tkitemconfigure(canvas, itemsVertices[[i]]$numbers, 
                      text = paste(i, strata, block, sep = "."))
                    tkitemconfigure(canvas, itemsVertices[[i]]$numbers, 
                      fill = color)
                  }
                  if (vertex.type == "ClosedBlock") {
                    moveVerticesInBlock(i, dxy, move.vertices = FALSE)
                    moveEdgesToVertex(posTo, -i, edge.type = "BlockEdge")
                  }
                  else if (vertex.type == "Vertex") {
                    if (closedVertex[i]) 
                      posTo <- retVertexPos(i, "Vertex")
                    moveEdgesToVertex(posTo, i, edge.type = "VertexEdge")
                  }
                  else if (vertex.type == "Extra") {
                    posTo <- retVertexPos(i, vertex.type)
                    moveEdgesToVertex(posTo, -i, edge.type = "ExtraEdge")
                  }
                  else if (vertex.type == "Factor") {
                    moveEdgesToVertex(posTo, i, edge.type = "FactorEdge")
                  }
                  else {
                    posTo <- retVertexPos(i, vertex.type)
                    moveEdgesToVertex(posTo, i, edge.type = "VertexEdge")
                  }
                }
            }
            "setFactorVertexPosition" <- function(i, vertex.indices) {
                if (!dg@factorVertexList[[abs(i)]]@fixed.positions) {
                  posFrom <- retVertexPos(-i, "Factor")
                  positions <- NULL
                  for (j in seq(along = vertex.indices)) if (vertex.indices[j] > 
                    0) 
                    positions <- rbind(positions, positionsVertices[vertex.indices[j], 
                      ])
                  position <- apply(positions, 2, mean)
                  posTo <- positionsCanvas(project(position))
                  subMoveVertex(-i, "Factor", posFrom, posTo)
                }
            }
            "moveFactorVertex" <- function(vertex) {
                if (!is.null(itemsFactors)) 
                  for (i in seq(along = itemsFactors)) if (!is.null(itemsFactors[[i]]) && 
                    !is.null(itemsFactors[[i]][[1]])) {
                    vertex.indices <- dg@factorVertexList[[i]]@vertex.indices
                    if (is.element(vertex, vertex.indices)) 
                      setFactorVertexPosition(i, vertex.indices)
                  }
            }
            "moveVertex" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                force(i)
                force(vertex.type)
                function(x, y) {
                  deActivateVertex(i, retVertexColor(i, vertex.type), 
                    vertex.type)
                  posFrom <- retVertexPos(i, vertex.type)
                  posTo <- replaceXY(x, y, posFrom)
                  subMoveVertex(i, vertex.type, posFrom, posTo)
                  moveFactorVertex(i)
                  setUpdatePositions("")
                }
            }
            "moveVertexLabel" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                force(i)
                force(vertex.type)
                function(x, y) {
                  deActivateVertex(i, retVertexColor(i, vertex.type), 
                    vertex.type)
                  if (vertex.type == "ClosedBlock") {
                    posFrom <- retVertexPos(i, vertex.type) + 
                      retBlockLabelPos(i)
                    X <- replaceXY(x, y, posFrom)
                    dxy <- findDifference(X, posFrom)
                    tkcoords(canvas, vertexItem(i, vertex.type)$label, 
                      x, y)
                    setLabelPos(i, X, dxy, vertex.type)
                  }
                  else {
                    posFrom <- retLabelPos(i, vertex.type)
                    X <- replaceXY(x, y, posFrom)
                    dxy <- findDifference(X, posFrom)
                    tkmove(canvas, vertexItem(i, vertex.type)$label, 
                      dxy[1], dxy[2])
                    setLabelPos(i, X, dxy, vertex.type)
                  }
                  setUpdatePositions("")
                }
            }
            "changeVertexLabel" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                force(i)
                force(vertex.type)
                function(...) {
                  ReturnVal <- modalDialog("Label Entry", "Enter new label", 
                    retVertexLabel(i, vertex.type), top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  setVertexLabel(i, ReturnVal, vertex.type)
                  setUpdatePositions("")
                }
            }
            "deleteVertexLabel" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                force(i)
                force(vertex.type)
                function(...) {
                  setVertexLabel(i, "", vertex.type)
                  setUpdatePositions("")
                }
            }
            "changeEdgeLabel" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(...) {
                  if (retActivatedEdge() == i) {
                    edges <- edgeItem(f, edge.type = edge.type)
                    if (length(edges) > 0) 
                      for (e in edges) if (!(is.null(e))) 
                        if (e$nr == i) 
                          if (e$to == t) {
                            label <- paste(edge.names(i, edge.type = edge.type), 
                              collapse = "%")
                            ReturnVal <- modalDialog("Edge label Entry", 
                              "Enter new label", label, top = GW.top)
                            if (ReturnVal == "ID_CANCEL") 
                              return()
                            setEdgeLabel(e, ReturnVal, e$label.number, 
                              f = f, permanent = TRUE)
                          }
                  }
                }
            }
            "subDropVertex" <- function(i, vertex.type = "Vertex", 
                slave = TRUE, upd = TRUE) {
                if (control$debug.position) 
                  print(paste("subDropVertex", i, vertex.type))
                tkconfigure(canvas, cursor = "watch")
                tkfocus(GW.top)
                redraw <- FALSE
                if (!(dg@viewType == "Simple")) 
                  redraw <- TRUE
                vertexEdges <- copyCurrentEdges(omitEdges = vertex.in.edge(i, 
                  edge.type = "VertexEdge"), edge.type = "VertexEdge", 
                  copyProperties = FALSE)
                factorEdges <- copyCurrentEdges(omitEdges = vertex.in.edge(i, 
                  edge.type = "FactorEdge"), edge.type = "FactorEdge", 
                  copyProperties = FALSE)
                extraEdges <- copyCurrentEdges(omitEdges = vertex.in.edge(i, 
                  edge.type = "ExtraEdge"), edge.type = "ExtraEdge", 
                  copyProperties = FALSE)
                blockEdges <- copyCurrentEdges(omitEdges = vertex.in.edge(i, 
                  edge.type = "BlockEdge"), edge.type = "BlockEdge", 
                  copyProperties = FALSE)
                if (vertex.type == "selected") {
                  newEdges <- list(vertexEdges = vertexEdges, 
                    factorEdges = factorEdges, extraEdges = extraEdges, 
                    blockEdges = blockEdges)
                }
                else if (vertex.type == "Vertex") {
                  newEdges <- list(vertexEdges = vertexEdges, 
                    factorEdges = factorEdges, extraEdges = extraEdges, 
                    blockEdges = blockEdges)
                }
                else if (vertex.type == "Factor") 
                  newEdges <- list(vertexEdges = vertexEdges, 
                    factorEdges = factorEdges, extraEdges = NULL, 
                    blockEdges = NULL)
                else if (vertex.type == "Extra") 
                  newEdges <- list(vertexEdges = vertexEdges, 
                    factorEdges = NULL, extraEdges = extraEdges, 
                    blockEdges = NULL)
                else if (vertex.type == "closedBlock") 
                  newEdges <- list(vertexEdges = vertexEdges, 
                    factorEdges = NULL, extraEdges = NULL, blockEdges = NULL)
                if (vertex.type == "selected") {
                  redraw <- TRUE
                  message("Selected vertices not added to 'newEdges';")
                  message("Resulting edges should be returned from modifyModel!")
                }
                R <- NULL
                if (is.null(object)) 
                  R <- TRUE
                Arguments <- Args()
                visible.Vertices <- returnVisibleVertices()
                if (i != 0) 
                  visible.Vertices <- visible.Vertices[visible.Vertices != i]
                if (!is.null(object) && (control$hasMethods || 
                  hasMethod("modifyModel", class(object)))) {
                  if (i == 0) 
                    name <- ""
                  else name <- retVertexName(i, vertex.type)
                  R <- modifyModel(object, action = "dropVertex", 
                    name = name, index = i, type = vertex.type, 
                    newEdges = newEdges, visibleVertices = visible.Vertices, 
                    selectedNodes = selectedNodes, selectedEdges = selectedEdges, 
                    Arguments = Arguments)
                }
                if (!is.null(R)) {
                  objectAssign(R)
                  if (slave || redraw) 
                    drawResult(newEdges, R, slave, "dropVertex")
                  else {
                    if (any(slotNames(R$object) == ".title"))
                      tktitle(GW.top) <- R$object@.title
                    setVisibleVertices(visible.Vertices)
                    if (i != 0) 
                      subSubUndisplayVertex(i)
                    updateSelectedEdges(R, vertexEdges, edge.type = NULL)
                  }
                }
                else message("Null result in dropVertex")
                updateBlockEdges()
                tkconfigure(canvas, cursor = "arrow")
                if (control$variableFrame) {
                  index <- i
                  if ((get("type", GW.top$env$box$env) == "variableList")) {
                    if (!(slave || redraw)) {
                      tkdelete(GW.top$env$box, index - 1)
                      tkinsert(GW.top$env$box, index - 1, tdv(namesVertices[[index]]))
                    }
                    if (tkselectionForVisibleVertices) {
                      for (i in 1:length(vertexList)) {
                        tkselection.clear(GW.top$env$box, i - 1)
                      }
                      for (i in returnVisibleVertices()) {
                        tkselection.set(GW.top$env$box, i - 1)
                      }
                    }
                  }
                  else updateVertexInBlock(index, k = retStratum(index), 
                    visibleBefore = TRUE, visibleAfter = FALSE)
                }
            }
            "undisplayVertex" <- function(i, vertex.type = ifelse(i > 
                0, "Vertex", "Factor"), slave = TRUE) {
                force(i)
                force(slave)
                function(...) {
                  deActivateVertex(i, retVertexColor(i, vertex.type), 
                    vertex.type)
                  subDropVertex(i, vertex.type, slave)
                }
            }
            "undisplayBlock" <- function(i, descendants = TRUE, 
                slave = TRUE, update = TRUE) {
                "subUndisplayBlock" <- function(i) {
                  deleteBlock(i)
                  visible.Blocks <- returnVisibleBlocks()
                  visible.Blocks <- visible.Blocks[visible.Blocks != 
                    i]
                  setVisibleBlocks(visible.Blocks)
                }
                if (descendants) 
                  for (j in blockList[[i]]@descendants) if ((j != 
                    i) && (j != 0)) {
                    subUndisplayBlock(j)
                    if ((closedBlock[j] || hiddenBlock[j])) 
                      openBlock(j, update = FALSE)()
                  }
                subUndisplayBlock(i)
                if (update) 
                  subUpdateGraphWindow("undisplayBlock", raiseEdges = TRUE, 
                    updateEdges = TRUE)
            }
            "removeBlock" <- function(i, descendants = TRUE, slave = TRUE, 
                update = TRUE) {
                "subRemoveBlock" <- function(i) {
                  deleteBlock(i)
                  visible.Blocks <- returnVisibleBlocks()
                  visible.Blocks <- visible.Blocks[visible.Blocks != i]
                  setVisibleBlocks(visible.Blocks)
                  positionsBlocks[i, , 1] <<- rep(-1000, local.N)
                  positionsBlocks[i, , 2] <<- rep(-1000, local.N) + 1e-04
                  positionsClosedBlocks[i, ] <<- rep(-1000, local.N)
                }
                if (descendants) 
                  for (j in blockList[[i]]@descendants) if ((j != i) && (j != 0)) {
                    subRemoveBlock(j)
                    if ((closedBlock[j] || hiddenBlock[j])) 
                      openBlock(j, update = FALSE)()
                  }
                subRemoveBlock(i)
                V <- sinkVertexList()
                if (updateAllBlockIndices() || TRUE) 
                  setUpdateBlockEdges("removeBlock")
                subUpdateGraphWindow("removeBlock", redrawVertices = TRUE, 
                  raiseEdges = TRUE, updateEdges = TRUE, all.blockframes = TRUE)
                setUpdateAll("removeBlock")
            }
            "markVerticesOfBlock" <- function(i, descendants = TRUE, 
                slave = TRUE, update = TRUE) {
                "subMarkVerticesOfBlock" <- function(i) {
                  subActivateVertex(i, "ClosedBlock", hit.type = "Mark-2", 
                    color = "GreenYellow")
                  for (v in seq(along = vertexList)) if (is.element(v, 
                    dg@visibleVertices)) {
                    bi <- retBlockIndex(v, vertex.type = "Vertex")
                    if ((bi == i)) {
                      subActivateVertex(v, "Vertex", hit.type = "Mark-2", 
                        color = "GreenYellow")
                      vertices <<- c(vertices, v)
                    }
                  }
                }
                vertices <- numeric(0)
                if (descendants) 
                  for (j in blockList[[i]]@descendants) if ((j != 
                    i) && (j != 0)) {
                    subMarkVerticesOfBlock(j)
                  }
                subMarkVerticesOfBlock(i)
                return(vertices)
            }
            "moveEdge" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(x, y) {
                  if (control$debug.position) 
                    print(paste("moveEdge", i, f, t, edge.type))
                  activateEdge(0, from = f, to = t, edge.type = edge.type)()
                  from.type <- vertexTypeOfEdge(f, edge.type)
                  to.type <- vertexTypeOfEdge(t, edge.type)
                  posFrom <- retVertexPos(f, from.type)
                  posTo <- retVertexPos(t, to.type)
                  pos <- (posFrom + posTo)/2
                  X <- replaceXY(x, y, pos)
                  dxy <- findDifference(X, pos)
                  if (sum(dxy^2) < 0.25^2 * sum((posFrom - posTo)^2)) {
                    "fun" <- function(ii, vertex.type, position) {
                      if (control$debug.position) 
                        print(paste("moveEdge, f", ii, vertex.type))
                      if (vertex.type == "ClosedBlock") {
                        if (setVertexPos(-ii, position + dxy, 
                          dxy, vertex.type = vertex.type)) 
                          if ((closedBlock[-ii]) || hiddenBlock[-ii]) 
                            tkmove(canvas, vertexItem(-ii, vertex.type = vertex.type)$tag, 
                              dxy[1], dxy[2])
                      }
                      else if (vertex.type == "Extra") {
                        if (setVertexPos(ii, position + dxy, 
                          dxy, vertex.type = vertex.type)) 
                          tkmove(canvas, vertexItem(ii, vertex.type = vertex.type)$tag, 
                            dxy[1], dxy[2])
                      }
                      else if (!closedVertex[ii]) {
                        if (setVertexPos(ii, position + dxy, 
                          dxy, vertex.type = vertex.type)) 
                          tkmove(canvas, vertexItem(ii)$tag, 
                            dxy[1], dxy[2])
                      }
                      edges <- edgeItem(ii, edge.type = edge.type)
                      if (length(edges) > 0) 
                        for (e in edges) if (!(is.null(e))) 
                          setEdgeCoords(e, edge.type = edge.type, 
                            f = ii, posFrom = position, raise = FALSE, 
                            setEdgeLabel = FALSE)
                    }
                    fun(f, from.type, posFrom)
                    fun(t, to.type, posTo)
                  }
                  setUpdatePositions("")
                }
            }
            "moveEdgeLabel" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(x, y) {
                  activateEdge(0, from = f, to = t, edge.type = edge.type)()
                  edges <- edgeItem(f, edge.type = edge.type)
                  if (length(edges) > 0) 
                    for (e in edges) if (!(is.null(e))) 
                      if (e$nr == i) 
                        if ((e$to == t)) 
                          if (!(e$reverse) || !Oriented) {
                            from.type <- vertexTypeOfEdge(f, 
                              edge.type)
                            to.type <- vertexTypeOfEdge(t, edge.type)
                            pos <- (retVertexPos(t, to.type) + 
                              retVertexPos(f, from.type))/2
                            pos <- pos + retEdgeLabelPos(e$label.number, 
                              f, e$to)
                            X <- replaceXY(x, y, pos)
                            dxy <- findDifference(X, pos)
                            # tkcoords(canvas, e$label, x, y)
                            tkmove(canvas, e$label, dxy[1], dxy[2])
                            setEdgeLabelPos(e, e$label.number, 
                              X, dxy, f = f)
                          }
                }
            }
            "subSubAddEdge" <- function(f, t, from.type, to.type, 
                edge.type = "VertexEdge", slave = TRUE, edgeClass = NULL) {
                fun <- function(result, ff, tt, edge.type = "VertexEdge") {
                  if (!any(which.edge(c(ff, tt), edge.type = edge.type)) && 
                    (!any(which.edge(c(tt, ff), edge.type = edge.type)) || 
                      dg@oriented)) 
                    return(append(result, list(c(ff, tt))))
                  else return(result)
                }
                redraw <- FALSE
                if (!(dg@viewType == "Simple")) 
                  redraw <- TRUE
                result <- NULL
                if (edge.type == "selected") {
                  redraw <- TRUE
                  message("Selected edges not added to 'newEdges';")
                  message("Resulting edges should be returned from modifyModel!")
                }
                else if (edge.type == "VertexEdge") {
                  result <- fun(result, f, t, edge.type = "VertexEdge")
                }
                else if (edge.type == "extraEdge") {
                  message("No edges to extra vertices!")
                }
                else if (edge.type == "FactorEdge") {
                  message("Resulting edges should be returned from modifyModel!")
                  redraw <- TRUE
                }
                else if (edge.type == "ExtraEdge") {
                  message("Resulting edges should be returned from modifyModel!")
                  redraw <- TRUE
                }
                else if (edge.type == "factorBlockEdge") {
                  message("Resulting edges should be returned from modifyModel!")
                  message("Factoredges are only draw between vertices and factors!")
                  redraw <- TRUE
                }
                else if (edge.type == "BlockEdge") {
                  from.block <- f
                  if (!(from.type == "Vertex")) {
                    from.block <- unique(sort(c(from.block, blockList[[from.block]]@descendants)))
                    from.block <- from.block[from.block != 0]
                  }
                  to.block <- t
                  if (!(to.type == "Vertex")) {
                    to.block <- unique(sort(c(to.block, blockList[[to.block]]@descendants)))
                    to.block <- to.block[to.block != 0]
                  }
                  if (from.type == "Vertex") {
                    for (w in seq(along = vertexList)) if (is.element(w, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(w, vertex.type = "Vertex"), 
                        to.block)) 
                        result <- fun(result, f, w, edge.type = "VertexEdge")
                  }
                  else if (to.type == "Vertex") {
                    for (v in seq(along = vertexList)) if (is.element(v, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(v, vertex.type = "Vertex"), 
                        from.block)) 
                        result <- fun(result, v, t, edge.type = "VertexEdge")
                  }
                  else {
                    for (v in seq(along = vertexList)) if (is.element(v, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(v, vertex.type = "Vertex"), 
                        from.block)) 
                        for (w in seq(along = vertexList)) if (is.element(w, 
                          dg@visibleVertices)) 
                          if (is.element(retBlockIndex(w, vertex.type = "Vertex"), 
                            to.block)) 
                            result <- fun(result, v, w, edge.type = "VertexEdge")
                  }
                }
                if (redraw || !is.null(result) || (edge.type == "selected")) {
                  vertexEdges <- appendToCurrentEdges(omitEdges = FALSE, 
                    new.edge = result, edge.type = "VertexEdge", 
                    edgeClass = edgeClass)
                  newEdges <- list(vertexEdges = vertexEdges, 
                    extraEdges = NULL, factorEdges = NULL, blockEdges = NULL)
                  R <- NULL
                  if (is.null(object)) 
                    R <- TRUE
                  Arguments <- Args()
                  if (!is.null(object) && (control$hasMethods || 
                    hasMethod("modifyModel", class(object)))) {
                    if (f == 0) 
                      name.f <- ""
                    else name.f <- retVertexName(f, from.type)
                    if (t == 0) 
                      name.t <- ""
                    else name.t <- retVertexName(t, to.type)
                    R <- modifyModel(object, action = "addEdge", 
                      name.1 = name.f, name.2 = name.t, from = f, 
                      to = t, from.type = from.type, to.type = to.type, 
                      newEdges = newEdges, selectedNodes = selectedNodes, 
                      selectedEdges = selectedEdges, Arguments = Arguments)
                  }
                  if (!is.null(R)) {
                    objectAssign(R)
                    if (slave || redraw) 
                      drawResult(newEdges, R, slave, "addEdge")
                    else {
                      if (any(slotNames(R$object) == ".title"))
                        tktitle(GW.top) <- R$object@.title
                      for (i in seq(along = result)) {
                        E <- append.index.edge(result[[i]], edge.type = "VertexEdge", 
                          edgeClass = edgeClass)
                        drawEdge(E[[length(E)]], length(E), lower = TRUE, 
                          edge.type = "VertexEdge")
                      }
                      updateSelectedEdges(R, vertexEdges, edge.type = edge.type, 
                        setVertex = TRUE, updateBE = TRUE)
                    }
                  }
                  else message("Null result in addEdge")
                }
            }
            "subAddEdge" <- function(f, t, from.type, to.type, 
                edge.type = "VertexEdge", slave = TRUE, edgeClass = NULL) {
                tkconfigure(canvas, cursor = "watch")
                tkfocus(GW.top)
                if ((from.type != "Extra") && (to.type != "Extra")) 
                  subSubAddEdge(f, t, from.type, to.type, edge.type = edge.type, 
                    slave = slave, edgeClass = edgeClass)
                if (f != 0) 
                  setVertexColor(f, color = retVertexColor(f, 
                    from.type), from.type)
                if (t != 0) 
                  setVertexColor(t, color = retVertexColor(t, 
                    to.type), to.type)
                setActivatedVertex(0, "Vertex")
                tkconfigure(canvas, cursor = "arrow")
            }
            "subNewBlock" <- function(position, get.name = TRUE, 
                color = control$blockColors[1]) {
                positionCenter <- apply(position, c(2), mean)
                n <- length(blockList) + 1
                if (get.name) {
                  label <- modalDialog("Name Entry", "Enter name of new block", 
                    paste("Block", n, sep = ""), top = GW.top)
                  if (label == "ID_CANCEL") 
                    return()
                }
                else label <- paste("B", n, sep = "")
                stratum <- n
                if (is.null(color)) 
                  if (is.null(blockList)) 
                    color <- "grey"
                  else color <- color(blockList[[1]])
                ancestors <- 0
                for (j in seq(along = blockList)) if (inBlock(positionCenter, 
                  j)) 
                  ancestors <- unique(sort(c(ancestors, j, ancestors(blockList[[j]]))))
                parent <- max(ancestors)
                ancestors <- unique(sort(c(ancestors(blockList[[parent]]), 
                  parent)))
                ancestors <- ancestors[ancestors != 0]
                new.block <- new("dg.Block", stratum = stratum, 
                  index = -stratum, position = position, color = color, 
                  closed = FALSE, visible = TRUE, label = label, 
                  parent = parent, ancestors = ancestors)
                sinkBlockList()
                for (j in seq(along = blockList)) if (j %in% 
                  ancestors) {
                  x <- unique(sort(c(blockList[[j]]@descendants, 
                    n)))
                  blockList[[j]]@descendants <<- x[x != 0]
                }
                children <- unique(sort(c(blockList[[parent]]@children, 
                  n)))
                blockList[[parent]]@children <<- children[children != 
                  0]
                if (is.null(blockList)) {
                  blockList <<- list(new.block)
                  class(blockList) <<- "dg.BlockList"
                  dgm.frameModels@blocks <<- blockList
                  positionsBlocks <- Positions(blockList)
                  d <- dim(positionsBlocks)
                  positionsBlocks <<- array(positionsBlocks, 
                    dim = c(d[1], d[2]/2, 2))
                  positionsBlockLabels <<- matrix(labelPosition(new.block), 
                    nrow = 1)
                  positionsClosedBlocks <<- matrix(positionCenter, 
                    nrow = 1)
                  blockReferences <<- c(n)
                  blockLabels <<- c(label)
                  strataBlocks <<- c(n)
                  closedBlock <<- c(FALSE)
                  hiddenBlock <<- c(FALSE)
                  openTreeBlock <<- c(FALSE)
                  dg@visibleBlocks <<- c(n)
                }
                else {
                  blockList <<- append(blockList, list(new.block))
                  class(blockList) <<- "dg.BlockList"
                  dgm.frameModels@blocks <<- blockList
                  d <- dim(positionsBlocks)
                  pBs <- array(positionsBlocks, dim = c(d[1], 
                    d[2] * 2))
                  pBs <- rbind(pBs, c(position(new.block)))
                  positionsBlocks <<- array(pBs, dim = c(d[1] + 
                    1, d[2], 2))
                  positionsBlockLabels <<- rbind(positionsBlockLabels, 
                    labelPosition(new.block))
                  positionsClosedBlocks <<- rbind(positionsClosedBlocks, 
                    positionCenter)
                  blockReferences <<- c(blockReferences, n)
                  blockLabels <<- c(blockLabels, label)
                  strataBlocks <<- c(strataBlocks, n)
                  closedBlock <<- c(closedBlock, FALSE)
                  hiddenBlock <<- c(hiddenBlock, FALSE)
                  openTreeBlock <<- c(openTreeBlock, FALSE)
                  dg@visibleBlocks <<- c(dg@visibleBlocks, n)
                  itemsClosedBlocks <<- append(itemsClosedBlocks, 
                    list(NULL))
                  itemsBlockEdges <<- append(itemsBlockEdges, 
                    list(NULL))
                }
                if (control$variableFrame && !.IsEmpty(blockList)) {
                  if (!(get("type", GW.top$env$box$env) == "variableList")) {
                    if (control$debug.strata) 
                      print(n)
                    if (control$debug.strata) 
                      print(ancestors)
                    m <- max(ancestors)
                    if (m == 0) 
                      parent = "root"
                    else {
                      parent <- ubl(block = blockList[[m]])
                    }
                    tkinsert.block(GW.top$env$box, parent = parent, 
                      block = blockList[[n]], open = openTreeBlock[n], 
                      delete = FALSE)
                  }
                  else message("Redraw graph (make slave window) for panel with block tree!")
                }
            }
            "new.Block" <- function(position, get.name = FALSE) {
                tkconfigure(canvas, cursor = "watch")
                tkconfigure(GW.top$env$viewLabel, text = paste(dg@viewType, 
                  " | Working !!!"))
                subNewBlock(position, get.name = get.name)
                n <- length(blockList)
                drawBlock(blockList[[n]], n)
                if (updateAllBlockIndices()) 
                  setUpdateBlockEdges("createNewBlock")
                subUpdateGraphWindow("createNewBlock", redrawVertices = TRUE, 
                  raiseEdges = TRUE, updateEdges = TRUE)
                setUpdateBlocks("")
                subUpdatePositions()
                tkconfigure(GW.top$env$viewLabel, text = dg@viewType)
                tkconfigure(canvas, cursor = "arrow")
            }
            "newBlockXY" <- function() {
                function(x, y) {
                  X <- replaceXY(x, y, rep(50, local.N))
                  position <- c(inversProject(inversCanvasPosition(X)))
                  delta <- c(10, 10, rep(50, local.N - 2))
                  position <- matrix(c(position - delta, position + 
                    delta), nrow = 2, byrow = TRUE)
                  new.Block(position, get.name = FALSE)
                }
            }
            "newEdge" <- function(i, vertex.type = ifelse(i > 0, 
                "Vertex", "Factor"), slave = TRUE, selectClass = FALSE) {
                force(i)
                force(slave)
                force(vertex.type)
                force(selectClass)
                function(...) {
                  edgeClass = control$edgeClasses[, 1][[1]]
                  if (selectClass) {
                    ReturnVal <- selectDialog("Select class of edge", 
                      "Select class", control$edgeClasses[, 1], 
                      top = GW.top)
                    if ((length(ReturnVal) > 0) && (ReturnVal != 
                      "ID_CANCEL")) 
                      edgeClass <- control$edgeClasses[ReturnVal, 
                        1]
                  }
                  from.type <- retActivatedVertexVertex.Type()
                  f <- retActivatedVertex()
                  t <- i
                  if (!subActivateVertex(i, "green", vertex.type, 
                    new.edge = TRUE)) 
                    if ((f != 0) && !((f == i) && (from.type == 
                      vertex.type))) {
                      edge.type <- "VertexEdge"
                      if ((vertex.type == "Extra") || (from.type == 
                        "Extra")) 
                        edge.type <- "extraEdge"
                      if ((vertex.type == "Factor") || (from.type == 
                        "Factor")) 
                        if ((edge.type == "extraEdge")) 
                          edge.type <- "factorExtraEdge"
                        else edge.type <- "FactorEdge"
                      if ((vertex.type == "ClosedBlock") || (from.type == 
                        "ClosedBlock")) 
                        if ((edge.type == "extraEdge")) 
                          edge.type <- "blockExtraEdge"
                        else if ((edge.type == "FactorEdge")) 
                          edge.type <- "factorBlockEdge"
                        else edge.type <- "BlockEdge"
                      subAddEdge(f, t, from.type, vertex.type, 
                        edge.type = edge.type, slave = slave, 
                        edgeClass = edgeClass)
                    }
                }
            }
            "addLastEdge" <- function(vertex.type = ifelse(i > 
                0, "Vertex", "Factor")) {
                function() {
                  n <- length(vertexList)
                  f <- retActivatedVertex()
                  if ((f > 0) && (retActivatedVertexVertex.Type() == 
                    "Vertex")) 
                    setVertexColor(f, color = retVertexColor(i, 
                      vertex.type), vertex.type)
                  setActivatedVertex(n - 1, vertex.type)
                  newEdge(n, vertex.type = "Vertex", slave = FALSE)()
                }
            }
            "subSubDropEdge" <- function(i, f, t, from.type = vertexTypeOfEdge(f, 
                edge.type), to.type = vertexTypeOfEdge(t, edge.type), 
                edge.type = "VertexEdge", slave = TRUE) {
                redraw <- FALSE
                if (!(dg@viewType == "Simple")) 
                  redraw <- TRUE
                j.g <- FALSE
                if (edge.type == "selected") {
                  redraw <- TRUE
                  message("Selected edges not added to 'newEdges';")
                  message("Resulting edges should be returned from modifyModel!")
                  vertexEdges <- copyCurrentEdges(edge.type = "VertexEdge")
                  newEdges <- list(vertexEdges = vertexEdges, 
                    extraEdges = NULL, factorEdges = NULL, blockEdges = NULL)
                }
                else if (edge.type == "VertexEdge") {
                  j.g <- which.unordered.edge(c(t, f), edge.type = edge.type)
                  vertexEdges <- copyCurrentEdges(omitEdges = j.g, 
                    edge.type = edge.type)
                  newEdges <- list(vertexEdges = vertexEdges, 
                    extraEdges = NULL, factorEdges = NULL, blockEdges = NULL)
                }
                else if (edge.type == "FactorEdge") {
                  message("Resulting edges should be returned from modifyModel!")
                  redraw <- TRUE
                  newEdges <- list(vertexEdges = NULL, extraEdges = NULL, 
                    factorEdges = NULL, blockEdges = NULL)
                }
                else if (edge.type == "ExtraEdge") {
                  message("Resulting edges should be returned from modifyModel!")
                  redraw <- TRUE
                  vertexEdges <- copyCurrentEdges(edge.type = "VertexEdge")
                  newEdges <- list(vertexEdges = NULL, extraEdges = NULL, 
                    factorEdges = NULL, blockEdges = NULL)
                }
                else if ((edge.type == "blockExtraEdge") || (edge.type == 
                  "factorExtraEdge") || (edge.type == "factorBlockEdge")) {
                  message("Not possible: Factoredges are only draw between vertices and factors!")
                  newEdges <- list(vertexEdges = NULL, extraEdges = NULL, 
                    factorEdges = NULL, blockEdges = NULL)
                }
                else if (edge.type == "BlockEdge") {
                  from.block <- -f
                  if (from.block > 0) 
                    from.block <- c(from.block, blockList[[from.block]]@descendants)
                  to.block <- -t
                  if (to.block > 0) 
                    to.block <- c(to.block, blockList[[to.block]]@descendants)
                  if (from.type == "Vertex") {
                    for (w in seq(along = vertexList)) if (is.element(w, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(w, vertex.type = "Vertex"), 
                        to.block)) 
                        j.g <- j.g | which.unordered.edge(c(f, 
                          w), edge.type = "VertexEdge")
                  }
                  else if (to.type == "Vertex") {
                    for (v in seq(along = vertexList)) if (is.element(v, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(v, vertex.type = "Vertex"), 
                        from.block)) 
                        j.g <- j.g | which.unordered.edge(c(v, 
                          t), edge.type = "VertexEdge")
                  }
                  else {
                    for (v in seq(along = vertexList)) if (is.element(v, 
                      dg@visibleVertices)) 
                      if (is.element(retBlockIndex(v, vertex.type = "Vertex"), 
                        from.block)) 
                        for (w in seq(along = vertexList)) if (is.element(w, 
                          dg@visibleVertices)) 
                          if (is.element(retBlockIndex(w, vertex.type = "Vertex"), 
                            to.block)) 
                            if (v != w) 
                              j.g <- j.g | which.unordered.edge(c(v, 
                                w), edge.type = "VertexEdge")
                  }
                  vertexEdges <- copyCurrentEdges(omitEdges = j.g, 
                    edge.type = "VertexEdge")
                  newEdges <- list(vertexEdges = vertexEdges, 
                    extraEdges = NULL, factorEdges = NULL, blockEdges = NULL)
                }
                R <- NULL
                if (is.null(object)) 
                  R <- TRUE
                Arguments <- Args()
                if (!is.null(object) && (control$hasMethods || 
                  hasMethod("modifyModel", class(object)))) {
                  if (f == 0) 
                    name.f <- ""
                  else name.f <- retVertexName(f, from.type)
                  if (t == 0) 
                    name.t <- ""
                  else name.t <- retVertexName(t, to.type)
                  R <- modifyModel(object, action = "dropEdge", 
                    name.1 = name.f, name.2 = name.t, from = f, 
                    to = t, from.type = from.type, to.type = to.type, 
                    edge.index = i, newEdges = newEdges, selectedNodes = selectedNodes, 
                    selectedEdges = selectedEdges, Arguments = Arguments)
                }
                if (!is.null(R)) {
                  objectAssign(R)
                  if (slave || redraw) 
                    drawResult(newEdges, R, slave, "dropEdge")
                  else {
                    if (any(slotNames(R$object) == ".title"))
                      tktitle(GW.top) <- R$object@.title
                    subSubDeleteEdge(i, f, t, edge.type = edge.type)
                    clearEdge(i, edge.type = edge.type)
                    if (edge.type == "BlockEdge") {
                      for (k in seq(along = j.g)) if (j.g[k]) {
                        subSubDeleteEdge(k, f, t, edge.type = "VertexEdge")
                        clearEdge(k, edge.type = "VertexEdge")
                      }
                    }
                    updateSelectedEdges(R, vertexEdges, edge.type = edge.type, 
                      setVertex = FALSE, updateBE = TRUE)
                  }
                }
                else message("Null result in dropEdge")
            }
            "subDropEdge" <- function(i, f, t, from.type = vertexTypeOfEdge(f, 
                edge.type), to.type = vertexTypeOfEdge(t, edge.type), 
                from.all = FALSE, to.all = FALSE, edge.type = "VertexEdge", 
                slave = TRUE) {
                tkconfigure(canvas, cursor = "watch")
                tkfocus(GW.top)
                if (from.type == "ClosedBlock") 
                  from.block <- f
                else from.block <- retBlockIndex(f, vertex.type = from.type)
                if (to.type == "ClosedBlock") 
                  to.block <- t
                else to.block <- retBlockIndex(t, vertex.type = to.type)
                if (from.all && ((from.type == "ClosedBlock") || 
                  ((from.type == "Vertex")))) {
                  for (v in seq(along = vertexList)) if (is.element(v, 
                    dg@visibleVertices)) 
                    if ((retBlockIndex(v, vertex.type = "Vertex") == 
                      from.block)) 
                      subDropEdge(NULL, v, t, "Vertex", to.type, 
                        FALSE, to.all, edge.type = edge.type, 
                        slave = slave)
                }
                else if (to.all && ((to.type == "ClosedBlock") || 
                  ((to.type == "Vertex")))) {
                  for (w in seq(along = vertexList)) if (is.element(w, 
                    dg@visibleVertices)) 
                    if ((retBlockIndex(w, vertex.type = "Vertex") == 
                      to.block)) 
                      subDropEdge(NULL, f, w, "Vertex", "Vertex", 
                        FALSE, FALSE, edge.type = edge.type, 
                        slave = slave)
                }
                else {
                  j <- which.unordered.edge(c(t, f), edge.type = edge.type)
                  if (any(j)) {
                    i <- (1:(length(j)))[j]
                    for (j in seq(along = i)) subSubDropEdge(i[j], 
                      f, t, edge.type = edge.type, slave = slave)
                  }
                }
                tkconfigure(canvas, cursor = "arrow")
            }
            "deleteEdge" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(...) {
                  if (retActivatedEdge() == i) 
                    subDropEdge(i, f, t, from.all = FALSE, to.all = FALSE, 
                      edge.type = edge.type, slave = FALSE)
                }
            }
            "changeEdgeClass" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(...) {
                  ReturnVal <- selectDialog("Test selectDialog Entry", 
                    "Select name", control$edgeClasses[, 1], 
                    top = GW.top)
                  if ((length(ReturnVal) > 0) && (ReturnVal != 
                    "ID_CANCEL")) {
                    E <- getEdges(edge.type = edge.type)
                    newEdge <- E[[i]]
                    class(newEdge) <- control$edgeClasses[ReturnVal, 
                      2]
                    E <- append.edge(newEdge, edge.type = edge.type)
                    subSubDeleteEdge(i, f, t, edge.type = edge.type)
                    clearEdge(i, edge.type = edge.type)
                    drawEdge(E[[length(E)]], length(E), lower = TRUE, 
                      edge.type = edge.type)
                    setModel(object, txt = "changeEdgeClass", 
                      copyProperties = TRUE, RR = NULL)
                  }
                }
            }
            "propertyEdge" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(...) {
                  E <- sinkEdgeList()
                  E <- sinkBlockEdges()
                  E <- sinkFactorEdgeList()
                  E <- sinkExtraEdgeList()
                  E <- getEdges(edge.type = edge.type)
                  fixedSlots <- c("vertex.indices", "label.position")
                  if (length(selectedEdges) > 0) 
                    fixedSlots <- c(fixedSlots, "label", "oriented")
                  ReturnVal <- propertyDialog(E[[i]], control$edgeClasses, 
                    fixedSlots = fixedSlots, top = GW.top)
                  if (!is.null(ReturnVal$values)) {
                    E[[i]] <- ReturnVal$object
                    subSubDeleteEdge(i, f, t, edge.type = edge.type)
                    drawEdge(E[[i]], i, lower = TRUE, edge.type = edge.type, 
                      newE = TRUE)
                    if (length(selectedEdges) > 0) {
                      lapply(selectedEdges, function(k) setEdgeColor(k$index, 
                        color = NULL, edge.type = k$edge.type))
                      if (!is.null(ReturnVal$values$oriented)) {
                        message("Change of 'oriented' not implemented; ")
                      }
                      if (!is.null(ReturnVal$values$dash)) 
                        lapply(selectedEdges, function(k) if (k$edge.type == 
                          edge.type) {
                          dash(E[[k$index]]) <<- ReturnVal$values$dash
                          setEdgeDash(k$index, dash = ReturnVal$values$dash, 
                            edge.type = k$edge.type)
                        })
                      if (!is.null(ReturnVal$values$vertex.indices)) {
                        message("Change of 'vertex.indices' not possible; ")
                      }
                      if (!is.null(ReturnVal$values$width)) 
                        lapply(selectedEdges, function(k) if (k$edge.type == 
                          edge.type) {
                          width(E[[k$index]]) <<- ReturnVal$values$width
                          edges <- edgeItem(k$from, edge.type = edge.type)
                          if (length(edges) > 0) 
                            for (e in edges) if (!(is.null(e))) 
                              if (e$nr == k$index) 
                                if (e$to == k$to) 
                                  setEdgeWidth(e, width = ReturnVal$values$width, 
                                    edge.type = k$edge.type)
                        })
                      if (!is.null(ReturnVal$values$color)) 
                        lapply(selectedEdges, function(k) if (k$edge.type == 
                          edge.type) {
                          color(E[[k$index]]) <<- ReturnVal$values$color
                          setEdgeColor(k$index, color = ReturnVal$values$color, 
                            edge.type = k$edge.type)
                        })
                      if (!is.null(ReturnVal$values$label)) {
                        message("Change of 'label' not possible for group of edges; ")
                      }
                      if (!is.null(ReturnVal$values$label.position)) {
                        message("Change of 'label.position' not possible for group of edges; ")
                      }
                      if (!is.null(ReturnVal$values$class)) 
                        lapply(selectedEdges, function(k) if (k$edge.type == 
                          edge.type) {
                          class(E[[k$index]]) <<- ReturnVal$values$class
                          subSubDeleteEdge(k$index, k$from, k$to, 
                            edge.type = k$edge.type)
                          drawEdge(E[[k$index]], k$index, lower = TRUE, 
                            edge.type = k$edge.type)
                        })
                      if (is.null(ReturnVal$values$color)) 
                        clearSelectedEdges()
                      else selectedEdges <<- list()
                    }
                    if (edge.type == "VertexEdge") 
                      GraphWindow@dg@edgeList <<- E
                    else if (edge.type == "FactorEdge") 
                      GraphWindow@dg@factorEdgeList <<- E
                    else if (edge.type == "ExtraEdge") 
                      GraphWindow@dg@extraEdgeList <<- E
                    else if (edge.type == "BlockEdge") 
                      GraphWindow@dg@blockEdgeList <<- E
                    setModel(object, txt = "changeEdgeClass", 
                      copyProperties = TRUE, RR = NULL)
                  }
                }
            }
            "propertyNode" <- function(i, vertex.type = "Vertex") {
                force(i)
                force(vertex.type)
                function(...) {
                  "f" <- function() {
                    if ((length(selectedNodes) > 0) && (length(ReturnVal$values) > 
                      1) && ("color" %in% names(ReturnVal$values))) 
                      message("Only color changed for selected vertices; ")
                    if ("name" %in% names(ReturnVal$values)) {
                      if ((vertex.type == "ClosedBlock") || (vertex.type == 
                        "OpenBlock")) {
                        message("No names for blocks! ; ")
                      }
                      else {
                        message("Names can not be changed, use labels ; ")
                      }
                    }
                    if ("index" %in% names(ReturnVal$values)) {
                      message("The internal index not be changed! ; ")
                    }
                    if ("ancestors" %in% names(ReturnVal$values)) {
                      message("Ancestors should not be set manually; ")
                    }
                    if ("descendants" %in% names(ReturnVal$values)) {
                      message("Descendants should not be set manually; ")
                    }
                    if ("position" %in% names(ReturnVal$values)) {
                      if ((vertex.type == "ClosedBlock") || (vertex.type == 
                        "OpenBlock")) {
                        positionsBlocks[i, , ] <<- ReturnVal$values$position
                        subUpdateGraphWindow("Update from main menu", 
                          redrawVertices = TRUE, raiseEdges = TRUE, 
                          updateEdges = TRUE, all.blockframes = TRUE)
                        setUpdateBlocks("")
                      }
                      else {
                        posFrom <- retVertexPos(i, vertex.type)
                        posTo <- positionsCanvas(ReturnVal$values$position)
                        subMoveVertex(i, vertex.type, posFrom, 
                          posTo)
                        moveFactorVertex(i)
                        setUpdatePositions("")
                      }
                    }
                    if ("blockindex" %in% names(ReturnVal$values)) {
                      if ((vertex.type == "ClosedBlock") || (vertex.type == 
                        "OpenBlock")) {
                      }
                      else {
                        message(paste("Blockindex of vertices only used when positions", 
                          " of vertices relative blocks are ignored ; "))
                        # if (FALSE && .IsEmpty(blockList)) 
                          # if (vertex.type == "Vertex") {
                          #   blockindexVertices[i] <<- ReturnVal$blockindex
                          # 'blockindexVertices' !!!
                          # }
                          # else if (vertex.type == "Factor") {
                          #   blockindexFactorVertices[-i] <<- ReturnVal$blockindex
                          # 'blockindexFactorVertices' !!!
                          # }
                          # else if (vertex.type == "Extra") {
                          #   blockindexExtraVertices[i] <<- ReturnVal$blockindex
                          # 'blockindexExtraVertices' !!!
                          # }
                      }
                    }
                    if ("stratum" %in% names(ReturnVal$values)) {
                      if ((vertex.type == "ClosedBlock") || (vertex.type == 
                        "OpenBlock")) {
                        strataBlocks[abs(i)] <<- ReturnVal$values$stratum
                        if (updateAllBlockIndices()) 
                          setUpdateBlockEdges("propertyDialog: Block stratum")
                        subUpdateGraphWindow("propertyDialog: Block stratum", 
                          updateEdges = TRUE, blockframes = 0)
                        setUpdateBlocks("")
                      }
                      else {
                        message(paste("Stratas of vertices only used when positions", 
                          " of vertices relative blocks are ignored ; "))
                        if (FALSE && .IsEmpty(blockList)) 
                          if (vertex.type == "Vertex") {
                            strataVertices[i] <<- ReturnVal$stratum
                          }
                          else if (vertex.type == "Factor") {
                            strataFactorVertices[-i] <<- ReturnVal$stratum
                          }
                          else if (vertex.type == "Extra") {
                            strataExtraVertices[i] <<- ReturnVal$stratum
                          }
                      }
                    }
                    if ("closed" %in% names(ReturnVal$values)) {
                      message("Closed block by mouse interaction; ")
                    }
                    if ("visible" %in% names(ReturnVal$values)) {
                      message("Closed ancestor block by mouse interaction; ")
                    }
                    if ("color" %in% names(ReturnVal$values)) {
                      if (vertex.type != "OpenBlock") 
                        setVertexColor(i, ReturnVal$values$color, 
                          vertex.type, permanent = TRUE)
                      else message("Close and open block to activate color; ")
                      if (length(selectedNodes) > 0) {
                        lapply(selectedNodes, function(k) if ((k$node.type != 
                          "OpenBlock") && (k$node.type != "ClosedBlock")) 
                          setVertexColor(k$index, color = ReturnVal$values$color, 
                            vertex.type = k$node.type, permanent = TRUE))
                        selectedNodes <<- list()
                      }
                      setUpdatePositions("")
                    }
                    if ("constrained" %in% names(ReturnVal$values)) {
                      constrainedVertices[i] <<- ReturnVal$values$constrained
                    }
                    if ("label" %in% names(ReturnVal$values)) {
                      setVertexLabel(i, ReturnVal$values$label, 
                        vertex.type)
                      setUpdatePositions("")
                    }
                    if ("label.position" %in% names(ReturnVal$values)) 
                      if (!(vertex.type == "OpenBlock")) {
                        posFrom <- retLabelPos(i, vertex.type)
                        X <- positionsCanvas(ReturnVal$values$label.position)
                        dxy <- findDifference(X, posFrom)
                        tkmove(canvas, vertexItem(i, vertex.type)$label, 
                          dxy[1], dxy[2])
                        setLabelPos(i, X, dxy, vertex.type)
                        setUpdatePositions("")
                      }
                    if ("class" %in% names(ReturnVal$values)) {
                      message("Changing class of variable!")
                      subUpdateGraphWindow("Update from main menu", 
                        redrawVertices = TRUE, raiseEdges = TRUE, 
                        updateEdges = TRUE, all.blockframes = TRUE)
                      setUpdateVertices("")
                    }
                  }
                  if ((vertex.type == "ClosedBlock") || (vertex.type == 
                    "OpenBlock")) {
                    B <- sinkBlockList()
                    block <- blockList[[i]]
                    difficultSlots <- c("ancestors", "descendants", 
                      "closed", "visible")
                    if ((vertex.type == "OpenBlock")) 
                      difficultSlots <- c(difficultSlots, "color", 
                        "label", "label.position")
                    ReturnVal <- propertyDialog(block, NULL, 
                      okReturn = FALSE, fixedSlots = c("index"), 
                      difficultSlots = difficultSlots, top = GW.top)
                    if (!is.null(ReturnVal$values)) {
                      dgm.frameModels@blocks[[abs(i)]] <<- ReturnVal$object
                      blockList[[abs(i)]] <<- ReturnVal$object
                      f()
                    }
                  }
                  else if (vertex.type == "Vertex") {
                    V <- sinkVertexList()
                    vertex <- dgm.frameModels@vertices[[i]]
                    ReturnVal <- propertyDialog(vertex, control$vertexClasses, 
                      fixedSlots = c("index", "name"), difficultSlots = c("blockindex", 
                        "stratum"), okReturn = FALSE, top = GW.top)
                    if (!is.null(ReturnVal$values)) {
                      dgm.frameModels@vertices[[i]] <<- ReturnVal$object
                      vertexList[[i]] <<- ReturnVal$object
                      f()
                    }
                  }
                  else if (vertex.type == "Factor") {
                    V <- sinkFactorVertexList()
                    vertex <- GraphWindow@dg@factorVertexList[[abs(i)]]
                    ReturnVal <- propertyDialog(vertex, control$factorClasses, 
                      fixedSlots = c("name", "index", "vertex.indices"), 
                      difficultSlots = c("blockindex", "stratum"), 
                      okReturn = FALSE, top = GW.top)
                    if (!is.null(ReturnVal$values)) {
                      GraphWindow@dg@factorVertexList[[-i]] <<- ReturnVal$object
                      dg@factorVertexList[[-i]] <<- ReturnVal$object
                      f()
                    }
                  }
                  else if (vertex.type == "Extra") {
                    V <- sinkExtraVertexList()
                    vertex <- GraphWindow@dg@extraList[[i]]
                    ReturnVal <- propertyDialog(vertex, NULL, 
                      okReturn = FALSE, fixedSlots = c("name", 
                        "index"), difficultSlots = c("blockindex", 
                        "stratum"), top = GW.top)
                    if (!is.null(ReturnVal$values)) {
                      GraphWindow@dg@extraList[[i]] <<- ReturnVal$object
                      dg@extraList[[i]] <<- ReturnVal$object
                      f()
                    }
                  }
                }
            }
            "computeEdgeLabel" <- function(i, f, t, force = FALSE, 
                edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(force)
                force(edge.type)
                function(...) {
                  if (retActivatedEdge() == i) {
                    tkconfigure(canvas, cursor = "watch")
                    tkfocus(GW.top)
                    edges <- edgeItem(f, edge.type = edge.type)
                    if (length(edges) > 0) 
                      for (e in edges) if (!(is.null(e))) 
                        if (e$nr == i) 
                          if (e$to == t) 
                            if (!is.null(object) && (control$hasMethods || 
                              hasMethod("testEdge", class(object)))) {
                              from.type <- vertexTypeOfEdge(f, 
                                edge.type)
                              to.type <- vertexTypeOfEdge(t, 
                                edge.type)
                              R <- testEdge(object, action = "remove", 
                                name.1 = retVertexName(f, from.type), 
                                name.2 = retVertexName(t, to.type), 
                                from = f, to = t, from.type = from.type, 
                                to.type = to.type, edge.index = i, 
                                force = force, Arguments = Args())
                              if (!is.null(R)) {
                                if ((control$hasMethods || hasMethod("label", 
                                  class(R)))) 
                                  setEdgeLabel(e, label(R), e$label.number, 
                                    f = f, permanent = TRUE)
                                if ((control$hasMethods || hasMethod("width", 
                                  class(R)))) 
                                  setEdgeWidth(e, width(R), e$label.number, 
                                    f = f)
                                activateEdge(0, from = f, to = t, 
                                  edge.type = edge.type)()
                              }
                            }
                    tkconfigure(canvas, cursor = "arrow")
                  }
                }
            }
            "deleteEdgeLabel" <- function(i, f, t, edge.type = "VertexEdge") {
                force(i)
                force(f)
                force(t)
                force(edge.type)
                function(...) {
                  edges <- edgeItem(f, edge.type = edge.type)
                  E <- getEdges(edge.type = edge.type)[[i]]
                  if (length(edges) > 0) 
                    for (e in edges) if (!(is.null(e))) 
                      if (e$nr == i) {
                        setEdgeLabel(e, "", e$label.number, f = f, 
                          permanent = TRUE)
                        setEdgeWidth(e, E@width, e$label.number, 
                          f = f)
                      }
                }
            }
            "moveBlockPoint" <- function(i, A) {
                force(i)
                force(A)
                function(x, y) {
                  posTo <- retBlockPoints(i)[A, ]
                  X <- replaceXY(x, y, posTo)
                  dxy <- findDifference(X, posTo)
                  changeBlockCornerPos(i, A, dxy)
                  tkcoordsBlock(i, color = "black", lower = FALSE)
                  if (updateAllBlockIndices()) 
                    setUpdateBlockEdges("moveBlockPoint")
                  subUpdateGraphWindow("moveBlockPoint", blockframes = 0)
                  setUpdateBlocks("")
                }
            }
            "moveBlockLine" <- function(i, A, B) {
                force(i)
                force(A)
                force(B)
                function(x, y) {
                  if ((A == 1) && (B == 3)) 
                    direction <- 1
                  else if ((A == 4) && (B == 8)) 
                    direction <- 1
                  else if ((A == 1) && (B == 4)) 
                    direction <- 2
                  else if ((A == 3) && (B == 8)) 
                    direction <- 2
                  else if ((A == 5) && (B == 6)) 
                    direction <- 1
                  else if ((A == 7) && (B == 2)) 
                    direction <- 1
                  else if ((A == 5) && (B == 7)) 
                    direction <- 2
                  else if ((A == 6) && (B == 2)) 
                    direction <- 2
                  else if ((A == 1) && (B == 5)) 
                    direction <- 3
                  else if ((A == 3) && (B == 6)) 
                    direction <- 3
                  else if ((A == 4) && (B == 7)) 
                    direction <- 3
                  else if ((A == 8) && (B == 2)) 
                    direction <- 3
                  posA <- retBlockPoints(i)[A, ]
                  posB <- retBlockPoints(i)[B, ]
                  X <- replaceXY(x, y, (posA + posB)/2)
                  Y <- X
                  Y[direction] <- (posA[direction] + posB[direction])/2
                  dxy <- findDifference(X, Y)
                  changeBlockCornerPos(i, A, dxy)
                  tkcoordsBlock(i, color = "black", lower = FALSE)
                  if (updateAllBlockIndices()) 
                    setUpdateBlockEdges("moveBlockLine")
                  subUpdateGraphWindow("moveBlockLine", blockframes = 0)
                  setUpdateBlocks("")
                }
            }
            "moveBlock" <- function(i, A) {
                force(i)
                force(A)
                Y <- NULL
                function(x, y) {
                  posTo <- retBlockPoints(i)[A, ]
                  if (is.null(Y)) 
                    Y <<- replaceXY(x, y, posTo)
                  else {
                    X <- replaceXY(x, y, posTo)
                    dxy <- findDifference(X, Y)
                    Y <<- X
                    changeBlockPos(i, A, dxy)
                    moveVerticesInBlock(i, dxy)
                    for (j in blockList[[i]]@descendants) if ((j != 
                      0) && (j != i)) {
                      changeBlockPos(j, A, dxy)
                      moveVerticesInBlock(j, dxy)
                      if (!(hiddenBlock[j] || closedBlock[j])) 
                        tkcoordsBlock(j, color = "black", lower = FALSE)
                    }
                    if (updateAllBlockIndices()) 
                      setUpdateBlockEdges("moveBlock")
                    subUpdateGraphWindow("moveBlock", blockframes = 0)
                    tkcoordsBlock(i, color = "black", lower = FALSE)
                    setUpdateBlocks("")
                  }
                }
            }
            "hideBlock" <- function(i, ancestor, update = TRUE) {
                blockReferences[i] <<- ancestor
                for (v in seq(along = vertexList)) if (is.element(v, 
                  dg@visibleVertices)) 
                  if (retBlockIndex(v, vertex.type = "Vertex") == 
                    i) 
                    setCloseVertex(v, TRUE, "Vertex")
                setHiddenBlock(i, TRUE, update = update)
                for (j in blockList[[i]]@descendants) {
                  if ((j != 0) && !hiddenBlock[j]) 
                    hideBlock(j, ancestor, update = FALSE)
                }
            }
            "closeBlock" <- function(i, update = TRUE) {
                force(i)
                function(...) {
                  for (v in seq(along = vertexList)) if (is.element(v, 
                    dg@visibleVertices)) 
                    if (retBlockIndex(v, vertex.type = "Vertex") == 
                      i) 
                      setCloseVertex(v, TRUE, "Vertex")
                  for (j in blockList[[i]]@descendants) {
                    if ((j != i) && (j != 0) && !hiddenBlock[j]) 
                      hideBlock(j, i, update = FALSE)
                  }
                  setClosedBlock(i, TRUE, update = update)
                  drawVertex(i, w = 10, vertexcolor = "Black", 
                    vertex.type = "ClosedBlock")
                  if (update) 
                    subUpdateGraphWindow("closeBlock", blockframes = i, 
                      updateEdges = TRUE)
                }
            }
            "subOpenBlock" <- function(i, update = TRUE) {
                blockReferences[i] <<- i
                if ((hiddenBlock[i] && closedBlock[i])) {
                  drawVertex(i, w = 10, vertexcolor = "Black", 
                    vertex.type = "ClosedBlock")
                }
                else {
                  setClosedBlock(i, FALSE, update = update)
                  vertex.type <- "ClosedBlock"
                  deActivateVertex(i, retVertexColor(i, vertex.type), 
                    vertex.type)
                  if (is.element(i, dg@visibleBlocks)) 
                    drawBlock(blockList[[i]], i)
                  for (v in seq(along = vertexList)) if (is.element(v, 
                    dg@visibleVertices)) {
                    bi <- retBlockIndex(v, vertex.type = "Vertex")
                    if ((bi == i) || (bi == 0)) 
                      setCloseVertex(v, FALSE, "Vertex")
                  }
                  for (j in blockList[[i]]@descendants) if ((j != 
                    0) && (hiddenBlock[j])) 
                    if (!isInClosedBlock(j)) 
                      subOpenBlock(j, update = FALSE)
                }
                setHiddenBlock(i, FALSE, update = FALSE)
            }
            "openBlock" <- function(i, update = TRUE) {
                force(i)
                force(update)
                function(...) {
                  subOpenBlock(i, update)
                  if (update) 
                    subUpdateGraphWindow("openBlock", raiseEdges = TRUE, 
                      updateEdges = TRUE)
                }
            }
            "subNewVertex" <- function(position, get.name = TRUE, 
                get.vertex.type = TRUE, get.how.to.compute = TRUE) {
                n <- length(vertexList) + 1
                if (get.name) {
                  ReturnVal <- modalDialog("Name Entry", "Enter name of new variable", 
                    paste("V", n, sep = ""), top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                }
                else ReturnVal <- paste("V", n, sep = "")
                vertextypes <- control$vertexClasses[, 1]
                vertextypes <- paste(vertextypes)
                vertex.type <- vertextypes[1]
                if (get.vertex.type) {
                  Return.Val <- selectDialog("Variable vertex.type selection entry", 
                    "Select vertex.type", vertextypes, top = GW.top)
                  vertex.type <- vertextypes[Return.Val]
                }
                if (!is.na(vertex.type) && (vertex.type == "TextVertex")) 
                  prototype <- "dg.TextVertex"
                else prototype <- typeToPrototype(type = vertex.type, 
                  prototype = "dg.Vertex", classes = control$vertexClasses)
                vertex <- new(prototype, name = ReturnVal, label = ReturnVal, 
                  index = n, position = position, blockindex = 0, 
                  stratum = 0, color = control$vertexColor)
                if (get.how.to.compute) {
                  Expression <- modalDialog("Expression Entry", 
                    "Enter expression for computing new variable", 
                    paste("object <<- object"), top = GW.top)
                  if (Expression == "ID_CANCEL") 
                    return()
                  else eval(parse(text = Expression))
                }
                else Expression <- ReturnVal
                sinkVertexList()
                vertexList <<- append(vertexList, list(vertex))
                class(vertexList) <<- "dg.VertexList"
                dgm.frameModels@vertices <<- vertexList
                positionsVertices <<- rbind(positionsVertices, 
                  position(vertex))
                positionsLabels <<- rbind(positionsLabels, position(vertex))
                Labels <<- c(Labels, Expression)
                closedVertex <<- c(closedVertex, FALSE)
                constrainedVertices <<- c(constrainedVertices, 
                  control$constrained)
                colorsVertices <<- c(colorsVertices, control$vertexColor)
                blocksVertices <<- c(blocksVertices, 0)
                strataVertices <<- c(strataVertices, 0)
                namesVertices <<- c(namesVertices, ReturnVal)
                itemsVertices <<- append(itemsVertices, list(NULL))
                itemsEdges <<- append(itemsEdges, list(NULL))
                updateAllBlockIndices()
                if (control$variableFrame) {
                  if ((get("type", GW.top$env$box$env) == "variableList")) {
                    if (!is.element(n, dg@visibleVertices)) 
                      ReturnVal <- tdv(ReturnVal)
                    tkinsert(GW.top$env$box, "end", ReturnVal)
                    if (tkselectionForVisibleVertices) 
                      for (i in returnVisibleVertices()) {
                        tkselection.set(GW.top$env$box, i - 1)
                      }
                  }
                }
            }
            "subAddVertex" <- function(index, vertex.type = "Vertex", 
                slave = TRUE) {
                tkconfigure(canvas, cursor = "watch")
                tkfocus(GW.top)
                redraw <- FALSE
                if (!(dg@viewType == "Simple")) 
                  redraw <- TRUE
                vertexEdges <- appendToCurrentEdges(omitEdges = FALSE, 
                  new.edge = NULL, edge.type = "VertexEdge")
                newEdges <- list(vertexEdges = vertexEdges, 
                  factorEdges = getEdges(edge.type = "FactorEdge"), 
                  extraEdges = getEdges(edge.type = "ExtraEdge"), 
                  blockEdges = getEdges(edge.type = "BlockEdge"))
                if (vertex.type == "selected") {
                  redraw <- TRUE
                  message("Selected vertices not added to 'newEdges';")
                  message("Resulting edges should be returned from modifyModel!")
                }
                R <- NULL
                if (is.null(object)) 
                  R <- TRUE
                Arguments <- Args()
                visible.Vertices <- returnVisibleVertices()
                if (index != 0) 
                  visible.Vertices <- c(visible.Vertices, index)
                if (!is.null(object) && (control$hasMethods || 
                  hasMethod("modifyModel", class(object)))) {
                  if (index == 0) 
                    name <- ""
                  else name <- retVertexName(index, vertex.type)
                  R <- modifyModel(object, action = "addVertex", 
                    name = name, index = index, type = vertex.type, 
                    newEdges = newEdges, visibleVertices = visible.Vertices, 
                    selectedNodes = selectedNodes, selectedEdges = selectedEdges, 
                    Arguments = Arguments)
                }
                if (!is.null(R)) {
                  objectAssign(R)
                  if (slave || redraw) 
                    drawResult(newEdges, R, slave, "addVertex")
                  else {
                    if (any(slotNames(R$object) == ".title"))
                      tktitle(GW.top) <- R$object@.title
                    setVisibleVertices(visible.Vertices)
                    drawVertex(index, w = control$w, vertexcolor = control$vertexColor, 
                      vertex.type = "Vertex")
                    updateSelectedEdges(R, vertexEdges, edge.type = NULL)
                  }
                }
                else message("Null result in addVertex")
                if (control$variableFrame) {
                  if ((get("type", GW.top$env$box$env) == "variableList")) {
                    if (!(slave || redraw)) {
                      tkdelete(GW.top$env$box, index - 1)
                      tkinsert(GW.top$env$box, index - 1, namesVertices[[index]])
                    }
                    if (tkselectionForVisibleVertices) 
                      for (i in returnVisibleVertices()) {
                        tkselection.set(GW.top$env$box, i - 1)
                      }
                  }
                  else updateVertexInBlock(index, k = retStratum(index), 
                    visibleBefore = FALSE, visibleAfter = TRUE)
                }
                tkconfigure(canvas, cursor = "arrow")
            }
            "newVertexXY" <- function() {
                function(x, y) {
                  X <- replaceXY(x, y, rep(50, local.N))
                  position <- c(inversProject(inversCanvasPosition(X)))
                  subNewVertex(position, get.name = FALSE, get.vertex.type = FALSE, 
                    get.how.to.compute = FALSE)
                  n <- length(vertexList)
                  subAddVertex(n, slave = FALSE)
                  tkfocus(canvas)
                }
            }
            "zoom" <- function(xy, factor) {
                "f" <- function(title, x) print(paste(title, paste(x, 
                  collapse = ", ")))
                if (is.null(xy)) 
                  xy <- c(control$width, control$height, rep(100, 
                    local.N - 2))/2
                X <- xy
                if (is.null(xy)) {
                  xy <- c(control$width, control$height, rep(100, 
                    local.N - 2))/2
                }
                else if (!is.null(zoomCenter)) 
                  X <- (xy - zoomCenter)/Scale + zoomCenter
                Scale <<- Scale * factor
                xy <- round(xy)
                for (i in seq(along = GW.tags)) {
                  tkitemscale(canvas, GW.tags[[i]], xy[1], xy[2], 
                    factor, factor)
                }
                if (control$debug.position) {
                  f("xy: ", xy)
                  f("X: ", X)
                  f("Old zoomCenter: ", zoomCenter)
                  f("Old upperLeft: ", upperLeft)
                  f("Old lowerRight: ", lowerRight)
                }
                upperLeft <<- upperLeft - (upperLeft - X)/factor
                lowerRight <<- lowerRight - (lowerRight - X)/factor
                newCenter <- (lowerRight + upperLeft)/2
                if (!is.null(zoomCenter)) 
                  if ((sum(zoomCenter - newCenter)^2) > 0.1) 
                    zoomCenterSet <<- TRUE
                zoomCenter <<- newCenter
                if (control$debug.position) {
                  f("New zoomCenter: ", zoomCenter)
                  f("New upperLeft: ", upperLeft)
                  f("New lowerRight: ", lowerRight)
                }
                if (zoomCenterSet) 
                  subUpdateGraphWindow("Update from zoom", redrawVertices = TRUE, 
                    raiseEdges = TRUE, updateEdges = TRUE, all.blockframes = TRUE)
                setSrcLabel(GW.top$env$viewLabel)
            }
            "zoomIn" <- function() {
                function(x, y) {
                  xy <- replaceXY(x, y, rep(50, local.N))
                  zoom(xy, sqrt(sqrt(2)))
                }
            }
            "zoomOut" <- function() {
                function(x, y) {
                  xy <- replaceXY(x, y, rep(50, local.N))
                  zoom(xy, 1/sqrt(sqrt(2)))
                }
            }
            "resizeCanvas" <- function() {
                function(x, y, ...) {
                }
            }
            "initFactorVariables" <- function(dg.factorVertexList) {
                if (length(dg.factorVertexList) > 0) {
                  itemsFactors <<- vector("list", length(dg.factorVertexList))
                  itemsFactorEdges <<- vector("list", length(dg.factorVertexList))
                  namesFactorVertices <<- Names(dg.factorVertexList)
                  positionsFactorVertices <<- Positions(dg.factorVertexList)
                  positionsFactorLabels <<- positionsFactorVertices
                  positionsFactorLabels[, 1] <<- positionsFactorLabels[, 
                    1] + 0.1 * control$w
                  factorLabels <<- Labels(dg.factorVertexList)
                  colorsFactorVertices <<- Colors(dg.factorVertexList)
                  blocksFactorVertices <<- rep(0, length(dg.factorVertexList))
                  strataFactorVertices <<- rep(0, length(dg.factorVertexList))
                  if (!is.matrix(positionsFactorVertices)) 
                    warning("Positions of factor-vertices should have same number of coordinates")
                  else if (!(dim(positionsFactorVertices)[2] == 
                    dim(positionsVertices)[2])) 
                    warning("Factor-vertices should have same number of coordinates as vertices")
                }
            }
            "initExtraVariables" <- function(dg.extraList) {
                if (length(dg.extraList) > 0) {
                  itemsExtras <<- vector("list", length(dg.extraList))
                  itemsExtraEdges <<- vector("list", length(dg.extraList))
                  namesExtraVertices <<- Names(dg.extraList)
                  positionsExtraVertices <<- Positions(dg.extraList)
                  positionsExtraLabels <<- positionsExtraVertices
                  positionsExtraLabels[, 1] <<- positionsExtraLabels[, 
                    1] + 0.1 * control$w
                  extraLabels <<- Labels(dg.extraList)
                  colorsExtraVertices <<- Colors(dg.extraList)
                  strataExtraVertices <<- rep(0, length(dg.extraList))
                  blocksExtraVertices <<- rep(0, length(dg.extraList))
                  if (!is.matrix(positionsExtraVertices)) 
                    warning("Positions of extra-vertices should have same number of coordinates")
                  else if (!(dim(positionsExtraVertices)[2] == 
                    dim(positionsVertices)[2])) 
                    warning("Extra-vertices should have same number of coordinates as vertices")
                }
            }
            "replaceBlocks" <- function(blockList) {
                if (!is.null(blockList)) {
                  positionsBlocks <- Positions(blockList)
                  d <- dim(positionsBlocks)
                  positionsBlocks <<- array(positionsBlocks, 
                    dim = c(d[1], d[2]/2, 2))
                  positionsClosedBlocks <- matrix(rep(NA, local.N * 
                    length(blockList)), ncol = local.N)
                  positionsClosedBlocks <<- apply(positionsBlocks, 
                    c(1, 2), mean)
                  blockReferences <<- 1:length(blockList)
                  positionsBlockLabels <<- matrix(rep(0, local.N * 
                    length(blockList)), ncol = local.N)
                  blockLabels <<- Labels(blockList)
                  strataBlocks <<- Strata(blockList)
                  setUpdateAll()
                }
            }
            "replaceVertices" <- function(vertexList) {
                if (length(vertexList) > 0) {
                  namesVertices <<- Names(vertexList)
                  positionsVertices <<- Positions(vertexList)
                  if (is.matrix(positionsVertices)) {
                    positionsLabels <<- positionsVertices
                    positionsLabels[, 1] <<- positionsLabels[, 
                      1] + 0.1 * control$w
                    Labels <<- Labels(vertexList)
                    colorsVertices <<- Colors(vertexList)
                    blocksVertices <<- Blockindices(vertexList)
                    strataVertices <<- Strata(vertexList)
                    setUpdateVertices()
                    setUpdatePositions()
                  }
                  if (!is.matrix(positionsVertices)) 
                    warning("Positions of extra-vertices should have same number of coordinates;")
                  else if (!(dim(positionsVertices)[2] == local.N)) 
                    warning("New vertices should have same number of coordinates as old vertices;")
                }
            }
            "drawFactors" <- function(X.FactorEdges, X.FactorVertices) {
                dg@factorEdgeList <<- X.FactorEdges
                GraphWindow@dg@factorEdgeList <<- dg@factorEdgeList
                dg@factorVertexList <<- X.FactorVertices
                GraphWindow@dg@factorVertexList <<- dg@factorVertexList
                initFactorVariables(dg@factorVertexList)
                for (i in seq(along = dg@factorEdgeList)) {
                  f <- dg@factorEdgeList[[i]]@vertex.indices[1]
                  t <- dg@factorEdgeList[[i]]@vertex.indices[2]
                  E <- append.index.edge(c(f, t), edge.type = "FactorEdge")
                  drawEdge(E[[length(E)]], length(E), lower = TRUE, 
                    edge.type = "FactorEdge")
                }
                if (length(dg@factorVertexList) > 0) 
                  for (i in seq(along = dg@factorVertexList)) drawVertex(-i, 
                    w = control$w, vertexcolor = control$vertexColor, 
                    vertex.type = "Factor")
            }
            "setMainMenu" <- function() {
                topMenu <- tkmenu(GW.top)
                tkconfigure(GW.top, menu = topMenu)
                fileMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(fileMenu, "command", label = "Save as Postscript ... ", 
                  command = function() {
                    fileName <- tclvalue(tkgetSaveFile(initialfile = "dynamicGraph.ps", 
                      filetypes = "{{Postsript Files} {.ps}}"))
                    if (!nchar(fileName)) 
                      tkmessageBox(message = "No file was selected!")
                    else tkpostscript(canvas, file = fileName, 
                      width = control$width, height = control$height)
                  })
                tkadd(fileMenu, "command", label = "Quit", 
                  command = function() destroyView(txt = "Quit")())
                tkadd(fileMenu, "command", label = "Popup selected in panel", 
                  command = function() {
                    popupSelectedInPanel()
                  })
                tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
                userMenu <- tkmenu(topMenu, tearoff = FALSE)
                "UserMainMenu" <- function(i) {
                  force(i)
                  function(...) {
                    sinkView(control$UserMenus[[i]], blocks = TRUE)
                    arguments <- Args()
                    control$UserMenus[[i]]$command(object, Arguments = arguments)
                  }
                }
                if (length(control$UserMenus) > 0) 
                  for (i in seq(along = control$UserMenus)) 
                    if (names(control$UserMenus[i]) == "MainUser") 
                    tkadd(userMenu, "command", label = control$UserMenus[[i]]$label, 
                      command = UserMainMenu(i))
                tkadd(topMenu, "cascade", label = "User", menu = userMenu)
                variableMenu <- tkmenu(topMenu, tearoff = FALSE)
                "selectVariableDialog" <- function() {
                  vertexlabels <- Labels[dg@visibleVertices]
                  ReturnVal <- selectDialog("Variable vertex selection entry", 
                    "Select variable", vertexlabels, top = GW.top)
                  if ((length(ReturnVal) > 0) && (ReturnVal != 
                    "ID_CANCEL")) {
                    variableChoice <- vertexlabels[ReturnVal]
                    j <- (1:length(Labels))[Labels == variableChoice]
                    subActivateVertex(j, color = "green", vertex.type = "Vertex")
                    msg <- paste("Click other vertex to add edge to ", 
                      variableChoice, sep = "")
                    tkmessageBox(message = msg)
                  }
                }
                tkadd(variableMenu, "command", label = "Highlight vertex (for adding edge)", 
                  command = function() selectVariableDialog())
                "selectOtherDialog" <- function(slave = TRUE) {
                  "not.in" <- function(x, l = max(x)) {
                    y <- rep(TRUE, l)
                    y[x] <- FALSE
                    return((1:l)[y])
                  }
                  vertexnames <- Names(vertexList)[not.in(dg@visibleVertices, 
                    length(vertexList))]
                  ReturnVal <- selectDialog("Variable vertex selection entry", 
                    "Select variable", vertexnames, top = GW.top)
                  if ((length(ReturnVal) > 0) && (ReturnVal != 
                    "ID_CANCEL")) {
                    variableChoice <- vertexnames[ReturnVal]
                    index <- nameToVertexIndex(variableChoice, 
                      vertexList)
                    if (length(vertexnames) > 0) 
                      subAddVertex(index, slave = slave)
                  }
                }
                tkadd(variableMenu, "command", 
                  label = "Select vertex among variables not displayed (slave)", 
                  command = function() selectOtherDialog())
                tkadd(variableMenu, "command", 
                  label = "Select vertex among variables not displayed (here)", 
                  command = function() selectOtherDialog(slave = FALSE))
              # tkadd(variableMenu, "command", label = paste("'addVertex', selected vertices and edges"), 
              #   command = function() subAddVertex(0, vertex.type = "selected", 
              #     slave = FALSE))
                tkadd(variableMenu, "command", 
                  label = paste("'dropVertex', selected vertices and edges"), 
                  command = function() subDropVertex(0, vertex.type = "selected", 
                    slave = FALSE))
                tkadd(variableMenu, "command", 
                  label = "Create new variable (not displayed before selected)", 
                  accelerator = "(~ <F3>)", command = function() subNewVertex(rep(0, 
                    local.N)))
                tkadd(topMenu, "cascade", label = "Variables", 
                  menu = variableMenu)
                expVarMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(variableMenu, "cascade", label = "Export & show", 
                  menu = expVarMenu)
                edgesMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(edgesMenu, "command", label = " ... ", 
                  command = function() {
                  })
                tkadd(edgesMenu, "command", label = "Delete all edge labels (for good)", 
                  accelerator = "<F6>", command = function() {
                    deleteAllEdgeLabels()()
                    control$namesOnEdges <<- FALSE
                  })
                tkadd(edgesMenu, "command", label = "Delete all edge labels (temporary)", 
                  command = function() {
                    deleteAllEdgeLabels(permanent = FALSE)()
                  })
                tkadd(edgesMenu, "command", label = paste("'addEdge', selected vertices and edges"), 
                  command = function() subAddEdge(0, 0, "none", 
                    "none", edge.type = "selected", slave = FALSE, 
                    edgeClass = NULL)) # 'edgeClass' = NULL ???
                tkadd(edgesMenu, "command", label = paste("'dropEdge', selected vertices and edges"), 
                  command = function() subSubDropEdge(0, 0, 0, 
                    "none", "none", edge.type = "selected", slave = FALSE))
                expEdgesMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(edgesMenu, "cascade", label = "Export & show", 
                  menu = expEdgesMenu)
                tkadd(topMenu, "cascade", label = "Edges", menu = edgesMenu)
                generatorsMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(generatorsMenu, "command", label = " ... e.i. factors ", 
                  command = function() {
                  })
                expGenMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(generatorsMenu, "cascade", label = "Export & show", 
                  menu = expGenMenu)
                tkadd(topMenu, "cascade", label = "Generators", 
                  menu = generatorsMenu)
                blocksMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(blocksMenu, "command", label = " ... ", 
                  command = function() {
                  })
                tkadd(blocksMenu, "command", label = "Add new block", 
                  accelerator = "(~ <F4>)", command = function() {
                    n <- length(blockList)
                    position <- rep(0, local.N) + n
                    delta <- c(10, 10, rep(50, local.N - 2))
                    position <- matrix(c(position - delta, position + 
                      delta), nrow = 2, byrow = TRUE)
                    new.Block(position)
                  })
                expBlockMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(blocksMenu, "cascade", label = "Export & show", 
                  menu = expBlockMenu)
                tkadd(topMenu, "cascade", label = "Blocks", menu = blocksMenu)
                graphMenu <- tkmenu(topMenu, tearoff = FALSE)
                slaveMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(slaveMenu, "command", label = "Same model", 
                  command = function() {
                    makeSlave(sameModel = TRUE, local.Views = dm.frameViews, 
                      Object = object, label = control$label)
                  })
                tkadd(slaveMenu, "command", label = "Same model, switch 'variableFrame'", 
                  command = function() {
                    makeSlave(sameModel = TRUE, local.Views = dm.frameViews, 
                      Object = object, label = control$label, variableFrame = !control$variableFrame)
                  })
                tkadd(slaveMenu, "command", label = "Copy model", 
                  command = function() {
                    makeSlave(sameModel = FALSE, local.Views = NULL, 
                      Object = object, label = "Default")
                  })
                tkadd(graphMenu, "cascade", label = "Make slave window", 
                  menu = slaveMenu)
                tkadd(graphMenu, "command", label = "Set 'viewType', class of graph window", 
                  command = function() {
                    Arguments <- Args()
                    ReturnVal <- selectDialog("Test selectDialog Entry", 
                      "Select name", control$viewClasses[, 1], 
                      top = GW.top)
                    if ((length(ReturnVal) > 0) && (ReturnVal != 
                      "ID_CANCEL")) {
                      setModel(object, txt = "SetClassOfView", 
                        setUpdate = FALSE, RR = NULL)
                      class(GraphWindow) <<- control$viewClasses[ReturnVal, 
                        2]
                      dg@viewType <<- control$viewClasses[ReturnVal, 
                        1]
                      tkconfigure(GW.top$env$viewLabel, text = dg@viewType)
                      updateModel()
                    }
                  })
                refreshMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(refreshMenu, "command", label = "Refresh view (set positions as 'stored')", 
                  command = function() {
                    subUpdateGraphWindow("Update from main menu", 
                      redrawVertices = TRUE, raiseEdges = TRUE, 
                      updateEdges = TRUE, all.blockframes = TRUE)
                    sinkView(NULL, edges = TRUE, blocks = TRUE)
                  })
                tkadd(refreshMenu, "command", label = "Redraw graph window (more refreshing)", 
                  command = function() {
                    Arguments <- Args()
                    redrawView(frameModels = dgm.frameModels, 
                      frameViews = dm.frameViews, graphWindow = GraphWindow, 
                      dg = dg, control = control, Arguments = Arguments)
                  })
                tkadd(graphMenu, "cascade", label = "Refresh", 
                  menu = refreshMenu)
                transformationMenu <- tkmenu(topMenu, tearoff = FALSE, 
                  postcommand = function() nullTrans <<- tclVar(is.null(transformation)))
                nullTrans <- tclVar(is.null(transformation))
                tkadd(transformationMenu, "radiobutton", label = "Reset (enable) transformation", 
                  variable = nullTrans, value = 0, command = function() setTransformation(diag(1, 
                    local.N)))
                tkadd(transformationMenu, "radiobutton", label = "Disable rotation", 
                  variable = nullTrans, value = 1, command = function() setTransformation(NULL))
                tkadd(graphMenu, "cascade", label = "Rotation: The transformation", 
                  menu = transformationMenu)
                tkadd(graphMenu, "command", label = "Close", 
                  accelerator = "(~ [X])", command = function() destroyView(txt = "Close")())
                if (control$permitZoom) {
                  tkadd(graphMenu, "command", label = "Zoom in", 
                    accelerator = "(~ <F1>)", command = function() zoom(zoomCenter, 
                      2^0.25))
                  tkadd(graphMenu, "command", label = "Zoom out", 
                    accelerator = "(~ <F2>)", command = function() zoom(zoomCenter, 
                      2^-0.25))
                }
                nameLabelMenu <- tkmenu(topMenu, tearoff = FALSE, 
                  postcommand = function() useNames <<- tclVar(control$useNamesForLabels))
                useNames <- tclVar(control$useNamesForLabels)
                tkadd(nameLabelMenu, "radiobutton", label = "Labels is labels", 
                  variable = useNames, value = 0, command = function() {
                    control$useNamesForLabels <<- FALSE
                    subUpdateGraphWindow("useNamesForLabels", 
                      redrawVertices = TRUE, raiseEdges = FALSE, 
                      updateEdges = FALSE, all.blockframes = FALSE)
                  })
                tkadd(nameLabelMenu, "radiobutton", label = "Use names for labels", 
                  variable = useNames, value = 1, command = function() {
                    control$useNamesForLabels <<- TRUE
                    subUpdateGraphWindow("useNamesForLabels", 
                      redrawVertices = TRUE, raiseEdges = FALSE, 
                      updateEdges = FALSE, all.blockframes = FALSE)
                  })
                tkadd(graphMenu, "cascade", label = "Names for labels", 
                  menu = nameLabelMenu)
                tkadd(topMenu, "cascade", label = "Graph", menu = graphMenu)
                exportMenu <- tkmenu(topMenu, tearoff = FALSE)
                expGraphMenu <- tkmenu(topMenu, tearoff = FALSE)
                tkadd(graphMenu, "cascade", label = "Export & show", 
                  menu = expGraphMenu)
                "argsExport" <- function() {
                  ReturnVal <- modalDialog("Args Name Entry", 
                    "Enter name for Args", "Args", top = GW.top)
                  sinkView(NULL, edges = TRUE, blocks = TRUE)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, Args(), pos = 1)
                }
                tkadd(expGraphMenu, "command", label = "Assign 'Args' in .GlobalEnv", 
                  command = function() argsExport())
                tkadd(exportMenu, "command", label = "Assign 'Args' in .GlobalEnv", 
                  command = function() argsExport())
                "updateAllPositionsMenu" <- function() {
                  subSinkAllFrames("Positions", txt = "Update")
                  Str(dgm.frameModels, excludeSlots = TRUE)
                }
                tkadd(expGraphMenu, "command", label = "Update positions of all views of 'frameModels'", 
                  command = function() updateAllPositionsMenu())
                tkadd(exportMenu, "command", label = "Update positions of all views of 'frameModels'", 
                  command = function() updateAllPositionsMenu())
                "updateAllArgumentsMenu" <- function() {
                  subSinkAllFrames("Arguments", menuItem = NULL, 
                    vertices = TRUE, edges = TRUE, blocks = TRUE)
                  Str(dgm.frameModels, excludeSlots = TRUE)
                }
                tkadd(expGraphMenu, "command", label = "Update all slots of 'frameModels'", 
                  command = function() updateAllArgumentsMenu())
                tkadd(exportMenu, "command", label = "Update all slots of 'frameModels'", 
                  command = function() updateAllArgumentsMenu())
                "latticeModelsExport" <- function() {
                  ReturnVal <- modalDialog("Lattice for Models Name Entry", 
                    "Enter name for lattice for models", "frameModels", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  sinkView(NULL, edges = TRUE, blocks = TRUE)
                  Str(dgm.frameModels, excludeSlots = TRUE)
                  message("Consider: \"Update all slots of 'frameModels'\" before 'export'")
                  assign(ReturnVal, dgm.frameModels, pos = 1)
                }
                tkadd(expGraphMenu, "command", label = "Assign 'frameModels' in .GlobalEnv", 
                  command = function() latticeModelsExport())
                tkadd(exportMenu, "command", label = "Assign 'frameModels' in .GlobalEnv", 
                  command = function() latticeModelsExport())
                "latticeGraphsExport" <- function() {
                  ReturnVal <- modalDialog("Lattice for Views Name Entry", 
                    "Enter name for lattice for views", "frameViews", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  sinkView(NULL, edges = TRUE, blocks = TRUE)
                  Str(dm.frameViews, excludeSlots = c("class", 
                    "label", "label.position"))
                  message("Consider: \"Update all slots of 'frameModels'\" ")
                  assign(ReturnVal, dm.frameViews, pos = 1)
                }
                tkadd(expGraphMenu, "command", label = "Assign 'frameViews' in .GlobalEnv", 
                  command = function() latticeGraphsExport())
                tkadd(exportMenu, "command", label = "Assign 'frameViews' in .GlobalEnv", 
                  command = function() latticeGraphsExport())
                "graphWindowExport" <- function() {
                  ReturnVal <- modalDialog("GraphWindow Name Entry", 
                    "Enter name for graphWindow", "graphWindow", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  sinkView(NULL, edges = TRUE, blocks = TRUE)
                  Str(GraphWindow, excludeSlots = c("class", 
                    "label", "label.position"))
                  assign(ReturnVal, GraphWindow, pos = 1)
                }
                tkadd(expGraphMenu, "command", label = "Assign 'graphWindow' in .GlobalEnv", 
                  command = function() graphWindowExport())
                tkadd(exportMenu, "command", label = "Assign 'graphWindow' in .GlobalEnv", 
                  command = function() graphWindowExport())
                "objectExport" <- function() {
                  ReturnVal <- modalDialog("Object Name Entry", 
                    "Enter name for object", "object", top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, object, pos = 1)
                }
                tkadd(exportMenu, "command", label = "Assign 'object' in .GlobalEnv", 
                  command = function() objectExport())
                tkadd(expGraphMenu, "command", label = "Assign 'object' in .GlobalEnv", 
                  command = function() objectExport())
                "transformationExport" <- function() {
                  ReturnVal <- modalDialog("Transformation Name Entry", 
                    "Enter name for transformation", "transformation", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, transformation, pos = 1)
                }
                tkadd(exportMenu, "command", label = "Assign 'transformation' in .GlobalEnv", 
                  command = function() transformationExport())
                tkadd(expGraphMenu, "command", label = "Assign 'transformation' in .GlobalEnv", 
                  command = function() transformationExport())
                "topExport" <- function() {
                  ReturnVal <- modalDialog("Top Name Entry", 
                    "Enter name for top", "top", top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, GW.top, pos = 1)
                }
                tkadd(exportMenu, "command", label = "Assign 'top' in .GlobalEnv", 
                  command = function() topExport())
                tkadd(expGraphMenu, "command", label = "Assign 'top' in .GlobalEnv", 
                  command = function() topExport())
                "canvasExport" <- function() {
                  ReturnVal <- modalDialog("Canvas Name Entry", 
                    "Enter name for canvas", "canvas", top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, GW.top$env$canvas, pos = 1)
                }
                tkadd(exportMenu, "command", label = "Assign 'canvas' in .GlobalEnv", 
                  command = function() canvasExport())
                tkadd(expGraphMenu, "command", label = "Assign 'canvas' in .GlobalEnv", 
                  command = function() canvasExport())
                "viewLabelExport" <- function() {
                  ReturnVal <- modalDialog("ViewLabel Name Entry", 
                    "Enter name for viewLabel", "viewLabel", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  assign(ReturnVal, GW.top$env$viewLabel, pos = 1)
                }
                tkadd(exportMenu, "command", label = "Assign 'viewLabel' in .GlobalEnv", 
                  command = function() viewLabelExport())
                tkadd(expGraphMenu, "command", label = "Assign 'viewLabel' in .GlobalEnv", 
                  command = function() viewLabelExport())
                tkadd(expVarMenu, "command", label = "'print(selectedNodesMatrix(selectedNodes()))'", 
                  command = function() print(selectedNodesMatrix()))
                "vertexListExport" <- function() {
                  ReturnVal <- modalDialog("Vertices Name Entry", 
                    "Enter name for vertices", "vertexList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  vertices <- sinkVertexList()
                  assign(ReturnVal, vertexList, pos = 1)
                }
                tkadd(expVarMenu, "command", label = "Assign 'vertexList' in .GlobalEnv", 
                  command = function() vertexListExport())
                tkadd(exportMenu, "command", label = "Assign 'vertexList' in .GlobalEnv", 
                  command = function() vertexListExport())
                tkadd(expVarMenu, "command", label = "'print(asDataFrame(vertexList))'", 
                  command = function() print(asDataFrame(vertexList)))
                tkadd(expEdgesMenu, "command", label = "'print(selectedEdgesMatrix(selectedEdges()))'", 
                  command = function() print(selectedEdgesMatrix()))
                "edgeListExport" <- function() {
                  ReturnVal <- modalDialog("Edges Name Entry", 
                    "Enter name for Edges", "edgeList", top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  edge.List <- sinkEdgeList()
                  assign(ReturnVal, edge.List, pos = 1)
                }
                tkadd(expEdgesMenu, "command", label = "Assign 'edgeList' in .GlobalEnv", 
                  command = function() edgeListExport())
                tkadd(exportMenu, "command", label = "Assign 'edgeList' in .GlobalEnv", 
                  command = function() edgeListExport())
                tkadd(expEdgesMenu, "command", label = "'print(asDataFrame(dg@edgeList))'", 
                  command = function() print(asDataFrame(dg@edgeList)))
                "extraListExport" <- function() {
                  ReturnVal <- modalDialog("ExtraVertices Name Entry", 
                    "Enter name for the extravertices", "extraList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  extra.List <- sinkExtraVertexList()
                  assign(ReturnVal, extra.List, pos = 1)
                }
                tkadd(expEdgesMenu, "command", label = "Assign 'extraList' in .GlobalEnv", 
                  command = function() extraListExport())
                tkadd(exportMenu, "command", label = "Assign 'extraList' in .GlobalEnv", 
                  command = function() extraListExport())
                tkadd(expEdgesMenu, "command", label = "'print(asDataFrame(dg@extraList))'", 
                  command = function() print(asDataFrame(dg@extraList)))
                "extraEdgeListExport" <- function() {
                  ReturnVal <- modalDialog("ExtraEdges Name Entry", 
                    "Enter name for the extraedges", "extraEdgeList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  extra.Edge.List <- sinkExtraEdgeList()
                  assign(ReturnVal, extra.Edge.List, pos = 1)
                }
                tkadd(expEdgesMenu, "command", label = "Assign 'extraEdgeList' in .GlobalEnv", 
                  command = function() extraEdgeListExport())
                tkadd(exportMenu, "command", label = "Assign 'extraEdgeList' in .GlobalEnv", 
                  command = function() extraEdgeListExport())
                tkadd(expEdgesMenu, "command", label = "'print(asDataFrame(dg@extraEdgeList))'", 
                  command = function() print(asDataFrame(dg@extraEdgeList)))
                "factorVertexListExport" <- function() {
                  ReturnVal <- modalDialog("FactorVertices Name Entry", 
                    "Enter name for the factorvertices", "factorVertexList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  factorVertex.List <- sinkFactorVertexList()
                  assign(ReturnVal, factorVertex.List, pos = 1)
                }
                tkadd(expGenMenu, "command", label = "Assign 'factorVertexList' in .GlobalEnv", 
                  command = function() factorVertexListExport())
                tkadd(exportMenu, "command", label = "Assign 'factorVertexList' in .GlobalEnv", 
                  command = function() factorVertexListExport())
                tkadd(expGenMenu, "command", label = "'print(asDataFrame(dg@factorVertexList))'", 
                  command = function() print(asDataFrame(dg@factorVertexList)))
                "factorEdgeListExport" <- function() {
                  ReturnVal <- modalDialog("FactorEdges Name Entry", 
                    "Enter name for the factoredges", "factorEdgeList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  factor.Edge.List <- sinkFactorEdgeList()
                  assign(ReturnVal, factor.Edge.List, pos = 1)
                }
                tkadd(expGenMenu, "command", label = "Assign 'factorEdgeList' in .GlobalEnv", 
                  command = function() factorEdgeListExport())
                tkadd(exportMenu, "command", label = "Assign 'factorEdgeList' in .GlobalEnv", 
                  command = function() factorEdgeListExport())
                tkadd(expGenMenu, "command", label = "'print(asDataFrame(dg@factorEdgeList))'", 
                  command = function() print(asDataFrame(dg@factorEdgeList)))
                "blockListExport" <- function() {
                  ReturnVal <- modalDialog("Blocklist Name Entry", 
                    "Enter name for the blocklist", "blockList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  block.List <- sinkBlockList()
                  assign(ReturnVal, block.List, pos = 1)
                }
                tkadd(expBlockMenu, "command", label = "Assign 'blockList' in .GlobalEnv", 
                  command = function() blockListExport())
                tkadd(exportMenu, "command", label = "Assign 'blockList' in .GlobalEnv", 
                  command = function() blockListExport())
                tkadd(expBlockMenu, "command", label = "'print(asDataFrame(blockList))'", 
                  command = function() print(asDataFrame(blockList)))
                if (!is.null(blockTree)) {
                  "blockTreeExport" <- function() {
                    ReturnVal <- modalDialog("Blocktree Name Entry", 
                      "Enter name for the blocktree", "blockTree", 
                      top = GW.top)
                    if (ReturnVal == "ID_CANCEL") 
                      return()
                    sinkBlockList()
                    assign(ReturnVal, blockTree, pos = 1)
                  }
                  tkadd(expBlockMenu, "command", label = "Assign 'blockTree' in .GlobalEnv", 
                    command = function() blockTreeExport())
                  tkadd(exportMenu, "command", label = "Assign 'blockTree' in .GlobalEnv", 
                    command = function() blockTreeExport())
                }
                "blockEdgeListExport" <- function() {
                  ReturnVal <- modalDialog("BlockEdges Name Entry", 
                    "Enter name for the blockedges", "blockEdgeList", 
                    top = GW.top)
                  if (ReturnVal == "ID_CANCEL") 
                    return()
                  blockEdge.List <- sinkBlockEdges()
                  assign(ReturnVal, blockEdge.List, pos = 1)
                }
                tkadd(expBlockMenu, "command", label = "Assign 'blockEdgeList' in .GlobalEnv", 
                  command = function() blockEdgeListExport())
                tkadd(exportMenu, "command", label = "Assign 'blockEdgeList' in .GlobalEnv", 
                  command = function() blockEdgeListExport())
                tkadd(expBlockMenu, "command", label = "'print(asDataFrame(dg@blockEdgeList))'", 
                  command = function() print(asDataFrame(dg@blockEdgeList)))
                tkadd(topMenu, "cascade", label = "Export", menu = exportMenu)
            }
            top <- NULL
            box <- NULL
            canvas <- NULL
            viewLabel <- NULL
            tags <- NULL
            updateCountVertices <- updateCountVerticesMain
            updateCountPositions <- updateCountPositionsMain
            updateCountBlocks <- updateCountBlocksMain
            updateCountBlockEdges <- updateCountBlockEdgesMain
            if (setUpdateCountModelMain) 
                updateCountModelMain <<- updateCountModelMain + 
                  1
            updateCountModel <- updateCountModelMain
            updateWindow <- TRUE
            args <- list(...)
            "extract" <- function(x, y = x, z = "$", a = "graphComponents", 
                b = paste(a, z, y, collapse = "")) {
                eval(parse(text = paste(c("if (is.null(", x, 
                  ") && !is.null(", b, ")) {", x, " <<- ", b, 
                  " } "), collapse = "")))
            }
            Arguments <- args$Arguments
            if (!is.null(Arguments)) {
                X <- matrix(c("frameModels", "frameModels", "frameViews", 
                  "frameViews", "dg", "dg", "control", "control", 
                  "top", "top", "box", "box", "canvas", "canvas", 
                  "viewLabel", "viewLabel", "tags", "tags"), 
                  ncol = 2, byrow = TRUE)
                apply(X, 1, function(i) extract(i[1], i[2], a = "Arguments"))
                if ((is.null(dg@visibleBlocks)) && (!is.null(Arguments$visibleBlocks) || 
                  (length(Arguments$visibleBlocks) == 0))) 
                  dg@visibleBlocks <- Arguments$visibleBlocks
            }
            if (is.null(dg)) 
                dg <- .newDgGraphEdges(vertexList = vertexList, 
                  blockList = blockList, ...)
            if (is.null(dg@oriented)) 
                dg@oriented <- FALSE
            if (is.null(control)) 
                control <- dg.control(...)
            if (!is.null(args$vertexList)) {
                vertexList <<- args$vertexList
                replaceVertices(vertexList)
            }
            if (!is.null(args$blockList)) {
                blockList <<- args$blockList
                replaceBlocks(blockList)
            }
            if (is.null(args$blockList) && !is.null(args$blockTree)) {
                blockList <<- blockTreeToList(args$blockTree)
                replaceBlocks(blockList)
            }
            openTreeBlock <- rep(TRUE, length(blockList))
            if (is.null(dg@visibleVertices)) 
                dg@visibleVertices <- 1:length(vertexList)
            if (is.null(dg@visibleBlocks)) 
                dg@visibleBlocks <- 1:length(blockList)
            zoomPositions <- NULL
            Scale <- 1
            zoomCenterSet <- FALSE
            zoomCenter <- NULL
            upperLeft <- c(min.x, min.y, rep(0, local.N - 2))
            lowerRight <- c(max.x, max.y, rep(100, local.N - 
                2))
            if (is.null(dg@extraList)) 
                dg@extraList <- new("dg.VertexList")
            if (is.null(dg@edgeList)) 
                dg@edgeList <- new("dg.VertexEdgeList")
            if (is.null(dg@factorEdgeList)) 
                dg@factorEdgeList <- new("dg.FactorEdgeList")
            if (is.null(dg@extraEdgeList)) 
                dg@extraEdgeList <- new("dg.ExtraEdgeList")
            if (is.null(dg@blockEdgeList)) 
                dg@blockEdgeList <- new("dg.BlockEdgeList")
            if (!is.na(control$namesOnEdges) && !control$namesOnEdges) {
                for (i in seq(along = dg@edgeList)) label(dg@edgeList[[i]]) <- ""
                for (i in seq(along = dg@factorEdgeList)) label(dg@factorEdgeList[[i]]) <- ""
                for (i in seq(along = dg@extraEdgeList)) label(dg@extraEdgeList[[i]]) <- ""
                for (i in seq(along = dg@blockEdgeList)) label(dg@blockEdgeList[[i]]) <- ""
            }
            Oriented <- FALSE
            if (!is.na(dg@oriented)) 
                Oriented <- dg@oriented
            if (!Oriented) {
                for (i in seq(along = dg@edgeList)) {
                  if (!is.na(dg@edgeList[[i]]@oriented)) 
                    if (dg@edgeList[[i]]@oriented) {
                      dg@oriented <- TRUE
                      Oriented <- TRUE
                    }
                }
                if (Oriented) 
                  message("Oriented edge found!")
            }
            if (redraw) 
                graphWindow <- NULL
            GW.top <- (top)
            GW.tags <- tags
            GW.env <- NULL
            if (is.numeric(frameViews)) 
                dm.frameViews <- dgm.frameModels@models[[frameViews]]
            else dm.frameViews <- frameViews
            if (is.null(graphWindow)) {
                ArgWindow <- FALSE
                m <- dm.frameViews@index
                n <- 1
                if (!.IsEmpty(dgm.frameModels@models[[m]]@graphs)) 
                  n <- length(dgm.frameModels@models[[m]]@graphs) + 
                    1
                fv.env <- .get.env.frameViews(frameViews = dm.frameViews, 
                  frameModels = dgm.frameModels)
                GraphWindow <- newGraph(dg@viewType, dg, background = control$background, 
                  title = ifelse(control$label != "Default", 
                    control$label, paste(n, m, sep = "@")), parent = fv.env, 
                  index = n, width = control$width, height = control$height)
                if (control$saveTkReferences) {
                  assign("top", GW.top, envir = GW.env$env)
                  assign("tags", GW.tags, envir = GW.env$env)
                }
            }
            else {
                ArgWindow <- TRUE
                if (is.numeric(graphWindow)) 
                  GraphWindow <- dm.frameViews@graphs[[graphWindow]]
                else GraphWindow <- graphWindow
                m <- dm.frameViews@index
                gw.env <- .get.env.graphWindow(graphWindow = GraphWindow, 
                  frameViews = dgm.frameModels@models[[m]], frameModels = dgm.frameModels)
                GW.env <- gw.env
                if (is.element("top", ls(GW.env$env))) {
                  GW.top <- evalq(top, GW.env$env)
                  GW.tags <- evalq(tags, GW.env$env)
                }
                m <- dm.frameViews@index
                n <- graphWindow@index
                GraphWindow@id <- GraphWindow@id + 1
                if (!is.null(GW.top)) {
                  GraphWindow@dg@visibleVertices <- .nullToEmpty(dg@visibleVertices)
                  GraphWindow@dg@visibleBlocks <- .nullToEmpty(dg@visibleBlocks)
                  GraphWindow@dg@extraList <- .nullToList(dg@extraList, 
                    type = "dg.VertexList")
                  GraphWindow@dg@edgeList <- dg@edgeList
                  GraphWindow@dg@factorVertexList <- .nullToList(dg@factorVertexList, 
                    type = "dg.FactorVertexList")
                  GraphWindow@dg@factorEdgeList <- dg@factorEdgeList
                  GraphWindow@dg@extraEdgeList <- dg@extraEdgeList
                  GraphWindow@dg@blockEdgeList <- dg@blockEdgeList
                  deleteTags()
                  GW.tags <- list(NULL)
                  if (control$saveTkReferences) 
                    assign("tags", GW.tags, envir = GW.env$env)
                  bindBox(GW.top$env$box, label = "redrawView")
                }
                else {
                  GraphWindow <- newGraph(dg@viewType, dg, background = control$background, 
                    title = ifelse(control$label != "Default", 
                      control$label, paste(n, m, sep = "@")), 
                    id = GraphWindow@id, index = n, parent = dm.frameViews.env, 
                    width = control$width, height = control$height)
                }
            }
            activatedNode <- list(number = 0, vertex.type = "Null")
            selectedNodes <- list()
            activatedEdge <- list(number = 0, edge.type = NULL)
            selectedEdges <- list()
            itemsFactors <- NULL
            itemsFactorEdges <- NULL
            namesFactorVertices <- NULL
            itemsExtraEdges <- NULL
            positionsFactorVertices <- NULL
            positionsFactorLabels <- NULL
            factorLabels <- NULL
            colorsFactorVertices <- NULL
            strataFactorVertices <- NULL
            blocksFactorVertices <- NULL
            if (FALSE && !is.null(dg@factorVertexList)) {
                namesFactorVertices <- Names(dg@factorVertexList)
                positionsFactorVertices <- Positions(dg@factorVertexList)
                positionsFactorLabels <- NULL
                factorLabels <- Labels(dg@factorVertexList)
                colorsFactorVertices <- Colors(dg@factorVertexList)
                strataFactorVertices <- Strata(dg@factorVertexList)
                blocksFactorVertices <- Blockindices(dg@factorVertexList)
            }
            itemsExtras <- NULL
            itemsExtraEdges <- NULL
            namesExtraVertices <- NULL
            positionsExtraVertices <- NULL
            positionsExtraLabels <- NULL
            extraLabels <- NULL
            colorsExtraVertices <- NULL
            strataExtraVertices <- NULL
            blocksExtraVertices <- NULL
            if (FALSE && !is.null(dg@extraList)) {
                namesExtraVertices <- Names(dg@extraList)
                positionsExtraVertices <- Positions(dg@extraList)
                positionsExtraLabels <- NULL
                extraLabels <- Labels(dg@extraList)
                colorsExtraVertices <- Colors(dg@extraList)
                strataExtraVertices <- Strata(dg@extraList)
                blocksExtraVertices <- Blockindices(dg@extraList)
            }
            Angle <- 10
            canvas <- GW.top$env$canvas
            itemsEdges <- vector("list", length(vertexList))
            itemsVertices <- vector("list", length(vertexList))
            itemsOpenBlocks <- vector("list", length(blockList))
            itemsClosedBlocks <- vector("list", length(blockList))
            positionsEdgeLabels <- NULL
            closedVertex <- rep(FALSE, length(vertexList))
            closedBlock <- rep(FALSE, length(blockList))
            hiddenBlock <- rep(FALSE, length(blockList))
            if (!.IsEmpty(blockList)) {
                blockList <<- checkBlockList(blockList)
                if (!is.null(Arguments) && !is.null(Arguments$closedBlock)) 
                  closedBlock <- Arguments$closedBlock
                else closedBlock <- Closed(blockList)
                hiddenBlock <- closedBlock
                for (i in seq(along = blockList)) hiddenBlock[i] <- isInClosedBlock(i)
                if (!is.null(Arguments) && !is.null(Arguments$hiddenBlock)) {
                  hiddenBlock <- hiddenBlock | Arguments$hiddenBlock
                }
                dg@visibleBlocks <- unique(sort(dg@visibleBlocks))
                dg@visibleBlocks <- dg@visibleBlocks[dg@visibleBlocks != 
                  0]
                for (i in seq(along = vertexList)) {
                  s <- blocksVertices[i]
                  if (s %in% dg@visibleBlocks) 
                    closedVertex[i] <- closedBlock[s] || hiddenBlock[s]
                }
                if (!is.null(dg@visibleBlocks) && length(dg@visibleBlocks) == 
                  0) 
                  hiddenBlock <- hiddenBlock
            }
            itemsBlockEdges <- list(NULL)
            if (!.IsEmpty(blockList)) {
                itemsBlockEdges <- vector("list", length(blockList))
            }
            initFactorVariables(dg@factorVertexList)
            initExtraVariables(dg@extraList)
            if (!.IsEmpty(blockList)) 
                for (i in seq(along = blockList)) if (is.na(hiddenBlock[i])) 
                  message("hiddenBlock missing!!")
                else if (!hiddenBlock[i]) 
                  if (is.element(i, dg@visibleBlocks)) 
                    if (is.na(closedBlock[i])) 
                      message("closedBlock missing!!")
                    else if (closedBlock[i]) 
                      drawVertex(i, w = 10, vertexcolor = "Black", 
                        vertex.type = "ClosedBlock")
                    else drawBlock(blockList[[i]], i)
            for (i in seq(along = dg@edgeList)) drawEdge(dg@edgeList[[i]], 
                i, edge.type = "VertexEdge")
            for (i in seq(along = dg@factorEdgeList)) drawEdge(dg@factorEdgeList[[i]], 
                i, edge.type = "FactorEdge")
            for (i in seq(along = dg@extraEdgeList)) drawEdge(dg@extraEdgeList[[i]], 
                i, edge.type = "ExtraEdge")
            for (i in seq(along = dg@blockEdgeList)) drawEdge(dg@blockEdgeList[[i]], 
                i, edge.type = "BlockEdge")
            for (i in seq(along = vertexList)) if (is.element(i, 
                dg@visibleVertices)) 
                if (!closedVertex[i]) 
                  drawVertex(i, w = control$w, vertexcolor = control$vertexColor, 
                    vertex.type = "Vertex")
            if (length(dg@factorVertexList) > 0) 
                for (i in seq(along = dg@factorVertexList)) drawVertex(-i, 
                  w = control$w, vertexcolor = control$vertexColor, 
                  vertex.type = "Factor")
            if (length(dg@extraList) > 0) 
                for (i in seq(along = dg@extraList)) drawVertex(i, 
                  w = control$w, vertexcolor = control$vertexColor, 
                  vertex.type = "Extra")
            if (initialWindow) 
                update.edge.labels()
            setMainMenu()
            # tkbind(canvas, "<B1-Motion>", moveCanvas())
            tkbind(canvas, "<B2-Motion>", moveCanvas())
            tkbind(canvas, "<B3-Motion>", doHandRotate())
            tkbind(canvas, "<F12>", rockPlot(k = 1))
            tkbind(canvas, "<F1>", zoomIn())
            tkbind(canvas, "<F2>", zoomOut())
            tkbind(canvas, "<F3>", newVertexXY())
            tkbind(canvas, "<F4>", newBlockXY())
            tkbind(canvas, "<F5>", addLastEdge())
            tkbind(canvas, "<F6>", deleteAllEdgeLabels())
            tkbind(canvas, "<F7>", function(...) {
                print("<<F7>>")
            })
            tkbind(canvas, "<F12>", function(...) {
                print("<<F12>>")
            })
            tkbind(canvas, "<Up>", keyRotate(v = 1, sign = 1))
            tkbind(canvas, "<Down>", keyRotate(v = 1, sign = -1))
            tkbind(canvas, "<Left>", keyRotate(v = 2, sign = 1))
            tkbind(canvas, "<Right>", keyRotate(v = 2, sign = -1))
            tkbind(canvas, "<Prior>", keyRotate(v = 3, sign = 1))
            tkbind(canvas, "<Next>", keyRotate(v = 3, sign = -1))
            tkbind(canvas, "<Destroy>", destroyView(txt = "Canvas"))
            tkbind(canvas, "<Home>", function(...) {
                print("<<Fn-PgUp>>")
            })
            tkbind(canvas, "<End>", function(...) {
                print("<<Fn-PgDn>>")
            })
            tkbind(canvas, "<Alt_L>", function(...) {
                print("<<A>>")
            })
            tkbind(canvas, "<Delete>", function(...) {
                print("<<Delete>>")
            })
            tkbind(canvas, "<Pause>", function(...) {
                print("<<Pause>>")
            })
            tkbind(canvas, "<space>", function(...) {
                print("<<space>>")
            })
            tkbind(canvas, "<Tab>", function(...) {
                print("<<Tab>>")
            })
            tkbind(canvas, "<slash>", function(...) {
                print("<<slash>>")
            })
            tkbind(canvas, "<less>", function(...) {
                print("<<less>>")
            })
            tkbind(canvas, "<greater>", function(...) {
                print("<<greater>>")
            })
            tkbind(canvas, "<Escape>", function(...) {
                print("<<Escape>>")
            })
            tkbind(canvas, "<Shift-Up>", function(...) {
                print("<<Shift-Up>>")
            })
            tkbind(canvas, "<Shift-Down>", function(...) {
                print("<<Shift-Down>>")
            })
            tkbind(canvas, "<minus>", function(...) {
                print("<<minus>>")
            })
            tkbind(canvas, "<plus>", function(...) {
                print("<<plus>>")
            })
            tkbind(canvas, "<backslash>", function(...) {
                print("--backslash--")
            })
            tkbind(canvas, "<Alt-1>", function(...) {
                print("--A1#--")
            })
            if (control$debug.update) 
                tkbind(GW.top$env$viewLabel, "<Enter>", 
                       function() print(dg@viewType))
            tkbind(canvas, "<Configure>", resizeCanvas())
            if (control$enterLeaveUpdate) {
                tkbind(canvas, "<Enter>", updatePositions("Enter"))
                if (control$updateAllViews) 
                  tkbind(canvas, "<Leave>", sinkAllFrames("Positions", 
                    "Leave"))
                else tkbind(canvas, "<Leave>", updatePositions("Leave"))
            }
            m <- dm.frameViews@index
            if (control$saveFunctions) {
                assign("Update", update, envir = GW.env$env)
            }
            if (!redraw) 
                if (!ArgWindow) {
                  if (.IsEmpty(dgm.frameModels@models[[m]]@graphs)) {
                    dm.frameViews@graphs <<- list(GraphWindow)
                    dgm.frameModels@models[[m]] <<- dm.frameViews
                  }
                  else {
                    dm.frameViews@graphs <<- append(dgm.frameModels@models[[m]]@graphs, 
                      list(GraphWindow))
                    dgm.frameModels@models[[m]] <<- dm.frameViews
                  }
                }
                else {
                  dm.frameViews@graphs[[n]] <<- GraphWindow
                  dgm.frameModels@models[[m]]@graphs[[n]] <<- GraphWindow
                }
            if (!initialWindow) 
                updateBlockEdges()
            if (hasMethod("setGraphEdges", class(object))) {
                # message("Using 'setGraphEdges' for your model class.")
                object <<- setGraphEdges(object, dg = dg, ...)
            }
            else if (hasMethod("setGraphComponents", class(object))) {
                message("Please implement 'setGraphEdges' for your model class.")
                if (is.null(dg@edgeList)) 
                  dg@edgeList <- new("dg.VertexEdgeList")
                object <<- setGraphComponents(object, viewType = dg@viewType, 
                  visibleVertices = dg@visibleVertices, visibleBlocks = dg@visibleBlocks, 
                  extraVertices = dg@extraList, vertexEdges = dg@edgeList, 
                  blockEdges = dg@blockEdgeList, factorVertices = dg@factorVertexList, 
                  factorEdges = dg@factorEdgeList, extraEdges = dg@extraEdgeList)
            }
            if (control$returnNull) 
                return(NULL)
            else if (returnFrameModel) 
                invisible(dgm.frameModels)
            else invisible(GraphWindow)
        }
        updateCountModelMain <- 0
        if (is.numeric(frameViews)) 
            dm.frameViews <- dgm.frameModels@models[[frameViews]]
        else dm.frameViews <- frameViews
        if (is.null(dm.frameViews)) {
            ArgFrameViews <- FALSE
            m <- 1
            if (!.IsEmpty(dgm.frameModels@models)) 
                m <- length(dgm.frameModels@models) + 1
            Label <- control$label
            if (is.object(object) && ("name" %in% slotNames(object))) 
                Label <- paste(control$label, object@name)
            dm.frameViews.env <- .Dg.toplevel(parent = frameModelsEnv)
            dm.frameViews <- .newDynamicGraphModelObject(object, 
                label = Label, index = m, parent = frameModelsEnv, 
                env = dm.frameViews.env)
            if (control$saveFunctions) {
                assign("RedrawView", redrawView, envir = dm.frameViews.env$env)
            }
            else {
                txt <- "(( No 'update' and 'overwrite' since 'saveFunctions = FALSE'))"
                message(txt)
            }
        }
        else {
            ArgFrameViews <- TRUE
            m <- dm.frameViews@index
            dm.frameViews.env <- .get.env.frameViews(frameViews = dgm.frameModels@models[[m]], 
                frameModels = dgm.frameModels)
            if (is.null(dm.frameViews.env)) {
                dm.frameViews.env <- .Dg.toplevel(parent = frameModelsEnv)
                dm.frameViews@id.env <- dm.frameViews.env$ID
            }
            object <- dm.frameViews@model[[1]]
        }
        if (returnNewMaster || redraw) {
            graphs <- dm.frameViews@graphs
            if (!redraw) 
                dm.frameViews@graphs <- list()
            if (length(graphs) > 0) 
                for (i in (1:length(graphs))) if (!is.null(graphs[[i]])) {
                  ViewType <- control$viewClasses[control$viewClasses[, 
                    2] == class(graphs[[i]]), 1]
                  ldg <- .newDgGraphEdges(vertexList = vertexList, 
                    oriented = graphs[[i]]@dg@oriented, viewType = ViewType, 
                    visibleVertices = graphs[[i]]@dg@visibleVertices, 
                    visibleBlocks = graphs[[i]]@dg@visibleBlocks, 
                    edgeList = graphs[[i]]@dg@edgeList, blockList = blockList, 
                    blockEdgeList = graphs[[i]]@dg@blockEdgeList, 
                    factorVertexList = graphs[[i]]@dg@factorVertexList, 
                    factorEdgeList = graphs[[i]]@dg@factorEdgeList, 
                    extraList = graphs[[i]]@dg@extraList, extraEdgeList = graphs[[i]]@dg@extraEdgeList)
                  R <- redrawView(frameModels = dgm.frameModels, 
                    frameViews = dm.frameViews, graphWindow = graphs[[i]], 
                    dg = ldg, returnNewMaster = returnNewMaster, 
                    redraw = redraw, returnFrameModel = FALSE, 
                    control = control)
                }
        }
        if (!redraw) 
            if (!ArgFrameViews) {
                if (.IsEmpty(dgm.frameModels@models)) 
                  dgm.frameModels@models <<- list(dm.frameViews)
                else dgm.frameModels@models <<- append(dgm.frameModels@models, 
                  list(dm.frameViews))
            }
            else {
                dgm.frameModels@models[[m]] <<- dm.frameViews
            }
        if (!redraw && !returnNewMaster) {
            local.graphWindow <- redrawView(frameModels = dgm.frameModels, 
                frameViews = dm.frameViews, graphWindow = graphWindow, 
                dg = dg, initialWindow = initialWindow, returnNewMaster = FALSE, 
                redraw = FALSE, returnFrameModel = FALSE, control = control)
        }
        if (returnFrameModel) 
            invisible(dgm.frameModels)
        else invisible(dm.frameViews)
    }
    Arguments <- list(...)
    dgm.frameModels <- NULL
    if (!is.null(Arguments$dgm.frameModels)) 
        dgm.frameModels <- Arguments$dgm.frameModels
    if (!is.null(Arguments$frameModels)) 
        dgm.frameModels <- Arguments$frameModels
    returnNewMaster <- FALSE
    if (!is.null(Arguments$returnNewMaster)) 
        returnNewMaster <- Arguments$returnNewMaster
    redraw <- FALSE
    if (!is.null(Arguments$redraw)) 
        redraw <- Arguments$redraw
    updateCountVerticesMain <- 0
    updateCountPositionsMain <- 0
    updateCountBlocksMain <- 0
    updateCountBlockEdgesMain <- 0
    if (is.null(dgm.frameModels)) {
        dgm.frameModels.env <- .Dg.toplevel()
        dgm.frameModels <- .newDynamicGraphObject(vertexList, 
            blocks = blockList, blockTree = blockTree, label = control$label, 
            control = control, parent = .DgRoot, env = dgm.frameModels.env)
        if (control$saveFunctions) 
            assign("DrawModel", drawModel, envir = dgm.frameModels.env$env)
    }
    else {
        dgm.frameModels.env <- .Dg.toplevel()
        dgm.frameModels@id.env <- dgm.frameModels.env$ID
    }
    if (redraw) {
        vertexList <- dgm.frameModels@vertices
        blockList <- dgm.frameModels@blocks
        if (length(blockList) == 0) 
            blockList <- NULL
    }
    namesVertices <- Names(vertexList)
    positionsVertices <- Positions(vertexList)
    if (is.matrix(positionsVertices)) {
        indices <- Indices(vertexList)
        if (!all(seq(length(indices)) == indices)) {
            warning("Invalid indices of vertices replaces")
            Indices(vertexList) <- seq(length(indices))
            indices <- Indices(vertexList)
        }
        positionsLabels <- positionsVertices
        positionsLabels[, 1] <- positionsLabels[, 1] + 0.1 * 
            control$w
        Labels <- Labels(vertexList)
        constrainedVertices <- Constrained(vertexList)
        colorsVertices <- Colors(vertexList)
        blocksVertices <- Blockindices(vertexList)
        strataVertices <- Strata(vertexList)
        local.N <- ncol(positionsVertices)
        if (is.null(blockList) && !is.null(blockTree)) 
            blockList <- blockTreeToList(blockTree)
        positionsBlock <- NULL
        positionsBlockLabels <- NULL
        positionsClosedBlocks <- NULL
        blockLabels <- NULL
        strataBlocks <- NULL
        if (!.IsEmpty(blockList)) {
            positionsBlocks <- Positions(blockList)
            d <- dim(positionsBlocks)
            positionsBlocks <- array(positionsBlocks, dim = c(d[1], 
                d[2]/2, 2))
            positionsClosedBlocks <- matrix(rep(NA, local.N * 
                length(blockList)), ncol = local.N)
            positionsClosedBlocks <- apply(positionsBlocks, c(1, 
                2), mean)
            blockReferences <- 1:length(blockList)
            positionsBlockLabels <- matrix(rep(0, local.N * length(blockList)), 
                ncol = local.N)
            blockLabels <- Labels(blockList)
            strataBlocks <- Strata(blockList)
        }
        if (returnNewMaster || redraw) {
            models <- dgm.frameModels@models
            if (!redraw) 
                dgm.frameModels@models <- list()
            for (i in (1:length(models))) {
                frameViews <- models[[i]]
                object <- models[[i]]@model
                R <- drawModel(frameModels = dgm.frameModels, 
                  frameViews = frameViews, graphWindow = NULL, 
                  object = object, frameModelsEnv = dgm.frameModels.env, 
                  initialWindow = TRUE, returnNewMaster = returnNewMaster, 
                  redraw = redraw, returnFrameModel = FALSE, 
                  control = control)
            }
        }
        if (!redraw) {
            local.frameViews <- drawModel(frameModels = dgm.frameModels, 
                frameViews = NULL, graphWindow = NULL, frameModelsEnv = dgm.frameModels.env, 
                dg = dg, object = object, initialWindow = TRUE, 
                returnNewMaster = FALSE, redraw = FALSE, returnFrameModel = FALSE, 
                control = control)
        }
    }
    else {
        dgm.frameModels <- NULL
        warning("Positions of vertices should have same number of coordinates")
    }
    if (control$returnNull) 
        return(NULL)
    else invisible(dgm.frameModels)
}
