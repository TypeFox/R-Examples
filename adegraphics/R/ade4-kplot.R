"kplot.mcoa" <- function(object, xax = 1, yax = 2, which.tab = 1:nrow(object$cov2), option = c("points", "axis", "columns"), pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "mcoa")) 
    stop("Object of class 'mcoa' expected")
  if((xax == yax) || (object$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > object$nf)
    stop("Non convenient xax")
  if(yax > object$nf)
    stop("Non convenient yax")
  
  option <- match.arg(option)
  
  ## parameters management
  sortparameters <- sortparamADEg(...)
  
  if(option == "points") {
    params1 <- list()
    params1$adepar <- list(psub = list(text = "Reference"), plabels = list(cex = 1.25))
    sortparameters1 <- modifyList(params1, sortparameters, keep.null = TRUE)
    ref <- do.call("s.label", c(list(dfxy = substitute(object$SynVar), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters1$adepar, sortparameters1$trellis, sortparameters1$g.args))
    
    params2 <- list()
    params2$adepar <- list(plabels = list(cex = 0))
    params2$g.args <- list(samelimits = FALSE)
    sortparameters2 <- modifyList(params2, sortparameters, keep.null = TRUE)
    
    facets1 <- substitute(object$TL[,1])
    coolig <- call("as.data.frame", call("matrix", call("kronecker", rep(1,nrow(object$cov2)), substitute(as.matrix(object$SynVar))), nrow = nrow(object$Tl1), ncol = ncol(object$Tl1), dimnames = substitute(list(rownames(object$Tl1), colnames(object$Tl1)))))
    g1 <- do.call("s.match", c(list(dfxy1 = coolig, dfxy2 = substitute(object$Tl1), facets = facets1, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters2$adepar, sortparameters2$trellis, sortparameters2$g.args))[which.tab]
    
    ## ADEgS creation
    ADEglist <- c(list(ref), g1@ADEglist)
    nrow_lay  <- floor(sqrt(length(ADEglist))) + 1
    ncol_lay <- -floor(-(length(ADEglist)) / nrow_lay)
    lay <- matrix(c(seq(1, length(ADEglist)), rep(0, nrow_lay * ncol_lay - length(ADEglist))), nrow = nrow_lay, byrow = TRUE)
    obj <- new(Class = "ADEgS", ADEglist = ADEglist, positions = layout2position(lay), add = matrix(0, ncol = length(ADEglist), nrow = length(ADEglist)), Call = match.call())
    names(obj) <- c("ref", names(g1))
    
  } else if(option == "axis") {
    params <- list()
    params$adepar <- list(pbackground = list(box = FALSE), plabels = list(cex = 1.25))
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    facets2 <- substitute(object$T4[, 1])
    obj <- do.call("s.corcircle", c(list(dfxy = substitute(object$Tax), facets = facets2, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args))[which.tab]
    
  } else if(option == "columns") {
    params <- list()
    params$adepar <- list(plabels = list(cex = 1.25))
    params$g.args <- list(samelimits = FALSE)
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    facets3 <- substitute(object$TC[, 1])
    obj <- do.call("s.arrow", c(list(dfxy = substitute(object$Tco), facets = facets3, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args))[which.tab]
  }
  
  obj@Call <- match.call()
  if(plot) 
    print(obj)
  invisible(obj)
}


"kplot.mfa" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$blo), traject = FALSE, permute = FALSE, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "mfa")) 
    stop("Object of class 'mfa' expected")
  if((xax == yax) || (object$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > object$nf)
    stop("Non convenient xax")
  if(yax > object$nf)
    stop("Non convenient yax")
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "traj")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list(plabels = list(cex = 0), ppoints = list(cex = 1.5), samelimits = FALSE)
  params$col <- list(psub = list(cex = 0), plabels = list(cex = 1.25))
  params$traj <- list(plabels = list(cex = 0), psub = list(cex = 0))
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## prepare
  if(permute) {
    dfxy_row <- substitute(object$co)
    dfxy_col <- substitute(object$lisup)
    facets_row <- substitute(object$TC[,1])
    facets_col <- substitute(object$TL[,1])
  } else {
    dfxy_row <- substitute(object$lisup)
    dfxy_col <- substitute(object$co)
    facets_row <- substitute(object$TL[,1])
    facets_col <- substitute(object$TC[,1])
  }
  
  ## create g1
  g1 <- do.call("s.label", c(list(dfxy = dfxy_row, facets = facets_row, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))[which.tab]
  
  ## prepare and create g2
  if(permute)
    dcol <- object$lisup
  else
    dcol <- object$co
  k <- c(min(dcol[, xax]), max(dcol[, xax]), min(dcol[, yax]), max(dcol[, yax])) / c(g1[[1]]@g.args$xlim, g1[[1]]@g.args$ylim)
  dcol <- substitute(dfxy_col * 0.7 / max(k))
  g2 <- do.call("s.arrow", c(list(dfxy = dcol, facets = facets_col, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]
  obj <- do.call("superpose", list(g1, g2))
  obj@Call <- call("superpose", g1@Call, g2@Call)
  
  ## create g3
  if(traject) {
    g3 <- do.call("s.traject", c(list(dfxy = dfxy_row, facets = facets_row, xax = xax, yax = yax, plot = FALSE, storeData = FALSE, pos = pos - 2), sortparameters$traj))[which.tab]
    obj <- do.call("superpose", list(obj, g3))
    obj@Call <- call("superpose", obj@Call, g3@Call)
  }
  
  ## ADEgS creation
  names(obj) <- object$tab.names[which.tab]
  obj@Call <- match.call()
  if(plot) 
    print(obj)
  invisible(obj)
}


"kplot.pta" <- function(object, xax = 1, yax = 2, which.tab = 1:nrow(object$RV), which.graph = 1:4, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "pta")) 
    stop("Object of class 'pta' expected")
  if((xax == yax) || (object$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  if(!is.numeric(which.graph) || any(which.graph < 1) || any(which.graph > 4)) 
    stop("'which' must be in 1:4")
  
  if(xax > object$nf)
    stop("Non convenient xax")
  if(yax > object$nf)
    stop("Non convenient yax")
  
  ## sort parameters for each graph
  graphsnames <- c("axis", "row", "col", "components")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$axis <- list(pbackground = list(box = FALSE), plabels = list(alpha = 1, cex = 1.25))
  params$rows <- list(plabels = list(alpha = 1, cex = 1.25))
  params$columns <- list(plabels = list(cex = 1.25))
  params$components <- list(pbackground = list(box = FALSE), plabels = list(cex = 1.25))
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  g <- as.null()
  adeglist <- as.null()
  ## creation of each individual ADEg
  if(1 %in% which.graph) {
  	facets1 <- substitute(object$T4[, 1])
  	g1 <- do.call("s.corcircle", c(list(dfxy = substitute(object$Tax), labels = substitute(object$T4[, 2]), facets = facets1, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$axis))[which.tab]
  	names(g1) <- paste(graphsnames[1], "_", object$tab.names, sep = "")[which.tab]
    g <- c(g, g1)
    adeglist <- c(adeglist, g1@ADEglist)
  }
  
  if(2 %in% which.graph) {
    facets2 <- substitute(object$TL[, 1])
    g2 <- do.call("s.label", c(list(dfxy = substitute(object$Tli), labels = substitute(object$TL[,2]), facets = facets2, xax = 1, yax = 2, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$rows))[which.tab]
    names(g2) <- paste(graphsnames[2], "_", object$tab.names, sep = "")[which.tab]
    g <- c(g, g2)
    adeglist <- c(adeglist, g2@ADEglist)
  }
  
  if(3 %in% which.graph) {
    facets3 <- substitute(object$TC[, 1])
  	g3 <- do.call("s.arrow", c(list(dfxy = substitute(object$Tco), labels = substitute(object$TC[,2]), facets = facets3, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$columns))[which.tab]
  	names(g3) <- paste(graphsnames[3], "_", object$tab.names, sep = "")[which.tab]
    g <- c(g, g3)
    adeglist <- c(adeglist, g3@ADEglist)
  }
  
  if(4 %in% which.graph) {
  	facets4 <- substitute(object$T4[, 1])
  	g4 <- do.call("s.corcircle", c(list(dfxy = substitute(object$Tcomp), labels = substitute(object$T4[, 2]), facets = facets4, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$components))[which.tab]
  	names(g4) <- paste(graphsnames[4], "_", object$tab.names, sep = "")[which.tab]
    g <- c(g, g4)
    adeglist <- c(adeglist, g4@ADEglist)
  }
  
  ## ADEgS creation
  ng <- sum(sapply(g, function(x) length(x)))
  lay <- matrix(1:ng, ncol = length(which.graph))
  obj <- new(Class = "ADEgS", ADEglist = c(adeglist), positions = layout2position(lay), add = matrix(0, ncol = ng, nrow = ng), Call = match.call())
  if(plot) 
    print(obj)
  invisible(obj)
}


"kplot.sepan" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$blo), permute = FALSE, traject = FALSE, posieig = "bottomleft", pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "sepan")) 
    stop("Object of class 'sepan' expected")
  if((xax == yax) || (length(object$Eig) == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > length(object$Eig))
    stop("Non convenient xax")
  if(yax > length(object$Eig))
    stop("Non convenient yax")
  
  ## prepare
  if(permute) {
    dfxy_row <- substitute(object$Co)
    dfxy_col <- substitute(object$Li)
    names_row <- substitute(object$TC[,2])
    names_col <- substitute(object$TL[,2])
    facets_row <- substitute(object$TC[,1])
    facets_col <- substitute(object$TL[,1])
  } else {
    dfxy_row <- substitute(object$Li)
    dfxy_col <- substitute(object$Co)
    names_row <- substitute(object$TL[,2])
    names_col <- substitute(object$TC[,2])
    facets_row <- substitute(object$TL[,1])
    facets_col <- substitute(object$TC[,1])
  }
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "traj", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list(psub = list(position = "bottomright"), samelimits = FALSE)
  params$traj <- list(psub = list(position = "bottomright"), plabels = list(cex = 0), samelimits = FALSE)
  params$col <- list(psub = list(cex = 0, position = "bottomright"), plabels = list(cex = 1.25))
  params$eig <- list(psub = list(text = ""), pbackground = list(box = TRUE), samelimits = FALSE)
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## create g1
  if(!traject) 
  	g1 <- do.call("s.label", c(list(dfxy = dfxy_row, labels = names_row, facets = facets_row, xax = 1, yax = 2, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))[which.tab]
  else 
    g1 <- do.call("s.traject", c(list(dfxy = dfxy_row, facets = facets_row, xax = xax, yax = yax, plot = FALSE, storeData = FALSE, pos = pos - 2), sortparameters$traj))[which.tab]
  
  ## prepare and create g2
  if(permute)
    dcol <- object$Li
  else
    dcol <- object$Co
  k <- c(min(dcol[, xax]), max(dcol[, xax]), min(dcol[, yax]), max(dcol[, yax])) / c(g1[[1]]@g.args$xlim, g1[[1]]@g.args$ylim)
  dcol <- substitute(dfxy_col * 0.7 / max(k))
  g2 <- do.call("s.arrow", c(list(dfxy = dcol, labels = names_col, facets = facets_col, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]
  obj <- do.call("superpose", list(g1, g2))
  obj@Call <- call("superpose", g1@Call, g2@Call)
  
  ## prepare and create g3
  facets_eig <- reorder(as.factor(rep(levels(object$TL[, 1]), object$rank)), rep(1:length(object$rank), object$rank))
  if(!any(posieig == "none")) {
    g3 <- do.call("plotEig", c(list(eigvalue = substitute(object$Eig), nf = 1:ncol(object$Li), facets = facets_eig, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))[which.tab]
    obj <- do.call("insert", list(g3, obj, posi = posieig, plot = FALSE, ratio = 0.2, inset = 0, dispatch = TRUE))
  }
  
  ## ADEgS creation
  names(obj) <- object$tab.names[which.tab]
  obj@Call <- match.call()
  if(plot) 
    print(obj)
  invisible(obj)
} 


"kplotsepan.coa" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$blo), permute = FALSE, posieig = "bottomleft", pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "sepan")) 
    stop("Object of class 'sepan' expected")
  if((xax == yax) || (length(object$Eig) == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > length(object$Eig))
    stop("Non convenient xax")
  if(yax > length(object$Eig))
    stop("Non convenient yax")
  
  ## prepare
  if(permute) {
    dfxy_row <- substitute(object$C1)
    dfxy_col <- substitute(object$Li)
    names_row <- substitute(object$TC[,2])
    names_col <- substitute(object$TL[,2])
    facets_row <- substitute(object$TC[,1])
    facets_col <- substitute(object$TL[,1])
  } else {
    dfxy_row <- substitute(object$Li)
    dfxy_col <- substitute(object$C1)
    names_row <- substitute(object$TL[,2])
    names_col <- substitute(object$TC[,2])
    facets_row <- substitute(object$TL[,1])
    facets_col <- substitute(object$TC[,1])
  }
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$col <- list(psub = list(position = "bottomright"), plabels = list(cex = 1.25), samelimits = FALSE)
  params$row <- list(psub = list(cex = 0, position = "bottomright"))
  params$eig <- list(psub = list(text = ""), pbackground = list(box = TRUE))
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation and create g1 and g2
  g1 <- do.call("s.label", c(list(dfxy = dfxy_col, labels = names_col, facets = facets_col, xax = 1, yax = 2, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]
  g2 <- do.call("s.label", c(list(dfxy = dfxy_row, labels = names_row, facets = facets_row, xax = 1, yax = 2, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))[which.tab]
  obj <- do.call("superpose", c(list(g1, g2)))
  obj@Call <- call("superpose", g1@Call, g2@Call)
  
  ## prepare and create g3
  facets_eig <- reorder(as.factor(rep(levels(object$TL[, 1]), object$rank)), rep(1:length(object$rank), object$rank))
  if(!any(posieig == "none")) {
    g3 <- do.call("plotEig", c(list(eigvalue = substitute(object$Eig), nf = 1:ncol(object$Li), facets = facets_eig, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))[which.tab]
    obj <- do.call("insert", list(g3, obj, posi = posieig, plot = FALSE, ratio = 0.2, inset = 0, dispatch = TRUE))
  }
  
  ## ADEgS creation
  names(obj) <- object$tab.names[which.tab]
  obj@Call <- match.call()
  if(plot) 
    print(obj)
  invisible(obj)
}


"kplot.statis" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$tab.names), traject = FALSE, arrow = TRUE, class = NULL, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "statis")) 
    stop("Object of class 'statis' expected")
  if((xax == yax) || (object$C.nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > object$C.nf)
    stop("Non convenient xax")
  if(yax > object$C.nf)
    stop("Non convenient yax")
  
  ## sort parameters for each graph
  graphsnames <- c("col", "traj", "class")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$col <- list(plabels = list(cex = 1.25))
  params$traj <- list(plabels = list(cex = 0), psub = list(cex = 0))
  params$class <- list(plabels = list(cex = 1.5), ppoints = list(cex = 2), pellipses = list(alpha = 0, axes = list(draw = FALSE)), psub = list(cex = 0))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## prepare
  facets <- substitute(object$TC[, 1])
  
  ## creation of each individual ADEg
  if(arrow) 
    g1 <- do.call("s.arrow", c(list(dfxy = substitute(object$C.Co), facets = facets, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]    	
  else
    g1 <- do.call("s.label", c(list(dfxy = substitute(object$C.Co), xax = xax, yax = yax, facets = facets, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]
  
  if(traject) {
    g2 <- do.call("s.traject", c(list(dfxy = substitute(object$C.Co), xax = xax, yax = yax, facets = facets, plot = FALSE, storeData = FALSE, pos = pos - 2), sortparameters$traj))[which.tab]
    obj <- do.call("superpose", list(g1, g2))
    obj@Call <- call("superpose", g1@Call, g2@Call)
  } else
    obj <- g1
  
  if(!is.null(class)) {
    if(length(class) == 1) {
      if(class)
	    	g3 <- do.call("s.class", c(list(dfxy = substitute(object$C.Co), fac = object$TC[, 1], xax = xax, yax = yax, facets = facets, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$class))[which.tab]
    } else {
      if(length(class) == length(object$TC[, 1]))
        g3 <- do.call("s.class", c(list(dfxy = substitute(object$C.Co), fac = factor(class), xax = xax, yax = yax, facets = facets, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$class))[which.tab]
    }
    obj <- do.call("superpose", list(obj, g3))
  	obj@Call <- call("superpose", g3@Call, obj@Call)
  }
  
  
  ## ADEgS creation
  names(obj) <- object$tab.names[which.tab]
  obj@Call <- match.call()
  if(plot) 
    print(obj)
  invisible(obj)
}


"kplot.foucart" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$blo), pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "foucart")) 
    stop("Object of class 'foucart' expected")
  if((xax == yax) || (object$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > object$nf)
    stop("Non convenient xax")
  if(yax > object$nf)
    stop("Non convenient yax")
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## limits calcul
  df <- rbind(as.matrix(object$li), as.matrix(object$Tli), as.matrix(object$Tco))
  adegtot <- adegpar()
  lim.global <- setlimits2D(minX = min(df[, xax]), maxX = max(df[, xax]), minY = min(df[, yax]), maxY = max(df[, yax]), origin = adegtot$porigin$origin, aspect.ratio = adegtot$paxes$aspectratio, includeOr = adegtot$porigin$include)
  
  ## parameters management
  params <- list()
  params$row <- list(plabels = list(cex = 1), xlim = lim.global$xlim, ylim = lim.global$ylim, plabels = list(cex = 1.25))
  params$col <- list(plabels = list(cex = 1.25), psub = list(text = ""), xlim = lim.global$xlim, ylim = lim.global$ylim, plabels = list(cex = 1.25))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation of each individual ADEg
  g1 <- do.call("s.label", c(list(dfxy = substitute(object$Tli), facets =  substitute(object$TL[, 1]), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))[which.tab]
  g2 <- do.call("s.label", c(list(dfxy = substitute(object$Tco), facets = substitute(object$TC[, 1]), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))[which.tab]
  
  ## ADEgS creation
  obj <- do.call("superpose", list(g1, g2))
  names(obj) <- object$tab.names
  obj@Call <- match.call()
  if(plot)
    print(obj)
  invisible(obj)
}


"kplot.mbpcaiv" <- function(object, xax = 1, yax = 2, which.tab = 1:length(object$blo), pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(object, "mbpcaiv")) 
    stop("Object of class 'mbpcaiv' expected")
  if((xax == yax) || (object$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > object$nf)
    stop("Non convenient xax")
  if(yax > object$nf)
    stop("Non convenient yax")
  
  sortparameters <- sortparamADEg(...)
  
  obj <- do.call("s.label", c(list(dfxy = substitute(object$Tli), xax = xax, yax = yax, facets = substitute(object$TL[, 1]), plot = plot, storeData = storeData, pos = pos - 2), adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args))[which.tab]
  
  obj@Call <- match.call()
  if(plot)
    print(obj)
  invisible(obj)
}
