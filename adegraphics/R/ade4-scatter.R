"scatter.dudi" <- function(x, xax = 1, yax = 2, permute = FALSE, posieig = "topleft", prop = FALSE, 
  density.plot = ifelse(permute, ncol(x$tab) > 1000, nrow(x$tab) > 1000), plot = TRUE, storeData = TRUE, pos = -1, ...) {
  if(!inherits(x, "dudi")) 
    stop("Object of class 'dudi' expected")
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")  
  
  position <- match.arg(posieig[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list(plabels = list(cex = 0.75))
  params$col <- list()
  params$eig <- list(pbackground = list(box = TRUE), psub = list(text = "Eigenvalues"))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  if(prop) {
    id <- inertia.dudi(x, col.inertia = TRUE)
    if(is.null(sortparameters[[2]]$plabels$cex)) {
      sortparameters$col$plabels$cex <- id$col.cum[, 2] / (max(id$col.cum[, 2]) / 1.5)
    } else {
	    sortparameters$col$plabels$cex <- sortparameters$col$plabels$cex * id$col.cum[, 2] / (max(id$col.cum[, 2]) / 1.5)
    }
  }
  
  ## prepare and create g1
  if(permute)
    df1 <- substitute(x$co)
  else
    df1 <- substitute(x$li)
  g1 <- do.call(ifelse(density.plot, "s.density", "s.label"), c(list(dfxy = df1, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
  
  ## prepare and create g2
  if(permute) {
    colss <- x$l1
  } else {
    colss <- x$c1
  }
  knormali <- c(min(colss[, xax]), max(colss[, xax]), min(colss[, yax]), max(colss[, yax])) / c(g1@g.args$xlim, g1@g.args$ylim)
  csts <- 0.9 / max(knormali)
  if(permute) {
    df2 <- substitute(x$l1 * csts)
  } else {
    df2 <- substitute(x$c1 * csts)
  }
  g2 <- do.call("s.arrow", c(list(dfxy = df2, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))
  
  ## create the final ADEgS
  object <- do.call("superpose", list(g1, g2))
  object@Call <- call("superpose", g1@Call, g2@Call)
  if(position != "none") {
    g3 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = 1:x$nf, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))
    object <- do.call("insert", list(g3@Call, object@Call, posi = position, plot = FALSE, ratio = 0.25))
  }
  
  names(object) <- graphsnames[1:length(object)]
  object@Call <- match.call()
  if(plot) 
    print(object)
  invisible(object)
}


"scatter.coa" <- function(x, xax = 1, yax = 2, method = 1:3, posieig = "topleft", pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "dudi"))
    stop("Object of class 'dudi' expected")
  if(!inherits(x, "coa"))
    stop("Object of class 'coa' expected")
  
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")
  
  position <- match.arg(posieig[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
  method <- method[1]
  
  ## limits management
  if(method == 1)
    x.global <- rbind(as.matrix(x$li), as.matrix(x$co))
  else if(method == 2)
    x.global <- rbind(as.matrix(x$c1), as.matrix(x$li))
  else if(method == 3)
    x.global <- rbind(as.matrix(x$l1), as.matrix(x$co))
  adegtot <- adegpar()
  lim.global <- setlimits2D(minX = min(x.global[, xax]), maxX = max(x.global[, xax]), minY = min(x.global[, yax]), maxY = max(x.global[, yax]), origin = adegtot$porigin$origin, aspect.ratio = adegtot$paxes$aspectratio, includeOr = adegtot$porigin$include)
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list(plabels = list(cex = 0.75), xlim = lim.global$xlim, ylim = lim.global$ylim)
  params$col <- list(xlim = lim.global$xlim, ylim = lim.global$ylim)
  params$eig <- list(pbackground = list(box = TRUE), psub = list(text = "Eigenvalues"))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation of each individual ADEg and of the final ADEgS
  if(method == 1) {
    g1 <- do.call("s.label", c(list(dfxy = substitute(x$li), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
    g2 <- do.call("s.label", c(list(dfxy = substitute(x$co), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))
  } else if(method == 2) {
    g1 <- do.call("s.label", c(list(dfxy = substitute(x$c1), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))
    g2 <- do.call("s.label", c(list(dfxy = substitute(x$li), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
  } else if(method == 3) {
    g1 <- do.call("s.label", c(list(dfxy = substitute(x$l1), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
    g2 <- do.call("s.label", c(list(dfxy = substitute(x$co), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))
  }  
  object <- do.call("superpose", list(g1, g2))
  object@Call <- call("superpose", g1@Call, g2@Call)
  if(position != "none") {
    g3 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = 1:x$nf, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))
    object <- do.call("insert", list(g3@Call, object@Call, posi = position, plot = FALSE, ratio = 0.25))
  }
  
  object@Call <- match.call()
  names(object) <- graphsnames[1:length(object)]
  if(plot) 
    print(object)
  invisible(object)
}


"plot.acm"  <- function(x, xax = 1, yax = 2, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "dudi"))
    stop("Object of class 'dudi' expected")
  if(!inherits(x, "acm"))
    stop("Object of class 'acm' expected")
  
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")
  
  ## prepare
  oritab <- as.list(x$call)[[2]]
  
  ## parameter management
  sortparameters <- sortparamADEg(...)
  params <- list()
  params$g.args <- list(starSize = 0)
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  object <- do.call("s.class", c(list(dfxy = substitute(x$li), fac = oritab, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
  object@Call <- match.call()
  if(plot)
    print(object)
  invisible(object)
}


"plot.fca" <- function(x, xax = 1, yax = 2, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "dudi"))
    stop("Object of class 'dudi' expected")
  if(!inherits(x, "fca"))
    stop("Object of class 'fca' expected")
  
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")
  
  ## prepare
  oritab <- as.list(x$call)[[2]]
  evTab <- eval.parent(oritab)
  indica <- factor(rep(names(x$blo), x$blo))
  ng <- length(levels(indica))   
  
  ## parameter management
  graphsnames <- as.character(levels(indica))
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  params <- list()
  params <- lapply(1:length(graphsnames), function(i) {params[[i]] <- list(starSize = 0.5, ellipseSize = 0, plabels = list(cex = 1.25), psub = list(text = graphsnames[i]))})
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation of each individual ADEg
  l <- list()
  l <- sapply(1:length(levels(indica)), function(i) {do.call("s.distri", c(list(dfxy = substitute(x$l1, env = sys.frame(-3)), dfdistri = call("[", oritab, call(":", 1, nrow(evTab)), which(indica == levels(indica)[i])), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters[[i]]))})
  
  ## ADEgS creation
  object <- new(Class = "ADEgS", ADEglist = l, positions = layout2position(.n2mfrow(ng), ng  = ng), add = matrix(0, ncol = ng, nrow = ng), Call = match.call())
  names(object) <- graphsnames
  if(plot)
    print(object)
  invisible(object)
}


"scatter.pco" <- function(x, xax = 1, yax = 2, posieig = "topleft", pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "dudi"))
    stop("Object of class 'dudi' expected")
  if(!inherits(x, "pco"))
    stop("Object of class 'pco' expected")
  
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")
  
  position <- match.arg(posieig[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
  
  ## sort parameters for each graph
  graphsnames <- c("row", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list()
  params$eig <- list(pbackground = list(box = TRUE), psub = list(text = "Eigenvalues"))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation of each individual ADEg and of the final ADEgS
  object <- do.call("s.label", c(list(dfxy = substitute(x$li), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
  if(position != "none") {
    g2 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = 1:x$nf, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))
    object <- do.call("insert", list(g2@Call, object@Call, posi = position, plot = FALSE, ratio = 0.25))
    names(object) <- graphsnames[1:length(object)]
  }
  
  object@Call <- match.call()
  if(plot)
    print(object)
  invisible(object)
}


"scatter.nipals" <- function(x, xax = 1, yax = 2, posieig = "topleft", pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "nipals"))
    stop("Object of class 'nipals' expected")
  
  if((xax == yax) || (x$nf == 1))
    stop("One axis only : not yet implemented")
  if(length(xax) > 1 | length(yax) > 1)
    stop("Not implemented for multiple xax/yax")
  
  if(xax > x$nf)
    stop("Non convenient xax")
  if(yax > x$nf)
    stop("Non convenient yax")
  
  position <- match.arg(posieig[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
  
  ## sort parameters for each graph
  graphsnames <- c("row", "col", "eig")
  sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
  
  ## parameters management
  params <- list()
  params$row <- list(plabels = list(cex = 0.75))
  params$col <- list()
  params$eig <- list(pbackground = list(box = TRUE), psub = list(text = "Eigenvalues"))
  names(params) <- graphsnames
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## prepare and create g1
  g1 <- do.call("s.label", c(list(dfxy = substitute(x$li), xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$row))
  
  ## prepare and create g2
  knormali <- c(min(x$c1[, xax]), max(x$c1[, xax]), min(x$c1[, yax]), max(x$c1[, yax])) / c(g1@g.args$xlim, g1@g.args$ylim)
  csts <- 0.8 / max(knormali)
  df2 <- substitute(x$c1 * csts)
  g2 <- do.call("s.arrow", c(list(dfxy = df2, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$col))
  
  ## creation of each individual ADEg and of the final ADEgS
  object <- do.call("superpose", list(g1, g2))
  object@Call <- call("superpose", g1@Call, g2@Call)
  if(position != "none") {
    g3 <- do.call("plotEig", c(list(eigvalue = substitute(x$eig), nf = 1:x$nf, xax = xax, yax = yax, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$eig))
    object <- do.call("insert", list(g3@Call, object@Call, posi = position, plot = FALSE, ratio = 0.25))
  }
  
  names(object) <- graphsnames[1:length(object)]
  object@Call <- match.call()
  if(plot)
    print(object)
  invisible(object)
}


"score.acm" <- function (x, xax = 1, which.var = NULL, type = c("points", "boxplot"), pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "acm")) 
    stop("Object of class 'acm' expected")
  if(x$nf == 1) 
    xax <- 1
  if((xax < 1) || (xax > x$nf)) 
    stop("non convenient axe number")
  
  ## prepare
  oritab <- as.list(x$call)[[2]]
  evTab <- eval.parent(oritab)
  if(is.null(which.var))
    which.var <- 1:ncol(evTab)
  
  type <- match.arg(type)
  
  ## parameter management
  sortparameters <- sortparamADEg(...)
  params <- list()
  
  if(type == "boxplot") {
    ## parameter management
    params$adepar <- list(plabels = list(boxes = list(draw = FALSE)), p1d = list(rug = list(draw = TRUE)), paxes = list(draw = TRUE, y = list(draw = FALSE)), 
      plegend = list(drawKey = FALSE), pgrid = list(text = list(cex = 0)), psub = list(position = "topleft"))
    params$g.args <- list(samelimits = FALSE)
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    ## ADEgS creation
    scorecall <- substitute(x$l1[, xax])
    fac <- call("[", oritab, which.var)
    object <- do.call("s1d.boxplot", c(list(score = scorecall, fac = fac, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
    
  } else if(type == "points") {
    ## parameter management
    params$adepar <- list(ppoints = list(pch = "|"), porigin = list(draw = FALSE), pgrid = list(draw = FALSE), psub = list(position = "topleft"), paxes = list(draw = TRUE), plabels = list(cex = 1.25))
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    
    ## creation of each individual ADEg
    ADEglist <- list()
    score <- x$l1[, xax]
    scorecall <- substitute(x$l1[, xax])
    for(i in which.var) {
      ## data management
      fac <- evTab[, i]
      faccall <- call("[", oritab, 1:NROW(evTab), i)
      meangroup <- call("as.numeric", call("tapply", scorecall, faccall, mean))
      dfxy <- call("cbind", scorecall, call("as.numeric", call("[", meangroup, faccall)))
      
      ## ADEg creation
      g1 <- do.call("s.class", c(list(dfxy = dfxy, fac = faccall, plot = FALSE, storeData = storeData, pos = pos - 2), c(sortparameters$adepar, list(psub.text = colnames(evTab)[i])), sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
      xlimg1 <- g1@g.args$xlim
      ylimg1 <- g1@g.args$ylim
      g2 <- xyplot(score ~ fac, xlab = "", ylab = "", scales = list(x = list(tck = c(1, 0)), y = list(tck = c(1, 0))), xlim = xlimg1, ylim = ylimg1, 
                   aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(h = as.numeric(tapply(y, x, mean)), a = 0, b = 1, lty = 1)})
      g2$call <- call("xyplot", substitute(scorecall ~ faccall), xlab = "", ylab = "", scales = list(x = list(tck = c(1, 0)), y = list(tck = c(1, 0))), xlim = substitute(xlimg1), ylim = substitute(ylimg1),
                      aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(h = as.numeric(tapply(y, x, mean)), a = 0, b = 1, lty = 1)})
      ADEglist[[i]] <- superpose(g2, g1, plot = FALSE)
    }
    ADEglist <- ADEglist[which.var]
    
    ## ADEgS creation
    posmatrix <- layout2position(.n2mfrow(length(which.var)), ng = length(which.var))
    object <- new(Class = "ADEgS", ADEglist = ADEglist, positions = posmatrix, add = matrix(0, ncol = length(which.var), nrow = length(which.var)), Call = match.call())
  } 
  
  names(object) <- colnames(evTab)[which.var]
  object@Call <- match.call()
  if(plot) 
    print(object)
  invisible(object)
}


#"score.coa" <- function (x, xax = 1, dotchart = FALSE, pos = -1, storeData = TRUE, plot = TRUE, ...) {
#
#  if(!inherits(x, "coa")) 
#    stop("Object of class 'coa' expected")
#  if(x$nf == 1) 
#    xax <- 1
#  if((xax < 1) || (xax > x$nf)) 
#    stop("non convenient axe number")
#  
#  if(dotchart)
#    stop("TRUE 'dotchart' not yet implemented")
#  
#  
#
#  def.par <- par(mar = par("mar"))
#  on.exit(par(def.par))
#  par(mar = c(0.1, 0.1, 0.1, 0.1))
#  
#  sco.distri.class.2g <- function(score, fac1, fac2, weight, labels1 = as.character(levels(fac1)), labels2 = as.character(levels(fac2)), clab1, clab2, cpoi, cet) {
#    nvar1 <- nlevels(fac1)
#    nvar2 <- nlevels(fac2)
#    ymin <- scoreutil.base(y = score, xlim = NULL, grid = TRUE, cgrid = 0.75, include.origin = TRUE, origin = 0, sub = NULL, csub = 0)
#    ymax <- par("usr")[4]
#    ylabel <- strheight("A", cex = par("cex") * max(1, clab1, clab2)) * 1.4
#    xmin <- par("usr")[1]
#    xmax <- par("usr")[2]
#    xaxp <- par("xaxp")
#    nline <- xaxp[3] + 1
#    v0 <- seq(xaxp[1], xaxp[2], le = nline)
#    
#    ## dessine la grille
#    segments(v0, rep(ymin, nline), v0, rep(ymax, nline), col = gray(0.5), lty = 1)
#    
#    ## dessine le cadre
#    rect(xmin, ymin, xmax, ymax)
#    
#    
#    sum.col1 <- unlist(tapply(weight, fac1, sum))
#    sum.col2 <- unlist(tapply(weight, fac2, sum))
#    sum.col1[sum.col1 == 0] <- 1
#    sum.col2[sum.col2 == 0] <- 1
#    
#    weight1 <- weight/sum.col1[fac1]
#    weight2 <- weight/sum.col2[fac2]
#    
#    y.distri1 <- tapply(score * weight1, fac1, sum)
#    y.distri1 <- rank(y.distri1)
#    y.distri2 <- tapply(score * weight2, fac2, sum)
#    y.distri2 <- rank(y.distri2) + nvar1 + 2
#    y.distri <- c(y.distri1, y.distri2)
#    
#    ylabel <- strheight("A", cex = par("cex") * max(1, clab1, clab2)) * 1.4
#    y.distri1 <- (y.distri1 - min(y.distri))/(max(y.distri) - min(y.distri))
#    y.distri1 <- ymin + ylabel + (ymax - ymin - 2 * ylabel) * y.distri1
#    y.distri2 <- (y.distri2 - min(y.distri))/(max(y.distri) - min(y.distri))
#    y.distri2 <- ymin + ylabel + (ymax - ymin - 2 * ylabel) * y.distri2
#    
#    for (i in 1:nvar1) {
#      w <- weight1[fac1 == levels(fac1)[i]]
#      y0 <- y.distri1[i]
#      score0 <- score[fac1 == levels(fac1)[i]]
#      x.moy <- sum(w * score0)
#      x.et <- sqrt(sum(w * (score0 - x.moy)^2))
#      x1 <- x.moy - cet * x.et
#      x2 <- x.moy + cet * x.et
#      etiagauche <- TRUE
#      if ((x1 - xmin) < (xmax - x2)) 
#        etiagauche <- FALSE
#      segments(x1, y0, x2, y0)
#      if (clab1 > 0) {
#        cha <- labels1[i]
#        cex0 <- par("cex") * clab1
#        xh <- strwidth(cha, cex = cex0)
#        xh <- xh + strwidth("x", cex = cex0)
#        yh <- strheight(cha, cex = cex0) * 5/6
#        if (etiagauche) 
#          x0 <- x1 - xh/2
#        else x0 <- x2 + xh/2
#        rect(x0 - xh/2, y0 - yh, x0 + xh/2, y0 + yh, col = "white", border = 1)
#        text(x0, y0, cha, cex = cex0)
#      }
#      points(x.moy, y0, pch = 20, cex = par("cex") * cpoi)
#    }
#    for (i in 1:nvar2) {
#      w <- weight2[fac2 == levels(fac2)[i]]
#      y0 <- y.distri2[i]
#      score0 <- score[fac2 == levels(fac2)[i]]
#      x.moy <- sum(w * score0)
#      x.et <- sqrt(sum(w * (score0 - x.moy)^2))
#      x1 <- x.moy - cet * x.et
#      x2 <- x.moy + cet * x.et
#      etiagauche <- TRUE
#      if ((x1 - xmin) < (xmax - x2)) 
#        etiagauche <- FALSE
#      segments(x1, y0, x2, y0)
#      if (clab2 > 0) {
#        cha <- labels2[i]
#        cex0 <- par("cex") * clab2
#        xh <- strwidth(cha, cex = cex0)
#        xh <- xh + strwidth("x", cex = cex0)
#        yh <- strheight(cha, cex = cex0) * 5/6
#        if (etiagauche) 
#          x0 <- x1 - xh/2
#        else x0 <- x2 + xh/2
#        rect(x0 - xh/2, y0 - yh, x0 + xh/2, y0 + yh, col = "white", border = 1)
#        text(x0, y0, cha, cex = cex0)
#      }
#      points(x.moy, y0, pch = 20, cex = par("cex") * cpoi)
#    }
#  }
#  
#  if (inherits(x, "witwit")) {
#    y <- eval.parent(as.list(x$call)[[2]])
#    oritab <- eval.parent(as.list(y$call)[[2]])
#  } else 
#    oritab <- eval.parent(as.list(x$call)[[2]])
#  
#  l.names <- row.names(oritab)
#  c.names <- names(oritab)
#  oritab <- as.matrix(oritab)
#  a <- x$co[col(oritab), xax]
#  a <- a + x$li[row(oritab), xax]
#  a <- a/sqrt(2 * x$eig[xax] * (1 + sqrt(x$eig[xax])))
#  a <- a[oritab > 0]
#  aco <- col(oritab)[oritab > 0]
#  aco <- factor(aco)
#  levels(aco) <- c.names
#  ali <- row(oritab)[oritab > 0]
#  ali <- factor(ali)
#  levels(ali) <- l.names
#  aw <- oritab[oritab > 0]/sum(oritab)
#  
#  sco.distri.class.2g(a, aco, ali, aw, clab1 = clab.c, clab2 = clab.r, cpoi = cpoi, cet = cet)
#  scatterutil.sub("Rows", csub = csub, possub = "topleft")
#  scatterutil.sub("Columns", csub = csub, possub = "bottomright")
#}



"score.mix" <- function (x, xax = 1, which.var = NULL, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "mix")) 
    stop("Object of class 'mix' expected")
  if(x$nf == 1) 
    xax <- 1
  if((xax < 1) || (xax > x$nf)) 
    stop("non convenient axe number")
  
  ## internal function
  lm.pcaiv <- function(x, df, weights) {
    lm0 <- lm(as.formula(paste("reponse.generic ~ ", paste(names(df), collapse = "+"))), data = cbind.data.frame(x, df), weights = weights)
    return(predict(lm0))
  }
  
  ## data management
  oritab <- as.list(x$call)[[2]]
  evTab <- eval.parent(oritab)
  if(is.null(which.var)) 
    which.var <- 1:length(x$index)
  
  index <- as.character(x$index)
  score <- x$l1[, xax]
  scorecall <- substitute(x$l1[, xax])
  
  ADEglist <- list()
  for (i in which.var) {
    ## parameters management
    sortparameters <- sortparamADEg(...)
    params <- list()
    
    ## data management
    type.var <- index[i]
    col.var <- which(x$assign == i)
    y <- x$tab[, col.var]
    ycall <- substitute(x$tab[, col.var])
    
    ## type of variable : quantitative
    if(type.var == "q") {
      ## parameters management
      params$adepar <- list(psub = list(text = colnames(evTab)[i], position = "topleft"), paxes = list(aspectratio = "fill", draw = TRUE), porigin = list(include = FALSE), pgrid = list(draw = FALSE), plabels = list(cex = 0))
      sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
      
      if(length(col.var) == 1) {
        g1 <- do.call("s.label", c(list(dfxy = call("cbind", scorecall, ycall), plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
        g2 <- xyplot(y ~ score, panel = function(x, y) {panel.abline(lm(y ~ x), lty = 1)})
        g2$call <- call("xyplot", substitute(ycall ~ scorecall), panel = function(x, y) {panel.abline(lm(y ~ x), lty = 1)})
        ADEglist[[i]] <- superpose(g1, g2)
        
      } else {
        ## data management
        lm0 <- lm(as.formula(paste("reponse.generic ~ ", paste(names(y), collapse = "+"))), data = cbind.data.frame(reponse.generic = score, y), weights = rep(1, nrow(y))/nrow(y))
        lm0call <- substitute(lm(as.formula(paste("reponse.generic ~ ", paste(names(ycall), collapse = "+"))), data = cbind.data.frame(reponse.generic = scorecall, ycall), weights = rep(1, nrow(ycall))/nrow(ycall)))
        score.est <- predict(lm0)
        score.estcall <- substitute(predict(lm0call))
        ord0 <- order(y[, 1])
        ord0call <- substitute(order(ycall[, 1]))
        y1call <- call("[", ycall, ord0call, 1)
        x1call <- call("[", score.estcall, ord0call)
        
        ## ADEgS creation
        g1 <- do.call("s.label", c(list(dfxy = call("cbind", scorecall, call("[", ycall, 1:NROW(y), 1)), plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
        g2 <- xyplot(y[ord0, 1] ~ score.est[ord0], panel = function(x, y) {panel.lines(x, y, lty = 1)})
        g2$call <- call("xyplot", substitute(y1call ~ x1call), panel = function(x, y) {panel.lines(x, y, lty = 1)})
        ADEglist[[i]] <- superpose(g1, g2)
      }
    }
    
    ## type of variable : factor
    else if(type.var == "f") {
      ## parameters management
      params$adepar <- list(ppoints = list(pch = "|"), paxes = list(aspectratio = "fill", draw = TRUE), porigin = list(draw = FALSE), pgrid = list(draw = FALSE), psub = list(text = colnames(evTab)[i], position = "topleft"))
      sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
      
      ## data management
      fac <- evTab[, i]
      faccall <- call("[", oritab, 1:NROW(evTab), i)
      meangroup <- call("as.numeric", call("tapply", scorecall, faccall, mean))
      dfxy <- call("cbind", scorecall, call("as.numeric", call("[", meangroup, faccall)))
      
      ## ADEg creation
      g1 <- do.call("s.class", c(list(dfxy = dfxy, fac = faccall, plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
      xlimg1 <- g1@g.args$xlim
      ylimg1 <- g1@g.args$ylim
      g2 <- xyplot(score ~ fac, xlab = "", ylab = "", scales = list(x = list(tck = c(1, 0)), y = list(tck = c(1, 0))), xlim = xlimg1, ylim = ylimg1, 
                   aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(h = as.numeric(tapply(y, x, mean)), a = 0, b = 1, lty = 1)})
      g2$call <- call("xyplot", substitute(scorecall ~ faccall), xlab = "", ylab = "", scales = list(x = list(tck = c(1, 0)), y = list(tck = c(1, 0))), xlim = substitute(xlimg1), 
                      ylim = substitute(ylimg1), aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(h = as.numeric(tapply(y, x, mean)), a = 0, b = 1, lty = 1)})
      ADEglist[[i]] <- superpose(g1, g2)
    }
    
    ## type of variable : ordered
    else if(type.var == "o") {
      ## parameters management
      params$adepar <- list(ppoints = list(pch = 20), paxes = list(aspectratio = "fill", draw = TRUE), porigin = list(draw = FALSE), pgrid = list(draw = FALSE), psub = list(text = colnames(evTab)[i], position = "topleft"))
      sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
      
      ## data management
      lm0 <- lm(as.formula(paste("reponse.generic ~ ", paste(names(y), collapse = "+"))), data = cbind.data.frame(reponse.generic = score, y), weights = rep(1, nrow(y))/nrow(y))
      lm0call <- substitute(lm(as.formula(paste("reponse.generic ~ ", paste(names(ycall), collapse = "+"))), data = cbind.data.frame(reponse.generic = scorecall, ycall), weights = rep(1, nrow(ycall))/nrow(ycall)))
      score.est <- predict(lm0)
      score.estcall <- substitute(predict(lm0call))
      ord0 <- order(y[, 1])
      ord0call <- substitute(order(ycall[, 1]))
      y1call <- call("[", ycall, ord0call, 1)
      x1call <- call("[", score.estcall, ord0call)
      
      ## ADEgS creation
      g1 <- do.call("s.label", c(list(dfxy = call("cbind", scorecall, call("[", ycall, 1:NROW(y), 1)), plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
      g2 <- xyplot(y[ord0, 1] ~ score.est[ord0], panel = function(x, y) {panel.lines(x, y)})
      g2$call <- call("xyplot", substitute(y1call ~ x1call), panel = function(x, y) {panel.lines(x, y)})
      ADEglist[[i]] <- superpose(g1, g2)
    }
  }
  ADEglist <- ADEglist[which.var]
  
  ## ADEgS creation
  posmatrix <- layout2position(.n2mfrow(length(which.var)), ng = length(which.var))
  object <- new(Class = "ADEgS", ADEglist = ADEglist, positions = posmatrix, add = matrix(0, ncol = length(which.var), nrow = length(which.var)), Call = match.call())
  names(object) <- colnames(evTab)[which.var]
  object@Call <- match.call()
  if(plot) 
    print(object)
  invisible(object)
}


"score.pca" <- function (x, xax = 1, which.var = NULL, pos = -1, storeData = TRUE, plot = TRUE, ...) {
  if(!inherits(x, "pca")) 
    stop("Object of class 'pca' expected")
  if(x$nf == 1) 
    xax <- 1
  if((xax < 1) || (xax > x$nf))
    stop("non convenient axe number")
  
  ## prepare
  oritab <- as.list(x$call)[[2]]
  evTab <- eval.parent(oritab)
  if(is.null(which.var))
    which.var <- 1:ncol(evTab)
  
  ## parameter management
  sortparameters <- sortparamADEg(...)
  params <- list()
  params$adepar <- list(paxes = list(aspectratio = "fill", draw = TRUE), porigin = list(include = FALSE), pgrid = list(draw = FALSE), plabels = list(cex = 0))
  sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
  
  ## creation of each individual ADEg
  ADEglist <- list()
  for(i in which.var) {
    dfxy <- call("cbind", substitute(x$l1[, xax]), call("[", oritab, 1:NROW(evTab), i))
    
    g1 <- do.call("s.label", c(list(dfxy = dfxy, plot = FALSE, storeData = storeData, pos = pos - 2), c(sortparameters$adepar, list(psub.text = colnames(evTab)[i])), sortparameters$trellis, sortparameters$g.args, sortparameters$rest))
    g2 <- xyplot(eval(dfxy)[, 2] ~ eval(dfxy)[, 1], aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(lm(y ~ x))})
    g2$call <- call("xyplot", substitute(dfxy[, 2] ~ dfxy[, 1]), aspect = g1@adeg.par$paxes$aspectratio, panel = function(x, y) {panel.abline(lm(y ~ x))})
    ADEglist[[i]] <- superpose(g1, g2)
  }
  ADEglist <- ADEglist[which.var]
  
  ## ADEgS creation
  posmatrix <- layout2position(.n2mfrow(length(which.var)), ng = length(which.var))
  object <- new(Class = "ADEgS", ADEglist = ADEglist, positions = posmatrix, add = matrix(0, ncol = length(which.var), nrow = length(which.var)), Call = match.call())
  names(object) <- colnames(evTab)[which.var]
  object@Call <- match.call()
  if(plot) 
    print(object)
  invisible(object)
}
