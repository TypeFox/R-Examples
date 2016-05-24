
# updated by Raimon in August 2007:
gsi.ilrBase2signary = function(V){
  if(is.null(dim(V))){ dim(V)=c(length(V),1)  }
  signary = sign(V)
  Vzero = abs(V)>1.e-15
  signary = signary * Vzero
  if(ncol(signary)>1){
    if(!any(geometricmeanRow(signary)==0)){
      stop("this ilr base matrix does not correspond to a partition!")
    }
  }
  return(signary)
}


# re-order a compositional data set and its ilr coordinates to ease its plotting
gsi.OrderIlr = function(V){
 # first, coordinates are ordered in increasing number of zeroes
 Vzero = abs(V)>1.e-15
 aux = rep(1,nrow(V)) %*% Vzero
 aux = order(aux, decreasing=TRUE)
 V = V[,aux]
 # then compute the ordering for the parts so that non-zero parts are always together
 signary = gsi.ilrBase2signary(V)
 signary = data.frame(signary)
 aux = do.call("order",as.list(signary))
 output = list(ilrBase=V,order=aux) # return the matrix (with coordinates reordered, but not parts)
 return(output)
}


# updated by Raimon in April 2008
CoDaDendrogram = function (X, V = NULL, expr=NULL, mergetree = NULL, signary = NULL, 
    range = c(-4,4), ..., xlim = NULL, ylim = NULL, yaxt = NULL, box.pos = 0,
    box.space = 0.25, col.tree = "black", lty.tree = 1, lwd.tree = 1,
    col.leaf = "black", lty.leaf = 1, lwd.leaf = 1, add = FALSE,border=NULL,
    type = "boxplot"){
    if (is.na(match(class(X), c("acomp", "rcomp")))) {
        stop("CoDaDendrogram only valid for relative compositions, of class acomp or rcomp")
    }
    if (!add) {
        if (is.null(V) & is.null(expr) & is.null(mergetree) & is.null(signary)) {
            stop("a hierarchical basis is needed! give one and only one of mergetree, signary, expr or V")
        }
        if (!is.null(expr)) {
            V = balanceBase(X, expr)
        }
        if (!is.null(mergetree)) {
            V = gsi.buildilrBase(gsi.merge2signary(mergetree))
        }
        if (!is.null(signary)) {
            V = gsi.buildilrBase(signary)
        }
        Vo = gsi.OrderIlr(V = V)
        idtx = idt(X, V = Vo$ilrBase)
         if(is.null(range)){ range=range(idtx) }
        varx = diag(var(idtx))
        meanx = mean(idtx)
        signary = gsi.ilrBase2signary(Vo$ilrBase)
        binary = abs(signary)
        heights = varx * 0
        heights[1:length(heights)] = max(binary[, 1:ncol(binary)] %*%
            varx[1:ncol(binary)])
        maxheight = max(heights)
        for (i in 1:ncol(binary)) {
            heights[i:ncol(binary)] = heights[i:ncol(binary)] -
                varx[i] * sign(binary[, i] %*% binary[, i:ncol(binary)])
        }
        nodes = matrix(0, ncol = 2, nrow = ncol(binary))
        is.leaf = nodes
        for (i in ncol(binary):1) {
            nparts = c(sum(signary[, i] < 0), sum(signary[, i] >
                0))
            for (j in 1:2) {
                if (nparts[j] == 1) {
                  take = (signary[, i] == (-1)^j)
                  take = (c(1:length(Vo$order))[take] == Vo$order)
                  nodes[i, j] = c(1:length(Vo$order))[take]
                  is.leaf[i, j] = 1
                }
                else {
                  take = (signary[, i] == (-1)^j)
                  take = as.integer(take) * binary[, i]
                  take = as.logical(take)
                  aux = rep(1, sum(take)) %*% binary[take, ]
                  aux[1:i] = 0
                  i1 = c(1:length(aux))[aux == max(aux)]
                  nodes[i, j] = approx(x = range, y = nodes[i1,
                    ], xout = meanx[i1], rule = 2)$y
                }
            }
        }
        if (is.null(xlim)) {
            xlim = c(0, 1 + nrow(V))
        }
        if (is.null(ylim)) {
            ylim = c(0, maxheight)
        }
        plot.window(xlim = xlim, ylim = ylim)
        plot(x = xlim, y = ylim, ann = FALSE, bty = "n", xaxt = "n",
            yaxt = yaxt, col = "white")
        axis(side = 1, at = 1:ncol(X), labels = colnames(X)[Vo$order],
            las = 2, lty = 0)
#        segments(x0 = nodes[ncol(V), 1], y0 = 0, x1 = nodes[ncol(V),
#            2], y1 = 0, col = col.tree, lty = lty.tree, lwd = lwd.tree)
#        for (i in 1:(ncol(V) - 1)) {
        for (i in 1:(ncol(V))) {
            segments(x0 = nodes[i, 1], y0 = heights[i], x1 = nodes[i,
                2], y1 = heights[i], col = col.tree, lty = lty.tree,
                lwd = lwd.tree)
        }
        for (i in 1:nrow(is.leaf)) {
            for (j in 1:ncol(is.leaf)) {
                if (is.leaf[i, j]) {
                  segments(x0 = nodes[i, j], y0 = heights[i],
                    x1 = nodes[i, j], y1 = 0, lty = lty.leaf,
                    col = col.leaf, lwd = lwd.leaf)
                }
            }
        }
        gsi.setCoorInfo(basis = Vo$ilrBase, nodes = nodes, heights = heights,
            range = range)
    }
    aux = gsi.getCoorInfo()
    basis = aux$basis
    nodes = aux$nodes
    heights = aux$heights
    range = aux$range
    idtx = idt(X, V = basis)
    methods = c("boxplot", "density", "histogram", "lines", "nothing",
        "points")
    what = methods[pmatch(type, methods)]
    if (what == "lines") {
        varx = diag(var(idtx))
        meanx = mean(idtx)
        for (i in 1:length(varx)) {
            aux = approx(x = range, y = nodes[i, ], xout = meanx[i],
                rule = 2)$y
            segments(x0 = aux, y0 = heights[i], x1 = aux, y1 = heights[i] +
                varx[i], ...)
        }
    }
    if (what == "boxplot") {
        varx = diag(var(idtx))
        meanx = mean(idtx)
        minheight = min(varx)
        for (i in 1:ncol(idtx)) {
            aux = approxfun(x = range, y = nodes[i, ], rule = 2)
            segments(x0 = aux(meanx[i]), y0 = heights[i], x1 = aux(meanx[i]),
                y1 = heights[i] + varx[i], ...)
            qx = quantile(idtx[, i], probs = seq(0, 1, 0.25))
            mbx = box.space * minheight
            hx = heights[i] + mbx * (box.pos - 1)/2
            rect(xleft = aux(qx[2]), ybottom = hx, xright = aux(qx[3]),
                ytop = hx + mbx, ...,border=border)
            rect(xleft = aux(qx[3]), ybottom = hx, xright = aux(qx[4]),
                ytop = hx + mbx, ...,border=border)
            hx = heights[i] + mbx/2 * (box.pos - 1)/2
            rect(xleft = aux(qx[1]), ybottom = hx, xright = aux(qx[2]),
                ytop = hx + mbx/2, ...,border=border)
            rect(xleft = aux(qx[4]), ybottom = hx, xright = aux(qx[5]),
                ytop = hx + mbx/2, ...,border=border)
        }
    }
    if (what == "points") {
        for (i in 1:ncol(idtx)) {
            auxfun = approxfun(x = range, y = nodes[i, ], rule = 2)
            aux = auxfun(idtx[, i])
            points(x = aux, y = rep(heights[i], length(aux)),
                ...)
        }
    }
    if (what == "histogram") {
        stop("sorry, type=histogram still unimplemented")
    }
    if (what == "density") {
        stop("sorry, type=density still unimplemented")
    }
    replot(plot=match.call(),add=FALSE)  
  }
