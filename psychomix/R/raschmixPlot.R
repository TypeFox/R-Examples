setMethod("plot", signature(x = "raschmix", y = "missing"),
          function(x, y, component = NULL, difficulty = TRUE,
                   center = TRUE, index = TRUE, names = TRUE,
                   abbreviate = FALSE, ref = TRUE, col =  "black",
                   refcol = "lightgray", linecol = NULL, lty = 2, cex = 1,
                   pch = 19, type = NULL, ylim = NULL, xlab = "Items",
                   ylab = NULL, legend = TRUE, pos = "topright",
                   srt = 45, adj = c(1.1, 1.1), ...){
  ##  modified code from psychotree
         
  ## parameters to be plotted
  if (is.null(component)) component <- 1:x@k
  cf <- worth(x, difficulty = difficulty)
  cf <- cf[ ,component, drop = FALSE]
  cf_ident <- is.finite(cf)
  cf_inf <- cf >= Inf
  cf_ninf <- cf <= -Inf
  if(!center){
      ref <- sapply(component, function(x){
          ref.i <- min(which(cf_ident[,x]))
          cf[ref.i, x]
          })
      cf <- cf - matrix(rep(ref, each = nrow(cf)), ncol = ncol(cf))
  }
  cf_ref <- 0
  ncf <- nrow(cf)

  ## labeling
  if(is.null(names)) names <- !index 
  if(isTRUE(names) & length(unique(sapply(strsplit(rownames(cf), split = ".", fixed = TRUE), function(x) x[1]))) == 1){
    rownames(cf) <- sapply(strsplit(rownames(cf), split = ".", fixed = TRUE),
                           function(x) paste(x[-1], collapse = "."))
  }
  if(is.character(names)) {
    rownames(cf) <- names
    names <- TRUE
  }
  if(!names & index) {
    lab <- rep(NA, ncf)
    lab[c(1, ncf)] <- c(1, ncf)
    pr <- pretty(1:ncf)
    pr <- pr[pr > 1 & pr < ncf]
    lab[pr] <- pr    
    rownames(cf) <- lab
  }

  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(rownames(cf), type = "width"))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  rownames(cf) <- abbreviate(rownames(cf), abbreviate)

  ## graphical parameter processing  
  if(is.null(type)) type <- if(index) "b" else "p"
  if (is.vector(pch)){
    pch <- matrix(rep(rep(pch, length.out = ncol(cf)), each = nrow(cf)),
                  ncol = ncol(cf))
  } else if (is.matrix(pch)){
    if (!all(is.logical(all.equal(ncol(pch), ncol(cf))),
             is.logical(all.equal(nrow(pch), nrow(cf)))))
      stop("pch needs to be a either a vector or a matrix of dimesion number of item parameters x number of components")
  }
  if ((length(component) > 1) & (length(col) == 1))
    col <- qualitative_hcl(length(component))
  if (is.vector(col)){
    col <- matrix(rep(rep(col, length.out = ncol(cf)), each = nrow(cf)),
                  ncol = ncol(cf))
  } else if (is.matrix(col)){
    if (!all(is.logical(all.equal(ncol(col), ncol(cf))),
             is.logical(all.equal(nrow(col), nrow(cf)))))
      stop("col needs to be a either a vector or a matrix of dimesion number of item parameters x number of components")
  }
  if (is.vector(cex)){
    cex <- matrix(rep(rep(cex, length.out = ncol(cf)), each = nrow(cf)),
                  ncol = ncol(cf))
  } else if (is.matrix(cex)){
    if (!all(is.logical(all.equal(ncol(cex), ncol(cf))),
             is.logical(all.equal(nrow(cex), nrow(cf)))))
      stop("cex needs to be a either a vector or a matrix of dimesion number of item parameters x number of components")
  }
  pch[!cf_ident] <- NA
  if (is.null(linecol)) {linecol <- col[1, ]} else {
    linecol <- rep(linecol, length.out = ncol(cf))
  }
  lty <- rep(lty, length.out = ncol(cf))

  if(is.null(ylim)) ylim <- range(cf[cf_ident])
  ylim <- rep(ylim, length.out = 2)
  if(any(!is.finite(cf))) {
    ydiff <- diff(ylim) * 0.7
    if(index & any(cf_ninf)) ylim[1] <- ylim[1] - ydiff
    if(index & any(cf_inf))  ylim[2] <- ylim[2] + ydiff
  }

  ## substitute non-identified parameters with plottable values
  cf[is.na(cf)] <- cf_ref
  if(index) {
    cf[cf_ninf] <- ylim[1]
    cf[cf_inf] <- ylim[2]
  }

  if(is.null(ylab)) ylab <- paste(if(center) "Centered item" else "Item",
    if(difficulty) "difficulty" else "easiness", "parameters")

  ## raw plot
  ix <- seq(along = cf[,1]) 
  ix.ni <- if (legend){
    matrix(rep(seq(from = 2, to = nrow(cf)*(1 - 1/ncol(cf)),
                   length.out = ncol(cf)), each = ncf),
           ncol = ncol(cf))
  } else {
    matrix(rep(seq(from = 2, to = nrow(cf)-1, length.out = ncol(cf)),
               each = ncf), ncol = ncol(cf))
  }

  plot(ix, cf[,1], xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ylim = ylim, ...)
  if(ref) abline(h = cf_ref, col = refcol)
  axis(2)
  box()  

  ## actual data
  if(!index & names) {
    for (i in 1:ncol(cf)){
      text(rownames(cf), x = ix.ni[,i], y = cf[,i], col = col[,i], ...)
    }
    if(legend) {
      if (any(!cf_ident)){
        ni_items <- unlist(lapply(component, function(comp){
          if(any(!cf_ident[,comp])){
            paste("Comp. ", comp, ": ", paste(
              rownames(cf)[!cf_ident[,comp]], collapse = ", "), sep = "")
            } 
          }))
        legend(pos, c(paste("Comp.", component),
                      "Not identified:", ni_items),
               fill = c(col[1,], rep(NA, 1 + length(ni_items))),
               border = rep(NA, ncol(cf) + 1 + length(ni_items)), bty = "n")
      } else {
        legend(pos, paste("Comp.", component), fill = col[1,], bty = "n")
      }
    }
  } else {
    for (i in 1:ncol(cf)){
      if(type %in% c("l", "b", "o")) lines(ix, cf[,i], type = type,
                                           lty = lty[i], pch = NA,
                                           col = linecol[i])
      lines(ix, cf[,i], type = type, lty = 0, pch = pch[,i], col = col[,i],
            cex = cex[,i])
      if(type %in% c("l", "b", "o")) {
        if(any(cf_ninf[,i])) for(j in which(cf_ninf[,i])){
          lines(c(ix[j], ix[j]), c(ylim[1], ylim[1] - 10 * ydiff), type = "l",
                lty = lty[i])
        }
        if(any(cf_inf[,i])) for(j in which(cf_inf[,i])){
          lines(c(ix[j], ix[j]), c(ylim[2], ylim[2] + 10 * ydiff), type = "l",
                lty = lty[i])
        }
      }
    }
    if(index){
      if(names){
        text(ix, par("usr")[3], labels = rownames(cf), srt = srt, adj = adj, xpd = TRUE, cex = 0.9)
      } else {
        axis(1, at = ix, labels = rownames(cf))
      }
    }
    if (length(component) > 1 & legend){
      pch <- sapply(component, function(x){
          ref.i <- min(which(cf_ident[,x]))
          pch[ref.i, x]
      })
      legend(pos, legend = paste("Comp.", component), bty = "n",
             col = col[1,], pch = pch, lty = lty)
      ## <FIXME> sensible transformation of "plot-cex" to "legend-cex"? </FIXME>
    }
  }
})

## convenience function for choosing (slightly tweaked) qualitative HCL palette
qualitative_hcl <- function(n, c. = 80, l. = 60, alpha = 1, fixup = TRUE)  
{
  hclcol <- function(n) hcl(seq(0, 360 * (n - 1)/n, length = n),
    c = c., l = l., alpha = alpha, fixup = fixup)
  rval <- if(n <= 1) {
    hcl(0, 0, 0, alpha = alpha, fixup = fixup)
  } else if(n <= 6) {
    hclcol(12)[c(9, 1, 5, 3, 11, 7)[1:n]]
  } else {
    hclcol(n)
  }
  if(missing(alpha)) {
    structure(substr(rval, 1L, 7L), names = names(rval))
  } else {
    return(rval)
  }
}

## histogram method reusing flexmix's plot() method
histogram.raschmix <- function(x, data, root = TRUE, ...) {
  if(!missing(data)) warning("argument 'data' is ignored")
  ## just calling plot() within an S3 method does not seem to do the right dispatch
  getMethod("plot", c(x = "flexmix", y = "missing"))(as(x, "flexmix"), root = root, ...)
}

## xyplot method
xyplot.raschmix <- function(x, data,
  component = NULL, item = NULL,
  difficulty = TRUE, plot.type = c("multiple", "single"),
  auto.key = NULL, type = "b", lty = NULL, xlab = "Item", ylab = NULL,
  panel = NULL, scales = NULL, ...)
{
  ## process data
  y <- worth(x, difficulty = difficulty)
  if(!missing(data)) warning("'data' argument is ignored")

  ## select items/components
  if(is.null(component)) component <- 1:NCOL(y)
  if(is.null(item)) item <- 1:NROW(y)
  if(is.character(item)) item <- which(item %in% rownames(y))
  y <- y[item, component, drop = FALSE]

  ## set up auxiliary data.frame
  d <- data.frame(
    y = as.vector(y),
    x = rep(item, length(component)),
    z = factor(rep(component, each = length(item)),
      levels = component, labels = paste("Comp.", component))
  )
  
  ## graphical arguments
  plot.type <- match.arg(plot.type)
  if(plot.type == "single") {
    f <- y ~ x
    groups <- ~ z
    if(is.null(lty)) lty <- trellis.par.get("superpose.line")$lty
  } else {
    f <- y ~ x | z
    groups <- NULL
    if(is.null(lty)) lty <- trellis.par.get("superpose.line")$lty[2]
  }
  if(is.null(ylab)) ylab <- paste("Centered item",
    if(difficulty) "difficulty" else "easiness",
    "parameters")
  if(is.null(auto.key)) auto.key <- plot.type == "single"
  if(is.null(scales)) scales <- list(x = list(at = item, alternating = 1))
  if(is.null(panel)) panel <- function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.abline(h = 0, reference = TRUE)
  }

  ## call xyplot() formula method
  xyplot(f, groups = groups, data = d,
    type = type, lty = lty, xlab = xlab, ylab = ylab,
    auto.key = auto.key, scales = scales, panel = panel, ...)
}
