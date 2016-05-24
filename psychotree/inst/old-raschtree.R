## high-level convenience interface
raschtree <- function(formula, data, minsplit = 10, reltol = 1e-10,
  deriv = c("sum", "diff", "numeric"), hessian = TRUE, maxit = 100L, ...)
{
  ## transform formula
  stopifnot(length(formula) > 2)
  formula <-  formula(terms(formula, data = data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]]

  ## formula/data/model pre-processing
  raschmod <- RaschModel(reltol = reltol, deriv = deriv, hessian = hessian, maxit = maxit)
  ff <- attr(ParseFormula(ff), "formula")
  ff$input[[3]] <- ff$input[[2]]
  ff$input[[2]] <- ff$response[[2]]
  ff <- dpp(raschmod, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), 
    data = data, na.action = na.pass)

  ## data sanity checks
  y <- as.matrix(ff@get("response"))
  if(ncol(y) < 3) stop("need at least three items")
  if(!all(as.vector(y) %in% c(0:1, NA))) stop("y must be a binary 0/1 matrix (potentially with NAs)")
  if(!all(apply(y, 1, function(x) all(0:1 %in% x))))
    stop("each row of y must have at least one 0 and one 1 entry")

  ## call mob()
  rval <- mob(ff, model = raschmod, control = mob_control(minsplit = minsplit,
      objfun = function(object) - as.vector(logLik(object)), ...))

  ## add class and return
  structure(list(mob = rval), class = "raschtree")
}

## convenience plotting
plot.raschtree <- function(x, terminal_panel = node_raschplot, tnex = 2,
  pval = TRUE, id = TRUE, ...) {
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex,
    tp_args = list(id = id, ...), ip_args = list(pval = pval, id = id))
}

## hand-crafted "Next()" to bridge to
## un-exported S4 classes "mob"/"BinaryTree", argh!
logLik.raschtree <- function(object, ...) logLik(object$mob, ...)
sctest.raschtree <- function(x, ...) sctest(x$mob, ...)
weights.raschtree <- function(object, ...) weights(object$mob, ...)
summary.raschtree <- function(object, ...) summary(object$mob, ...)
print.raschtree <- function(x, ...) {
  print(x$mob, ...)
  invisible(x)
}

## parameters for Rasch trees
coef.raschtree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  nam <- names(object@tree$model$coefficients)
  rval <- sapply(nodes(object, node), function(z) coef(z$model, ...)[nam])
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
    colnames(rval) <- nam
  }
  return(rval)
}


itempar.raschtree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) itempar.RaschModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

threshold.raschtree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) threshold.RaschModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

worth.raschtree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- sapply(nodes(object, node), function(z) worth(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}

## visualization function
node_raschplot <- function(mobobj, id = TRUE, difficulty = TRUE,
  center = TRUE, index = TRUE, names = NULL, abbreviate = FALSE, ref = TRUE,
  col = cbind("lightgray", "black"), refcol = "lightgray", linecol = "black",
  lty = 2, cex = 0.5, pch = cbind(19, 1), xscale = NULL, yscale = NULL,
  xaxis = TRUE, yaxis = TRUE, ylines = 1.5)
{
    ## extract parameter of interest
    node <- 1:max(terminal_nodeIDs(mobobj@tree))
    cf <- t(sapply(nodes(mobobj, node), function(z) if(center) {
        worth(z$model, difficulty = difficulty)
      } else {
        worth(z$model, difficulty = difficulty) - worth(z$model, difficulty = difficulty)[1]
      }))
    rownames(cf) <- node
    ncf <- NCOL(cf)

    ## labeling
    if(is.null(names)) names <- !index
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }
    if(!names & index) {
      lab <- rep("", ncf)
      lab[c(1, ncf)] <- c(1, ncf)
      colnames(cf) <- lab
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)

    ## graphical parameter processing  
    if(NCOL(pch) == 2) {
      pch2 <- pch[,2]
      pch <- pch[,1]
    } else {
      pch2 <- NA
    }
    if(NCOL(col) == 2) {
      col2 <- col[,2]
      col <- col[,1]
    } else {
      col2 <- NULL
    }
    pch <- rep(pch, length.out = ncf)
    col <- rep(col, length.out = ncf)
    cex <- rep(cex, length.out = ncf)
    pch2 <- rep(pch2, length.out = ncf)
    if(!is.null(col2)) col2 <- rep(col2, length.out = ncf)
  
    
    if(index) {
      x <- 1:NCOL(cf)
      if(is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    } else {
      x <- rep(0, length(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)      
    }
    yrange <- range(cf[is.finite(cf)], na.rm = TRUE)
    if(is.null(yscale)) yscale <- yrange + c(-0.1 - 0.2 * any(cf <= -Inf), 0.1 + 0.2 * any(cf >= Inf)) * diff(yrange)
         
    ## panel function for bt plots in nodes
    rval <- function(node) {
    
        ## dependent variable setup
	cfi <- cf[node$nodeID,]
	cf_ref <- mean(cfi)        

        cf_ident <- is.finite(cfi) & !is.na(cfi)
        cf_inf <- cfi >= Inf
        cf_ninf <- cfi <= -Inf
        cfi[is.na(cfi)] <- cf_ref
        if(index) {
          cfi[cf_ninf] <- yscale[1]
          cfi[cf_inf] <- yscale[2]
        }

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
			   heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_raschplot", node$nodeID, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()

        ## actual plot	
        plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
	    xscale = xscale, yscale = yscale, 
	    name = paste("node_raschplot", node$nodeID, "plot", sep = ""))
        pushViewport(plot_vpi)
	
        grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = refcol), default.units = "native")
	if(index) {
	  grid.lines(x, cfi, gp = gpar(col = linecol, lty = lty), default.units = "native")
  	  grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = ifelse(cf_ident, pch, NA), default.units = "native")
          if(!is.null(col2)) grid.points(x, cfi, gp = gpar(col = col2, cex = cex), pch = ifelse(cf_ident, pch2, NA), default.units = "native")
	  if(xaxis) grid.xaxis(at = x, label = colnames(cf))
	} else {	  
  	  if(names) grid.text(names(cfi), x = x, y = cfi, default.units = "native")
	    else grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
	}
        if(yaxis) grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
        grid.rect(gp = gpar(fill = "transparent"))

	upViewport(2)
    }
	    
    return(rval)
}
class(node_raschplot) <- "grapcon_generator"

