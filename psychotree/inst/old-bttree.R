## infrastructure (copied from party)
terminal_nodeIDs <- function(node) {
  if(node$terminal) return(node$nodeID)
  ll <- terminal_nodeIDs(node$left)
  rr <- terminal_nodeIDs(node$right)
  return(c(ll, rr))
}

## high-level convenience interface
bttree <- function(formula, data, na.action = na.pass,
  type = "loglin", ref = NULL, undecided = NULL, position = NULL, minsplit = 10, ...)
{
  ## transform formula
  stopifnot(length(formula) > 2)
  formula <-  formula(terms(formula, data = data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]]

  ## call mob()
  rval <- mob(ff, data = data, model = btReg(type = type, ref = ref,
      undecided = undecided, position = position),
    control = mob_control(minsplit = minsplit, ...), na.action = na.action)

  ## add class and return
  structure(list(mob = rval), class = "bttree")
}

## convenience plotting
plot.bttree <- function(x, terminal_panel = node_btplot, tnex = 2, ...) {
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex, tp_args = list(...))
}

## hand-crafted "Next()" to bridge to
## un-exported S4 classes "mob"/"BinaryTree", argh!
deviance.bttree <- function(object, ...) deviance(object$mob, ...)
logLik.bttree <- function(object, ...) logLik(object$mob, ...)
sctest.bttree <- function(x, ...) sctest(x$mob, ...)
weights.bttree <- function(object, ...) weights(object$mob, ...)
summary.bttree <- function(object, ...) summary(object$mob, ...)
print.bttree <- function(x, ...) {
  print(x$mob, ...)
  invisible(x)
}

## parameters for BT trees
coef.bttree <- function (object, node = NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- sapply(nodes(object, node), function(z) coef(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
  }
  return(rval)
}

worth.bttree <- function (object, node = NULL, ...) 
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

predict.bttree <- function(object, newdata = NULL,
  type = c("worth", "rank", "node"), ...)
{
  type <- match.arg(type)
  
  ## get nodes
  nodes <- predict(object$mob, newdata = newdata, type = "node", ...)
  if(type == "node") return(nodes)
  
  ## get worth
  w <- worth(object)
  w <- w[as.character(nodes), , drop = FALSE]
  rownames(w) <- NULL
  if(type == "worth") return(w)
  
  ## get order
  o <- t(apply(-w, 1, rank))
  return(o)
}


## visualization function
node_btplot <- function(mobobj, id = TRUE,
  worth = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## extract parameter of interest
    node <- 1:max(terminal_nodeIDs(mobobj@tree))
    cf <- t(sapply(nodes(mobobj, node), function(z)
      if(worth) worth(z$model) else coef(z$model, all = FALSE, ref = TRUE)))
    rownames(cf) <- node

    if(!worth) {
      if(is.character(ref) | is.numeric(ref)) {
        reflab <- ref
        ref <- TRUE
      } else {
        reflab <- mobobj@tree$model$ref
      }
      if(is.character(reflab)) reflab <- match(reflab, mobobj@tree$model$labels)
      cf <- cf - cf[,reflab]
    }

    ## reference
    if(worth) {
      cf_ref <- 1/ncol(cf)
    } else {
      cf_ref <- 0
    }

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)
    
    if(index) {
      x <- 1:NCOL(cf)
      if(is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    } else {
      x <- rep(0, length(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)      
    }
    if(is.null(yscale)) yscale <- range(cf) + c(-0.1, 0.1) * diff(range(cf))
         
    ## panel function for bt plots in nodes
    rval <- function(node) {
    
        ## dependent variable setup
	cfi <- cf[node$nodeID,]

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
			   heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_btplot", node$nodeID, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()

        ## actual plot	
        plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
	    xscale = xscale, yscale = yscale, 
	    name = paste("node_btplot", node$nodeID, "plot", sep = ""))
        pushViewport(plot_vpi)
	
        grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = linecol), default.units = "native")
	if(index) {
	  grid.lines(x, cfi, gp = gpar(col = col, lty = 2), default.units = "native")
	  grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
	  grid.xaxis(at = x, label = if(names) names(cfi) else x)
	} else {	  
  	  if(names) grid.text(names(cfi), x = x, y = cfi, default.units = "native")
	    else grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
	}
        grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
        grid.rect(gp = gpar(fill = "transparent"))

	upViewport(2)
    }
	    
    return(rval)
}
class(node_btplot) <- "grapcon_generator"

