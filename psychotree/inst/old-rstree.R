### high-level convenience interface for creating rstrees
rstree <- function (formula, data, minsplit = 20, reltol = 1e-10,
                    deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L, ...)
{
  ## transform formula
  stopifnot(length(formula) > 2)
  formula <-  formula(terms(formula, data=data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]]

  ## formula/data/model pre-processing
  rsmod <- RSModel(reltol = reltol, deriv = deriv, hessian = hessian, maxit = maxit)
  ff <- attr(ParseFormula(ff), "formula")
  ff$input[[3]] <- ff$input[[2]]
  ff$input[[2]] <- ff$response[[2]]
  ff <- dpp(rsmod, as.formula(ff$input), other=list(part=as.formula(ff$blocks)), 
            data=data, na.action=na.pass)

  ## data sanity checks
  y <- as.matrix(ff@get("response"))
  if(ncol(y) < 3) stop("need at least three items")
  rmax <- sum(apply(y, 2, max, na.rm = TRUE))
  if(any(apply(y, 1, function(x) sum(x, na.rm = TRUE)) %in% c(0, rmax)))
    stop("Please remove subjects who only scored zero or the highest category.")

  ## call mob()
  res <- mob(ff, model = rsmod, control=mob_control(minsplit = minsplit,
                                 objfun = function (object) -as.vector(logLik(object)), ...))

  ## add class and return
  structure(list(mob = res), class = c("rstree", "pctree"))
}

logLik.rstree <- function (object, ...) logLik(object$mob, ...)

sctest.rstree <- function (x, ...) sctest(x$mob, ...)

weights.rstree <- function (object, ...) weights(object$mob, ...)

summary.rstree <- function (object, node=NULL, ...) summary(object$mob, ...)

print.rstree <- function (x, ...) { 
  print(x$mob, ...)
  invisible(x)
}

plot.rstree <- function (x, terminal_panel = node_effects, tnex = 2, ...)
{
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex, tp_args = list(...))
}

coef.rstree <- function (object, node=NULL, ...) 
{
  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  nam <- names(object@tree$model$coefficients)
  rval <- sapply(nodes(object, node), function (z) coef(z$model, ...))
  if (!is.null(dim(rval))) {
    rval <- t(rval)
    rownames(rval) <- node
    colnames(rval) <- nam
  }
  return (rval)
}

itempar.rstree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) itempar.RSModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

threshold.rstree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) threshold.RSModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

node_rsmplot <- function(mobobj, names = NULL, ref = NULL,
                         ylab = "Latent trait", ylim = NULL,
                         col = c("gray", "gray"), pch = c(21, 22), lty = c(2, 3),
                         refline = TRUE, reflinecol = "lightgray")
{
  ## check input
  stopifnot(!is.null(mobobj))
  stopifnot(length(col) == 2)
  stopifnot(length(pch) == 2)
  stopifnot(length(lty) == 2)

  ## setup item parameters from terminal nodes
  tnodes <- terminal_nodeIDs(mobobj@tree)
  nodes_lst <- nodes(mobobj, tnodes)
  fun <- if (inherits(nodes_lst[[1]]$model, "RaschModel")) threshold.RaschModel else if (inherits(nodes_lst[[1]]$model, "RSModel")) threshold.RSModel else threshold.PCModel
  cf_lst <- lapply(nodes_lst, function (node) fun(node$model, ref = ref, simplify = FALSE, vcov = FALSE))
  ipar_lst <- lapply(cf_lst, "[[", 1)
  cpar_lst <- lapply(cf_lst, "[[", 2)
  names(ipar_lst) <- names(cpar_lst) <- tnodes
  
  ## maximum number of items and categories
  m_max <- max(sapply(ipar_lst, length))
  o_max <- max(sapply(cpar_lst, length))

  ## setup y-axis and labels
  if (is.null(ylim)) ylim <- extendrange(c(unlist(ipar_lst, use.names = FALSE), unlist(cpar_lst, use.names = FALSE)), f = 0.25)
  if (is.null(names)) names <- c(paste("I", 1:m_max, sep = ""), paste("T", 1:o_max, sep = "")) else stopifnot(length(names) == m_max + o_max)

  ## finally, the panel function.
  panelfun <- function(node) {
    
    ## get id, parameters and basic data
    nid <- as.character(node$nodeID)
    ipar <- ipar_lst[[nid]]
    cpar <- cpar_lst[[nid]]
    m <- length(ipar)
    o <- length(cpar)
    
    ## setup x positions, xlim and label
    ix <- 1:(m + o)
    xlim <- c(0, ix[m + o] + 1)
    lab <- paste("node", nid, sep = "")

    ## terminal panel viewport setup
    top.vp <- viewport(layout = grid.layout(nrow = 2, ncol = 1, widths = unit(1, "null"), heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"), name = paste(lab, "_rsmplot", sep = ""))
    pushViewport(top.vp)
    grid.rect(gp = gpar(fill = "white", col = 0), name = paste(lab, "_border", sep = ""))

    ## main title
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1, name = paste(lab, "_title_vp", sep = "")))
    grid.text(paste("Node ", nid, " (n = ", sum(node$weights), ")", sep = ""), name = paste(lab, "_title", sep = ""))
    upViewport()
    
    ## finally the actual plot
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2, name = lab))
    lab <- paste(lab, "_plot", sep = "")

    ## setup plotting area (3x3)
    wcol <- if (is.null(ylab)) c(2.5, 1, 1) else c(4, 1, 1)
    hrow <- c(0.5, 1, 1)

    top.vp <- viewport(layout = grid.layout(nrow = 3, ncol = 3, widths = unit(wcol, c("lines", "null", "lines")),
                         heights = unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top_vp", sep = ""))
    bmargin.vp <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
    lmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
    rmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
    plot.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
    pushViewport(top.vp)

    ## setup plotting window and (if requested) reference line
    pushViewport(plot.vp)
    if (refline) grid.lines(x = c(0, ix[m] + 0.49), y = rep(mean(ipar), 2), default.units = "native",
                        name = paste(lab, "_reference_line", sep = ""), gp = gpar(col = reflinecol))

    ## plot item parameters
    grid.lines(x = ix[1:m], y = ipar, default.units = "native", gp = gpar(lty = lty[1]), name = paste(lab, "_ipar_line", sep = ""))
    grid.points(x = ix[1:m], y = ipar, default.units = "native", pch = pch[1], gp = gpar(fill = col[1]), name = paste(lab, "_ipar_points", sep = ""))

    ## plot seperation line and threshold parameters
    grid.lines(x = ix[m] + c(0.49, 0.51), y = ylim, default.units = "native", name = paste(lab, "_seperation_line", sep = ""))
    grid.lines(x = ix[(m + 1):(m + o)], y = cpar, default.units = "native", gp = gpar(lty = lty[2]), name = paste(lab, "_cpar_line", sep = ""))
    grid.points(x = ix[(m + 1):(m + o)], y = cpar, default.units = "native",
                pch = pch[2], gp = gpar(fill = col[2]), , name = paste(lab, "_cpar_points", sep = ""))

        
    ## add box and axis
    grid.rect(name = paste(lab, "_plot-box", sep = ""))
    grid.xaxis(at = ix, label = names, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""))
    grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
    upViewport()
    
    ## add descriptions
    pushViewport(lmargin.vp)
    grid.text(ylab, x = 0.2, rot = 90, name = paste(lab, "_ylab-left", sep = ""))
    upViewport(2)
    
    ## go back to uper vp
    upViewport(2)
  }
  
  panelfun
}
class(node_rsmplot) <- "grapcon_generator"
