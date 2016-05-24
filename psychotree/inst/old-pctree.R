### high-level convenience interface for creating pctrees
pctree <- function (formula, data, minsplit = 30, nullcats = c("keep", "downcode", "ignore"),
                    reltol = 1e-10, deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L, ...)
{
  ## transform formula
  stopifnot(length(formula) > 2)
  formula <-  formula(terms(formula, data = data))
  ff <- y ~ 1 | x
  ff[[2]] <- formula[[2]]
  ff[[3]][[3]] <- formula[[3]]

  ## formula/data/model pre-processing
  pcmmod <- PCModel(nullcats = nullcats, reltol = reltol, deriv = deriv, hessian = hessian, maxit = maxit)
  ff <- attr(ParseFormula(ff), "formula")
  ff$input[[3]] <- ff$input[[2]]
  ff$input[[2]] <- ff$response[[2]]
  ff <- dpp(pcmmod, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), 
            data = data, na.action = na.pass)

  ## data sanity checks
  y <- as.matrix(ff@get("response"))
  if(ncol(y) < 3) stop("need at least three items")
  rmax <- sum(apply(y, 2, max, na.rm = TRUE))
  if(any(apply(y, 1, function(x) sum(x, na.rm = TRUE)) %in% c(0, rmax)))
    stop("Please remove subjects who only scored zero or the highest category.")

  ## call mob()
  res <- mob(ff, model = pcmmod, control = mob_control(minsplit = minsplit,
                                 objfun = function (object) -as.vector(logLik(object)), ...))

  ## add class and return
  structure(list(mob = res), class = "pctree")
}

logLik.pctree <- function (object, ...) logLik(object$mob, ...)

sctest.pctree <- function (x, ...) sctest(x$mob, ...)

weights.pctree <- function (object, ...) weights(object$mob, ...)

summary.pctree <- function (object, node = NULL, ...) summary(object$mob, ...)

print.pctree <- function (x, ...) {
  print(x$mob, ...)
  invisible(x)
}

plot.pctree <- function (x, terminal_panel = node_effects, tnex = 2, ...)
{
  plot(x$mob, terminal_panel = terminal_panel, tnex = tnex, tp_args = list(...))
}

coef.pctree <- function (object, node = NULL, ...) 
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

itempar.pctree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) itempar.PCModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

threshold.pctree <- function (object, node = NULL, ...) {

  object <- object$mob
  if(is.null(node)) node <- terminal_nodeIDs(object@tree)
  rval <- lapply(nodes(object, node), function(z) threshold.PCModel(z$model, ...))
  names(rval) <- node

  return(rval)

}

## terminal panel function for category characteristic curves
node_ccc <- function(mobobj, names = NULL, ref = NULL, ylab = "Latent trait",
                     ylim = NULL, col = NULL, n = 101, ccc_lwd = 1.25)
{
  ## check input
  stopifnot(!is.null(mobobj))

  ## get terminal nodes and thresholds
  tnodes <- terminal_nodeIDs(mobobj@tree)
  nodes_lst <- nodes(mobobj, tnodes)
  fun <- if (inherits(nodes_lst[[1]]$model, "RaschModel")) threshold.RaschModel else if (inherits(nodes_lst[[1]]$model, "RSModel")) threshold.RSModel else threshold.PCModel
  thresh_lst <- lapply(nodes_lst, function (node) fun(node$model, ref = ref, type = "unmodified", simplify = FALSE))
  item_lst <- lapply(thresh_lst, length)
  m_max <- max(unlist(item_lst))
  o_lst <- lapply(thresh_lst, function (thresh) sapply(thresh, length))
  o_max <- max(unlist(o_lst))
    
  ## setup plotting parameters
  if (is.null(names)) names <- paste0("I", 1:m_max) else stopifnot(length(names) == m_max)
  if (is.null(ylim)) ylim <- extendrange(unlist(thresh_lst, use.names = FALSE), f = 1)
  if (is.null(col)) col <- rainbow(o_max + 1) else stopifnot(length(col) == o_max + 1)

  ## setup theta and get probabilities for each node via psychotools:::ppcm()
  theta <- seq(from = ylim[1], to = ylim[2], length.out = n)
  probs_lst <- ppcm(theta = theta, delta = thresh_lst)

  ## setup label for extraction
  names(probs_lst) <- names(o_lst) <- names(item_lst) <- tnodes

  ## finally, panel function.
  panelfun <- function(node) {
    
    ## select node id, probabilities, identified items and number of categories of current node
    nid <- as.character(node$nodeID)
    probs <- probs_lst[[nid]]
    m <- item_lst[[nid]]
    o <- o_lst[[nid]] + 1

    ## setup xi, xlim, label
    xi <- 1:(m-1)
    xlim <- c(0, m)
    lab <- paste("node", nid, sep = "")

    ## terminal panel viewport setup
    top.vp <- viewport(layout = grid.layout(nrow = 2, ncol = 1, widths = unit(1, "null"), heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"), name = paste(lab, "_ccc", sep = ""))
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

    top.vp <- viewport(layout=grid.layout(nrow = 3, ncol = 3, widths = unit(wcol, c("lines", "null", "lines")),
                         heights = unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top-vp", sep = ""))
    bmargin.vp <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
    lmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
    rmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
    plot.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
    pushViewport(top.vp)
    pushViewport(plot.vp)

    ## plot ccc curves for each item
    for (i in seq_len(m)) {
      for (j in seq_len(o[i])) {
        grid.lines(x = - probs[[i]][, j] + i, y = theta, default.units = "native",
                   gp = gpar(col = col[j], lwd = ccc_lwd), name = paste(lab, "_item", i, "_cat", j, "_ccc", sep = ""))
      }
    }
    
    ## plot vertical seperation lines
    grid.polyline(x = rep(xi, each = 2), y = rep(ylim, m - 1), id.length = rep(2, m - 1),
                  default.units = "native", name = paste(lab, "_seperation_line", sep = ""))

    ## add box, labels and y-axis
    grid.rect(name = paste(lab, "_plot-box", sep = ""))
    grid.xaxis(at = c(0, xi) + 0.5, label = names, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""), gp = gpar(lty = 0))
    grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
    upViewport()
    
    ## add descriptions
    pushViewport(lmargin.vp)
    grid.text(ylab, x = 0.2, rot = 90, name = paste(lab, "_ylab-left", sep = ""))
    upViewport(2)
    
    ## go back to uper vp
    upViewport(2)
  }

  ## return
  return(panelfun)
}
class(node_ccc) <- "grapcon_generator"


## terminal panel function for effect displays
node_effects <- function(mobobj, names = NULL, type = c("mode", "median", "mean"),
                         ref = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, col_fun = gray.colors,
                         uo_show = TRUE, uo_col = "red", uo_lty = 2, uo_lwd = 1.25)
{
  ## check input
  stopifnot(!is.null(mobobj))
  stopifnot(off >= 0)
  type <- match.arg(type)

  ## get number of terminal nodes, corresponding nodes and thresholds
  tnodes <- terminal_nodeIDs(mobobj@tree)
  nodes_lst <- nodes(mobobj, tnodes)
  fun <- if (inherits(nodes_lst[[1]]$model, "RaschModel")) threshold.RaschModel else if (inherits(nodes_lst[[1]]$model, "RSModel")) threshold.RSModel else threshold.PCModel
  delta_lst <- lapply(nodes_lst, function (node) fun(node$model, type = type, ref = ref, simplify = FALSE))

  ## if requested and type = 'mode' check for unordered thresholds
  if (uo_show && type == "mode") {
    if (inherits(nodes_lst[[1]]$model, "RSModel")) {
      ip_lst <- lapply(nodes_lst, function (node) itempar.RSModel(node$model, ref = ref, vcov = FALSE, simplify = FALSE))
      ip_lst <- lapply(ip_lst, function (ip) lapply(as.list(ip[[1]]), function (beta) diff(0:length(ip[[2]]) * beta + c(0, ip[[2]]))))
    } else if (inherits(nodes_lst[[1]]$model, "PCModel")) {
      ip_lst <- lapply(nodes_lst, function (node) lapply(itempar.PCModel(node$model, ref = ref, vcov = FALSE, simplify = FALSE), function (j) diff.default(c(0, j))))
    } else {
      ip_lst <- lapply(nodes_lst, function (node) lapply(itempar.RaschModel(node$model, ref = ref, vcov = FALSE, simplify = FALSE), function (j) diff.default(c(0, j))))
    }

    names(ip_lst) <- tnodes
  }

  ## setup plotting parameters
  m <- max(sapply(delta_lst, length))
  xi <- 0:m + c(0:(m - 1), m - 1) * off
  xlim <- c(xi[1], xi[m + 1])

  ## setup axis range and labels
  if (is.null(ylim)) ylim <- extendrange(unlist(delta_lst, use.names = FALSE), f = 0.25)
  if (is.null(names)) names  <- paste("I", 1:m, sep = "") else stopifnot(length(names) == m)

  ## label for extraction
  names(delta_lst) <- tnodes

  ## finally, panel function, just select and paint.
  panelfun <- function(node) {
    
    ## select node id, coefficients, identified items, x-position vector and label
    nid <- as.character(node$nodeID)
    cf <- delta_lst[[nid]]
    lab <- paste("node", nid, sep = "")
 
    ## terminal panel viewport setup
    top.vp <- viewport(layout = grid.layout(nrow = 2, ncol = 1, widths = unit(1, "null"), heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"), name = paste(lab, "_effects", sep = ""))
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

    top.vp <- viewport(layout=grid.layout(nrow = 3, ncol = 3, widths = unit(wcol, c("lines", "null", "lines")), heights = unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top_vp", sep = ""))
    bmargin.vp <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
    lmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
    rmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
    plot.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
    pushViewport(top.vp)
    pushViewport(plot.vp)

    ## plot rectangles per item
    for (j in seq_along(cf)) {
      
      ncat <- length(cf[[j]]) + 1
      grid.rect(x = rep.int(xi[j], ncat), y = c(ylim[1], cf[[j]]), width = rep.int(1, ncat), height = diff.default(c(ylim[1], cf[[j]], ylim[2])),
                just = c("left", "bottom"), gp = gpar(fill = col_fun(ncat)), default.units = "native",
                name = paste(lab, "_item", j, "_rect", sep = ""))

    }

    ## if requested: indicate unordered parameters
    if (uo_show && type == "mode") {
      ip <- ip_lst[[nid]]
      uo_items <- which(!sapply(mapply(all.equal, cf, ip, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
      for (j in uo_items) {
        uo_pars <- setdiff(ip[[j]], cf[[j]])
        grid.polyline(x = rep(c(xi[j], xi[j] + 1), length(uo_pars)), y = rep(uo_pars, each = 2), 
                      default.units = "native", id = rep(1:length(uo_pars), each = 2),
                      name = paste(lab, "item", j, "_uolines", sep = ""),
                      gp = gpar(col = uo_col, lwd = uo_lwd, lty = uo_lty))
      }
    }

    ## add box and axis
    grid.rect(name = paste(lab, "_plot-box", sep = ""))
    grid.xaxis(at = (xi[-(m+1)] + 0.5), label = names, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""))
    grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
    upViewport()
    
    ## add descriptions
    pushViewport(lmargin.vp)
    grid.text(ylab, x = 0.2, rot = 90, name = paste(lab, "_ylab-left", sep = ""))
    upViewport(2)
    
    ## go back to uper vp
    upViewport(2)
  }

  ## return
  return(panelfun)
}
class(node_effects) <- "grapcon_generator"

### plotCCC.pctree() plots the category characteristic curves of a number of specified items for a given number of nodes
plotCCC.pctree <- function (object, nodes = NULL, items = NULL, byrow = TRUE, sum0 = TRUE, names = NULL, ylim = c(0, 1), xlim = NULL,
                            col = NULL, n = 101, main = "Category Characteristic Curves",  xlab = "Latent Trait", ylab = "Probability", ...)
{
  ## get available nodes if necessary, elsewise check input within nodes(...)
  nodes <- if (is.null(nodes)) terminal_nodeIDs(object$mob@tree) else nodes
  n_nodes <- length(nodes)
  nodes_lst <- nodes(object$mob, nodes)
  
  ##  set items to plot: if arguments items = NULL, then items is set to the maximum number of items (within all terminal nodes)
  n_items <- numeric(n_nodes)
  for (i in seq_len(n_nodes)) n_items[i] <- length(nodes_lst[[i]][["model"]][["items"]])
  n_items <- max(n_items)
  mvec <- 1:n_items
  if (is.null(items)) items <- mvec  else {
    stopifnot(is.numeric(items) && all(items %in% mvec))
    n_items <- length(items)
  }

  ## check custom item lab
  if (is.null(names)) names <- paste("Item", items) else stopifnot(length(names) == n_items)
  
  ## setup threshold parameters and number of categories per node
  oj_lst <- delta_lst <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    model <- nodes_lst[[i]]$model
    if (inherits(model, "RSModel")) {
      delta_lst[[i]] <- coef(model, type = "PCM", sum0 = sum0, as.list = TRUE)[items]
      oj_lst[[i]] <- model$categories[items]  + 1
    } else {
      delta_lst[[i]] <- coef(model, type = "threshold", sum0 = sum0, as.list = TRUE)[items]
      oj_lst[[i]] <- (sapply(model$categories, length, USE.NAMES = FALSE) + 1)[items]
    }
  }

  ## setup plotting parameters
  if (is.null(xlim)) xlim <- extendrange(unlist(delta_lst, use.names = FALSE), f = 1.75)
  if (is.null(col)) col <- rainbow(max(unlist(oj_lst, use.names = FALSE)))

  ## setup theta and get probabilities for each node via psychotools:::ppcm()
  theta <- seq(from = xlim[1], to = xlim[2], length.out = n)
  probs_lst <- ppcm(theta = theta, delta = delta_lst)

  ## setup plotting area
  top.vp <- viewport(layout = grid.layout(nrow = 4, ncol = 3,
                       heights = unit(c(3, 1, 5, 1), c("lines", "null", "lines", "lines")),
                       widths = unit(c(4, 1, 1), c("lines", "null", "lines"))), name = "top.vp")
  lyt <- if (byrow) grid.layout(nrow = n_nodes, ncol = 1) else grid.layout(nrow = 1, ncol = n_nodes)
  node.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = "node.vp", layout = lyt)
  grid.newpage()
  pushViewport(top.vp)
  pushViewport(node.vp)

  ## plot CCCs for all nodes
  for (i in seq_len(n_nodes)) {
    lab <- paste("node.vp.node", i, sep = "")
    if (byrow) {
      vp <- viewport(layout.pos.row = i, name = paste(lab, ".vp", sep = ""),
                     layout = grid.layout(1, n_items), xscale = xlim)
    } else {
      vp <- viewport(layout.pos.col = i, name = paste(lab, ".vp", sep = ""),
                     layout = grid.layout(n_items, 1) , xscale = xlim)
    }
    pushViewport(vp)
    
    ## plot all ccc for current node i
    for (j in seq_len(n_items)) {
      ## plotting area
      nlab <- paste(lab, ".item", items[j], sep = "")
      if (byrow) vp <- viewport(layout.pos.col = j, name = nlab, layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 1), c("lines", "null"))))
      else vp <- viewport(layout.pos.row = j, name = nlab, layout = grid.layout(nrow = 3, ncol = 1, heights = unit(c(1, 1), c("lines", "null"))))
      pushViewport(vp)
      grid.rect(name = paste(nlab, ".box", sep = ""))

      ## strip
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, name = paste(nlab, ".strip.vp", sep = "")))
      grid.rect(name = paste(nlab, ".strip", ".box", sep = ""), gp = gpar(fill = gray(0.8), ...))
      grid.text(paste("Node", nodes[i], "/", names[j]), name = paste(nlab, ".strip", ".text", sep = ""), gp = gpar(...))
      upViewport()

      ## main ccc plot
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(nlab, ".ccc.vp", sep = ""), xscale = xlim, yscale = ylim))
      
      ## axis
      if (byrow) {
        if (i == n_nodes) grid.xaxis(main = TRUE, name = paste(nlab, ".ccc", ".xaxis.bottom", sep = ""), gp = gpar(...))
        if (j == 1) grid.yaxis(main = TRUE, name = paste(nlab, ".ccc", ".yaxis.left", sep = ""), gp = gpar(...))
      } else {
        if (i == 1) grid.yaxis(main = TRUE, name = paste(nlab, ".ccc", ".yaxis.left", sep = ""), gp = gpar(...))
      }

      ## finally, plot all CCC's for item j and node i.
      ncat <- oj_lst[[i]][j]
      grid.polyline(x = rep(theta, ncat), y = probs_lst[[i]][[j]], id.lengths = rep(n, ncat),
                    default.units = "native", gp = gpar(col = col, ...), name = paste(nlab, ".ccc", sep = ""))

      ## go up again
      upViewport(2)
    }
    ## box and x-axis
    grid.rect(name = paste(lab, ".box", sep = ""))
    if (!byrow) grid.xaxis(main = TRUE, name = paste(lab, ".xaxis.top", sep = ""))
    upViewport()
  }
  upViewport()

  ## add main
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, name = "top.vp.title.vp"))
  grid.text(main, name = "top.vp.title.text", gp = gpar(...))
  upViewport()

  ## add ylab
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, name = "top.vp.ylab.vp"))
  grid.text(ylab, x = 0.25, rot = 90, name = "top.vp.ylab.text", gp = gpar(...))
  upViewport()

  ## add xlab
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2, name = "top.vp.xlab.vp"))
  grid.text(xlab, y = 0.4, name = "top.vp.xlab.text", gp = gpar(...))
  upViewport()

  ## add legend
  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2, name = "top.vp.legend.vp"))
  max_oj <- max(unlist(oj_lst, use.names = FALSE))
  x <- seq(0, 1, length.out = 2 * max_oj + 1)
  y <- rep(1, 2 * max_oj)
  id <- rep(1:(max_oj), each = 2)
  grid.polyline(x = x[1:(2 * max_oj)], y = y, id = id, gp = gpar(col = col, ...), name = "top.vp.legend.lines")
  grid.text(paste("Cat.", 1:max_oj), x = (x[seq(from = 2, by = 2, length.out = max_oj)]+x[2] * 0.5), y = rep(1, max_oj), name = "top.vp.legend.text")
  upViewport()
}
