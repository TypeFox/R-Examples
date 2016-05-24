## convenience function for rainbow in HCL space
hclrainbow <- function(n) hcl(h = seq(0, 360 * (n - 1)/n, length = n), c = 80, l = 60)

## region plot graphic function
regionplot <- function (object, parg = list(type = NULL, ref = NULL, alias = TRUE),
  names = TRUE, main = NULL, xlab = "", ylab = "Latent trait", ylim = NULL,
  off = 0.1, col = NULL, linecol = 2, srt = 45, adj = c(1.1, 1.1), axes = TRUE, ...)
{
  ## process parg list
  type <- parg$type
  ref <- parg$ref
  alias  <- if (is.null(parg$alias)) TRUE else parg$alias

  ## fetch requested threshold parameters
  tp <- threshpar(object, type = type, ref = ref, alias = alias, vcov = FALSE)

  ## process argument 'col'
  oj <- sapply(tp, function (j) sum(!is.na(j)))
  if (is.null(col)) {
      cols <- lapply(oj, function (j) gray.colors((j + 1)))
  } else if (is.character(col)) {
      stopifnot(length(col) == (max(oj) + 1))
      cols <- lapply(oj, function (j) col[1:(j + 1)])
  } else if (is.list(col)) {
      stopifnot(length(col) == length(oj))
      oj2 <- lapply(col, function (j) length(j))
      stopifnot(all(oj2 >= oj + 1))
      cols <- col
  } else if (is.function(col)) {
      cols <- lapply(oj, function (j) col(j + 1))
  } else stop("Argument 'col' is misspecified (see ?regionplot for possible values).")

  ## setup par, number of items, axis labels
  ## opar <- par(mar = c(4.25, 4.25, if (is.null(main)) 1 else 3, 2.5))
  m <- length(tp)
  if (isTRUE(names)) {
    nms <- names(tp)
    if(is.null(nms)) nms <- paste0("Item", formatC(1:m, width = nchar(m), digits = 0, flag = "0"))
  }
  if (is.character(names)) {
    stopifnot(length(names) == m)
    nms <- names
    names <- TRUE
  }
  if(!names) {
    lab <- rep(NA, m)
    lab[c(1, m)] <- c(1, m)
    pr <- pretty(1:m)
    pr <- pr[pr > 1 & pr < m]
    lab[pr] <- pr    
    nms <- lab
  }
  
  ## check if all threshold parameters are in order, if not, calculate sorted ones
  us <- sapply(tp, is.unsorted)
  if (any(us)) {
    usj <- which(us)
    for (j in usj) {
      tpj <- tp[[j]]
      nj <- length(tpj)
      
      ## check if there is a point with a biggest parameter, if yes, take mean
      for (i in 1:nj) {
        if (all(tpj[i] > tpj[(i+1):nj])) {
          tpj[i] <- mean(tpj[i:nj])
          tpj <- tpj[-(i+1:nj)]
          break
        }
      }
      
      ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
      while(is.unsorted(tpj)) {
        uo_pos <- which(diff(tpj) < 0)                     # locate unordered parameters, returns position of the first
        tpj[uo_pos] <- (tpj[uo_pos] + tpj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
        tpj <- tpj[-(uo_pos + 1)]                          # remove second
      }
      
      tp[[j]] <- tpj
    }
  }

  ## setup axis range and positions
  if (is.null(ylim)) ylim <- extendrange(unlist(tp, use.names = FALSE), f = 0.25)
  xi <- 0:m + c(0:(m-1), m-1) * off
  xlim <- c(xi[1], xi[m+1])

  ## setup graphical window
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", xaxs = "i", yaxs = "i", axes = FALSE, ylab = ylab, xlab = xlab, main = main, ...)

  ## plot items 
  for (j in seq_along(tp)) {
    rect(xleft = xi[j], xright = xi[j] + 1, ybottom = c(ylim[1], tp[[j]]), ytop = c(tp[[j]], ylim[2]), col = cols[[j]])
  }
  
  ## indicate unordered parameters
  if ((is.null(type) || type == "mode")) {
    orgtp <- threshpar(object, type = type, ref = ref, alias = alias, vcov = FALSE)
    uo_items <- which(!sapply(mapply(all.equal, tp, orgtp, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
    for (j in uo_items) {
      uo_pars <- setdiff(orgtp[[j]], tp[[j]])
      lines(x = rep(c(xi[j], xi[j] + 1, NA), length(uo_pars)), y = rep(uo_pars, each = 3), col = linecol, lty = 2)
    }
  }

  ## add axes
  if(axes) {
    axis(2)
    axis(4)
    if(names) {
      text(xi[-(m+1)] + 0.5, par("usr")[3], labels = nms, srt = srt, adj = adj, xpd = TRUE, cex = 0.9)
    } else {
      axis(1, at=(xi[-(m+1)] + 0.5), labels = nms)
    }
  }
  box()
  ## on.exit(par(opar))
}

## profiles of item, threshold or discrimination parameters
profileplot <- function(object, what = c("items", "thresholds", "discriminations"),
  parg = list(type = NULL, ref = NULL, alias = TRUE), index = TRUE, names = TRUE, main = NULL,
  abbreviate = FALSE, ref = TRUE, col = "lightgray", border = "black", pch = NULL, cex = 1, refcol = "lightgray",
  linecol = "black", lty = 2, ylim = NULL, xlab = NULL, ylab = NULL, add = FALSE,
  srt = 45, adj = c(1.1, 1.1), axes = TRUE, ...)
{
  ## check input
  what <- match.arg(what)
  if (what == "thresholds") type <- parg$type
  refpar <- parg$ref
  alias <- if (is.null(parg$alias)) TRUE else parg$alias
  addargs <- list(...)
  if ("difficulty" %in% names(addargs)) {
    warning("The argument 'difficulty' is deprecated and not longer used. All plotted parameters are difficulty parameters.")
  }
  if ("center" %in% names(addargs)) {
    warning("The argument 'center' is deprecated and not longer used. Centered parameters can be plotted when setting the 'ref' argument of 'parg' to NULL (default).")
  }

  ## parameters to be plotted
  if (what == "items") {
    cf <- itempar(object, ref = refpar, alias = alias, vcov = FALSE)
    ncf <- length(cf)
  } else if (what == "thresholds") {
    cf <- coef(threshpar(object, ref = refpar, alias = alias, type = type, vcov = FALSE), type = "matrix")
    ncf <- nrow(cf)
    ncat <- ncol(cf)
  } else {
    cf <- discrpar(object, ref = refpar, alias = alias, vcov = FALSE)
    ncf <- length(cf)
  }
  cf_ref <- mean(cf)

  ## labeling
  if (is.null(names)) names <- !index
  if (isTRUE(names)) nms <- if (what == "thresholds") rownames(cf) else names(cf)
  if (is.character(names)) {
    nms <- names
    names <- TRUE
  }
  if(!names & index) {
    lab <- rep(NA, ncf)
    lab[c(1, ncf)] <- c(1, ncf)
    pr <- pretty(1:ncf)
    pr <- pr[pr > 1 & pr < ncf]
    lab[pr] <- pr    
    nms <- lab
  }

  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(nms))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  nms <- abbreviate(nms, abbreviate)

  ## graphical parameter processing  
  type <- if(index) "b" else "p"
  missing_col_border <- !missing(col) && missing(border)
  if(what != "thresholds") {
    if(is.null(pch)) pch <- 21
    pch <- rep(pch, length.out = ncf)
    col <- rep(col, length.out = ncf)
    border <- rep(border, length.out = ncf)
  } else {
    col <- matrix(rep(t(col), length.out = ncf * ncat), ncol = ncat, byrow = TRUE)
    border <- matrix(rep(t(border), length.out = ncf * ncat), ncol = ncat, byrow = TRUE)
    linecol <- rep(linecol, length.out = ncat)
    if(!is.null(pch)) pch <- matrix(rep(t(pch), length.out = ncf * ncat), ncol = ncat, byrow = TRUE)
  }
  if((!index | is.null(pch)) && missing_col_border) border <- col
  cex <- rep(cex, length.out = ncf)
  if(is.null(ylim)) ylim <- range(cf, na.rm = TRUE)
  ylim <- rep(ylim, length.out = 2)
  xlim <- if (!index && what == "thresholds" && ncat > 1) extendrange(c(1, ncat), f = 0.25) else NULL
  lb <- if (what == "items") "difficulty" else if (what == "thresholds") "threshold" else "discrimination"
  if(is.null(xlab)) xlab <- if (index) "" else {
    if(what != "thresholds") "Items" else "Categories"
  }
  if(is.null(ylab)) ylab <- paste("Item", lb, "parameters")

  ## raw plot
  if(index) ix <- 1:ncf else if (what == "thresholds") {
    ix <- matrix(rep(1:ncat, each = ncf), nrow = ncf)
  } else ix <- rep(0, ncf)
  if (!add) {
      plot(ix, rep(0, length(ix)), xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ylim = ylim, xlim = xlim, main = main, ...)
      if(ref & what != "thresholds") abline(h = cf_ref, col = refcol)
      if(axes) axis(2)
      if (!index & what == "thresholds" & axes) axis(1, at = unique(ix), labels = paste("Category", 1:ncat))
      box()
  }

  ## actual data
  if(!index & names) {
    if (what == "thresholds") {
      for (i in 1:ncat) text(nms, x = ix[, i], y = cf[, i], col = border[, i])
    } else {
      text(nms, x = ix, y = cf, col = border)
    }
  } else {
    if (what == "thresholds") {
      for (i in 1:ncat) {
        lines(ix, cf[, i], type = type, lty = lty, pch = NA, col = linecol[i], cex = cex)
        if (is.null(pch)) {
          text(x = ix, y = cf[, i], labels = paste0("C", i), col = border[, i], cex = cex * 0.8, font = 2)
        } else {
	  lines(x = ix, y = cf[, i], type = type, lty = 0, pch = pch[, i], col = border[, i], bg = col[, i], cex = cex)
	}
      }
    } else {
      lines(ix, cf, type = type, lty = lty, pch = NA, col = linecol, cex = cex)
      lines(ix, cf, type = type, lty = 0, pch = pch, col = border, bg = col, cex = cex)
    }
    if(index && !add && axes) {
      if(names) {
        text(ix, par("usr")[3], labels = nms, srt = srt, adj = adj, xpd = TRUE, cex = 0.9)
      } else {
        axis(1, at = ix, labels = nms)
      }
    }
  }
}

## ICC/CCC/response curve functions plot
curveplot <- function (object, ref = NULL, items = NULL, names = NULL, layout = NULL,
                       xlim = NULL, ylim = c(0, 1), col = NULL, lty = NULL, main = NULL,
                       xlab = "Latent trait", ylab = "Probability", add = FALSE, ...)
{
  ## setup relevant informations and process input
  tp <- coef(threshpar(object, type = "mode", vcov = FALSE), type = "matrix")
  oj <- apply(tp, 1, function (j) sum(!is.na(j)) + 1)
  m0 <- nrow(tp)
  idx <- rep(1:m0, oj)
  nms <- if (!is.null(colnames(object$data))) colnames(object$data) else  paste0("Item", formatC(1:m0, width = nchar(m0), digits = 0, flag = "0"))
  if (is.null(items)) {
    items <- 1:m0
  } else if (is.numeric(items)) {
    stopifnot(all(items %in% 1:m0))
  } else if (is.character(items)) {
    stopifnot(all(items %in% nms))
    items <- which(nms %in% items)
  }
  m <- length(items)

  ## setup plotting parameters
  if (is.null(names)) {
    pnms <- nms
  } else {
    stopifnot(length(names) == m)
    pnms <- rep(NULL, m0)
    pnms[items] <- names
  }
  if (is.null(layout)) {
    nrow <- ceiling(sqrt(m))
    ncol <- ceiling(m/nrow)
    layout <- matrix(1:(nrow*ncol), nrow = nrow, ncol = ncol)
  } else stopifnot(prod(dim(layout)) >= m)
  if (is.null(xlim)) xlim <- extendrange(unclass(tp), f = 0.25)
  if (is.null(col)) col <- hclrainbow(max(oj))
  if (is.null(lty)) lty <- 1 else stopifnot(length(lty) %in% c(1, max(oj)))

  ## get probabilities
  theta <- seq(from = xlim[1], to = xlim[2], by = 0.1)
  pr <- predict(object, newdata = theta, type = "probability", ref = ref)

  ## setup plotting area, setup par
  if (!add) {
    layout(layout)
    opar <- par(mar = c(4.25, 4.25, 3, 2.5))
  } else {
    if (m > 1) stop("Overlaying response curves is only possible with a single item plotted.")
  }

  ## loop through items and plot CCC
  for (j in items) matplot(x = theta, y = pr[, j == idx], type = "l", main = if (is.null(main)) pnms[j] else main, xaxs = "i",
                           lty = lty, col = col, xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim, add = add, ...)

  if (!add) {
      layout(matrix(1, nrow = 1, ncol = 1))
      on.exit(par(opar))
  }
}

## person item plot
piplot <- function (object, ref = NULL, items = NULL, xlim = NULL, names = NULL,
                    labels = TRUE, main = "Person-Item Plot", xlab = "Latent trait",
                    abbreviate = FALSE, cex.axis = 0.8, cex.text = 0.5, cex.points = 1.5, ...)
{
  ## setup parameters, number of items and item labels
  tp <- threshpar(object, ref = ref, vcov = FALSE)
  ip <- sapply(tp, mean)
  pp <- personpar(object, ref = ref, vcov = FALSE)
  ppr <- range(as.numeric(names(pp)))
  rs <- rowSums(object$data, na.rm = TRUE)
  rs <- rs[rs >= ppr[1] && rs <= ppr[2]]
  ppt <- table(pp[rs])
  pptx <- as.numeric(names(ppt))
  m <- length(ip)
  nms <- names(ip)

  ## process argument items
  if (is.null(items)) {
    items <- 1:m
  } else if (is.numeric(items)) {
    stopifnot(all(items %in% 1:m))
  } else if (is.character(items)) {
    stopifnot(all(items %in% nms))
    items <- which(items %in% nms)
  } else stop("Argument 'items' is misspecified (see ?piplot for possible values).")
  
  ## subset to requested items
  m <- length(items)
  ip <- ip[items]
  tp <- tp[items]
  nms <- nms[items]

  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(nms))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  nms <- abbreviate(nms, abbreviate)
  w <- max(nchar(nms)) * 0.5

  ## setup x axis limits, backup par
  if (is.null(xlim)) xlim <- range(c(unlist(tp), pp))
  ylim <- c(1, m)

  ## setup graphic region
  layout(matrix(1:2, ncol = 1, nrow = 2), heights = c(1, 2))
  
  ## person parameter plot
  par(mar = c(0, w, 2.5, 1))
  plot(x = 0, xlim = xlim, ylim = c(0, max(ppt)), type = "n", axes = FALSE, xlab = "", ylab = "", main = main, cex.axis = cex.axis, ...)
  points(x = pptx, y = ppt, type = "h", col = "gray", lend = 2, lwd = 5, ...)
  box()
  
  ## item/threshold parameter plot
  opar <- par(mar = c(4.5, w, 0, 1))
  plot(x = 0, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = xlab, ylab = "", cex.axis = cex.axis)
  for (i in 1:m) {
    lines(y = rep(i, length(tp[[i]])), x = tp[[i]], type = "b", pch = 1, cex = cex.points, ...)
    if (labels) text(x = tp[[i]], y = rep(i, length(tp[[i]])), labels = gsub("C", "", names(tp[[i]])), cex = cex.text, ...)
  }
  points(x = ip, y = 1:m, pch = 16, cex = cex.points, ...)
  axis(side = 1, cex.axis = cex.axis)
  axis(side = 2, at = 1:m, labels = nms, las = 2, cex.axis = cex.axis)
  box()

  ## restore par
  on.exit(par(opar))
}

## information curve plot
infoplot <- function (object, what = c("categories", "items", "test"), ref = NULL, items = NULL, names = NULL,
                      layout = NULL, xlim = NULL, ylim = NULL, col = NULL, lty = NULL, lwd = NULL, main = NULL,
                      legend = TRUE, xlab = "Latent trait", ylab = "Information", add = FALSE, ...)
{
  ## process input
  what <- match.arg(what)

  ## setup items
  nms <- if (!is.null(colnames(object$data))) colnames(object$data) else paste0("Item", formatC(1:ncol(object$data), width = nchar(ncol(object$data)), digits = 0, flag = "0"))
  m <- length(nms)
  if (is.null(items)) {
    items <- 1:m
  } else if (is.numeric(items)) {
    stopifnot(all(items %in% 1:m))
  } else if (is.character(items)) {
    stopifnot(all(items %in% nms))
    items <- which(nms %in% items)
  }
  m <- length(items)

  ## setup plotting names
  if (is.null(names)) {
    pnms <- nms
  } else {
    stopifnot(length(names) == m)
    pnms <- rep(NULL, m)
    pnms[items] <- names
  }

  ## setup layout
  if (what == "test") {
    if (!is.null(layout)) warning("Argument 'layout' is not considered when visualizing test information.")
    layout <- matrix(1, nrow = 1, ncol = 1)
    overlay <- TRUE
  } else if (is.null(layout)) {
    if (what == "items") {
      layout <- matrix(1, nrow = 1, ncol = 1)
      overlay <- TRUE
    } else {
      nrow <- ceiling(sqrt(m))
      ncol <- ceiling(m/nrow)
      layout <- matrix(1:(nrow * ncol), nrow = nrow, ncol = ncol)
      overlay <- FALSE
    }
  } else {
    stopifnot(prod(dim(layout)) >= m)
    overlay <- FALSE
  }

  ## setup graphical stuff and information
  tp <- coef(threshpar(object, type = "mode", vcov = FALSE), type = "matrix")
  if (is.null(xlim)) xlim <- extendrange(unclass(tp), f = 0.5)
  theta <- seq(from = xlim[1], to = xlim[2], by = 0.05)
  if (is.null(lty)) lty <- 1
  if (what == "test") {
    type <- "test-information"
    if (is.null(col)) col <- "black"
    if (is.null(main)) main <- "Test Information"
  } else if (what == "items") {
    type <- "item-information"
    if (is.null(col)) col <- if (overlay) hclrainbow(m) else "black"
    if (is.null(main) & overlay) main <- "Item Information"
  } else {
    type <- "category-information"
    if (is.null(col)) col <- hclrainbow(ncol(tp) + 1)
  }
  info <- predict(object, ref = ref, type = type, newdata = theta)

  ## plot requested information
  if (!add) {
    layout(layout)
    opar <- par(mar = c(4.25, 4.25, 3, 2.5))
  } else {
    if (!overlay && m > 1 && what != "test") stop("Overlaying information curves is only possible for tests or single items.")
  }
  
  if (overlay) {
    if (what == "items") info <- info[, items]
    matplot(x = theta, y = info, xlim = xlim, ylim = ylim, type = "l", main = main, lty = lty, col = col, lwd = lwd,
            xlab = xlab, ylab = ylab, xaxs = "i", add = add, ...)
    if (what == "items" && !add && legend) legend(x = "topleft", legend = pnms[items], col = col, lwd = lwd, lty = lty, bty = "n")
  } else {
    for (j in items) {
      matplot(x = theta, y = info[, grepl(nms[j],colnames(info))], xlim = xlim, ylim = ylim, xaxs = "i",
              main = pnms[j], type = "l", lty = lty, lwd = lwd, col = col, xlab = xlab, ylab = ylab, add = add, ...)
    }
  }

  if (!add) {
    layout(matrix(1, nrow = 1, ncol = 1))
    on.exit(par(opar))
  }
}
