prepanel.parallel.horizontal <- 
function (x, y, z, horizontal = TRUE, ...) 
{
  if (horizontal) list(xlim = extend.limits(c(1, ncol(as.data.frame(z))), prop = 0.03), ylim = c(0, 1), dx = 1, dy = 1)
  else list(xlim = c(0, 1), ylim = extend.limits(c(1, ncol(as.data.frame(z))), prop = 0.03), dx = 1, dy = 1)
}

panel.parallel.horizontal <- 
function (x, y, z, subscripts, groups = NULL, col = superpose.line$col, 
    lwd = superpose.line$lwd, lty = superpose.line$lty, alpha = superpose.line$alpha, 
    common.scale = FALSE, lower = sapply(z, function(x) min(as.numeric(x), 
        na.rm = TRUE)), upper = sapply(z, function(x) max(as.numeric(x), 
        na.rm = TRUE)), horizontal = TRUE, ...) 
{
    superpose.line <- lattice::trellis.par.get("superpose.line")
    reference.line <- lattice::trellis.par.get("reference.line")
    n.r <- ncol(z)
    n.c <- length(subscripts)
    if (is.null(groups)) {
        col <- rep(col, length = n.c)
        lty <- rep(lty, length = n.c)
        lwd <- rep(lwd, length = n.c)
        alpha <- rep(alpha, length = n.c)
    }
    else {
        groups <- as.factor(groups)[subscripts]
        n.g <- nlevels(groups)
        gnum <- as.numeric(groups)
        col <- rep(col, length = n.g)[gnum]
        lty <- rep(lty, length = n.g)[gnum]
        lwd <- rep(lwd, length = n.g)[gnum]
        alpha <- rep(alpha, length = n.g)[gnum]
    }
    if (is.function(lower)) 
        lower <- sapply(z, lower)
    if (is.function(upper)) 
        upper <- sapply(z, upper)
    if (common.scale) {
        lower <- min(lower)
        upper <- max(upper)
    }
    lower <- rep(lower, length = n.r)
    upper <- rep(upper, length = n.r)
    dif <- upper - lower
    if (n.r > 1) {
      if (horizontal) lattice::panel.segments(y0 = 0, y1 = 1, x0 = seq_len(n.r), x1 = seq_len(n.r), 
                                     col = reference.line$col, lwd = reference.line$lwd, 
                                     lty = reference.line$lty)
      else lattice::panel.segments(x0 = 0, x1 = 1, y0 = seq_len(n.r), y1 = seq_len(n.r), 
                          col = reference.line$col, lwd = reference.line$lwd, 
                          lty = reference.line$lty)
    }else return(invisible())
    for (i in seq_len(n.r - 1)) {
        x0 <- (as.numeric(z[subscripts, i]) - lower[i])/dif[i]
        x1 <- (as.numeric(z[subscripts, i + 1]) - lower[i + 1])/dif[i + 
            1]
        if (horizontal) lattice::panel.segments(y0 = x0, x0 = i, y1 = x1, x1 = i + 1, 
            col = col, lty = lty, lwd = lwd, alpha = alpha, ...)
        else lattice::panel.segments(x0 = x0, y0 = i, x1 = x1, y1 = i + 1, 
                            col = col, lty = lty, lwd = lwd, alpha = alpha, ...)
    }
    invisible()
}

confidence.panel.boot <- function(x, y, z, subscripts, lwd = 1, SD = NULL, ..., lower, upper, range = c(0, 1)) {
  nc <- ncol(z)
  if (missing(lower)) lower <- sapply(z, function(x) quantile(x, range[1]))
  if (missing(upper)) upper <- sapply(z, function(x) quantile(x, range[2]))
  dif <- upper - lower
  if (!is.null(SD)) {
    SD <- lapply(SD, function(x) (x - lower)/dif)
    for (l in seq_along(SD)) {
      grid::grid.polygon(y = grid::unit(c(SD[[l]][,1], rev(SD[[l]][,3])), "native"),
                   x = grid::unit(c(seq_len(nc),rev(seq_len(nc))), "native"),
                   gp = grid::gpar(fill = rgb(190/225, 190/225, 190/225, 0.5), col = "darkgrey"))
    }
  }
  
  panel.parallel.horizontal(x, y, z, subscripts, ..., lower = lower, upper = upper) 
  if (!is.null(SD)) {
    for (l in seq_along(SD)) {
      lattice::llines(y = SD[[l]][,2], x = seq_len(nc), col="white", lwd=lwd, lty = 1)
    }
  }
}

setMethod("plot", signature(x = "FLXboot", y = "missing"), function(x, y, ordering = NULL, range = c(0, 1),
                                             ci = FALSE, varnames = colnames(pars), strip_name = NULL, ...) {
  k <- x@object@k
  pars <- parameters(x)
  if (ci) {
    x_refit <- refit(x@object)
    sd <- sqrt(diag(x_refit@vcov))
    CI <- x_refit@coef + qnorm(0.975) * cbind(-sd, 0, sd)
    indices_prior <- grep("alpha$", names(x_refit@coef))
    if (length(indices_prior)) {
      z <- rmvnorm(10000, x_refit@coef[indices_prior,drop=FALSE], x_refit@vcov[indices_prior,indices_prior,drop=FALSE])
      Priors <- t(apply(cbind(1, exp(z))/rowSums(cbind(1, exp(z))), 2, quantile, c(0.025, 0.5, 0.975)))
      indices <- lapply(seq_len(k), function(i) grep(paste("_Comp.", i, sep = ""), names(x_refit@coef[-indices_prior])))
      SD <- lapply(seq_len(k), function(i) rbind(CI[indices[[i]], ], prior = Priors[i,]))
    } else {
      indices <- lapply(seq_len(k), function(i) grep(paste("_Comp.", i, sep = ""), names(x_refit@coef)))
      SD <- lapply(seq_len(k), function(i) CI[indices[[i]], ])
      mnrow <- max(sapply(SD, nrow))
      SD <- lapply(SD, function(x) if (nrow(x) <  mnrow) do.call("rbind", c(list(x), as.list(rep(0, mnrow - nrow(x))))) else x)
    }
    if (any("gaussian" %in% sapply(x@object@model, function(x) if (is(x, "FLXMRglm")) x@family else ""))) {
      i <- grep("sigma$", colnames(pars))
      pars[,i] <- log(pars[,i])
      colnames(pars)[i] <- "log(sigma)"
    }
  } else SD <- NULL
  range_name <- vector(mode = "character", length=2)
  range_name[1] <- if (range[1] == 0) "Min" else paste(round(range[1]*100), "%", sep = "")
  range_name[2] <- if (range[2] == 1) "Max" else paste(round(range[2]*100), "%", sep = "")
  Ordering <- if (is.null(ordering)) NULL else factor(as.vector(apply(matrix(pars[,ordering], nrow = k), 2, function(x) order(order(x)))))
  if(is.null(strip_name)) formula =  ~ pars else {
    opt.old <- options(useFancyQuotes = FALSE)
    on.exit(options(opt.old))
    formula <- as.formula(paste("~ pars | ", sQuote(strip_name)))
  }
  pars <- na.omit(pars)
  if (!is.null(attr(pars, "na.action"))) 
    Ordering <- Ordering[-attr(na.omit(pars), "na.action")]
  parallel.plot <- lattice::parallelplot(formula, groups = Ordering, default.scales = list(y = list(at = c(0, 1), labels = range_name),
                                                                       x = list(alternating = FALSE, axs = "i", tck = 0, at = seq_len(ncol(pars)))), range = range,
                                         panel = confidence.panel.boot, prepanel = prepanel.parallel.horizontal, SD = SD, ...)
  parallel.plot$x.scales$labels <- varnames
  parallel.plot
})
