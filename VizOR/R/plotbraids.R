##' Braided streams plot for cohort visualization
##'
##' Displays 'paths' taken by individuals passing through a sequence of discrete states, such as a
##' sequence of treatments.
##' @title plotbraids
##' @param formula A formula of the form \code{trt ~ seq} or \code{trt ~ seq | cond}, where
##' \code{trt} is a treatment factor, \code{seq} is an integer sequence number, and the optional
##' \code{cond} is a conditioning factor used to trellis the braided stream plot. Probably only
##' two-valued conditioning factors will produce visually acceptable plots.
##' @param data A data frame with columns named on the LHS and RHS of argument \code{formula}
##' @param idvar A character vector naming columns of \code{data} that identify multiple records
##' from the same individual, used to reshape \code{data} into wide form
##' @param stratify A logical value, indicating whether whitespace should be introduced to stratify
##' the braids by initial treatment
##' @param steps Which values of the sequence number should be included in the plot. Presently, only
##' vectors of the form \code{1:n} (for some integer \code{n}) are supported.
##' @param outside A logical value determining whether the state labels are to appear outside the
##' panel.
##' @param xlab x-axis label
##' @param \dots Additional arguments passed to delegates
##' @param x.scales.labels Labels for the treatment states
##' @param x.scales Provides a hook for modifying the basically sensible default layout and labeling
##' of the x-axis
##' @param scales Provides a hook for modifying the basically sensible default layout of both x- and
##' y-axes
##' @param strip Allows user to provide a strip function if the default does not suffice
##' @param lattice.options Allows specification of plot-specific lattice options
##' @keywords hplot
##'
##' @examples
##' ## We demonstrate a simple braided stream plot based on the built-in occupationalStatus data set.
##' ## It should be noted that the semantics of these data are very slightly at odds with the intended
##' ## application of the braided stream plot, since the index cases in this data set were the _sons_,
##' ## rather than the fathers. Thus, although time goes left-to-right in this figure, the streamlines
##' ## run right-to-left epidemiologically. Notwithstanding this minor technicality, the figure gives
##' ## a lively, compelling and meaningful visualization of these data.
##' ## 1. Build a 'wide-form' data set from the table 'occupationalStatus'
##' df.wide <- data.frame(status.1=rep(1:8, 8),
##'                       status.2=rep(1:8, each=8),
##'                       N=as.vector(occupationalStatus))
##' df.wide <- df.wide[rep(1:64, times=df.wide$N),-3]
##' ## 2. Reshape this to the 'long-form' data set expected by 'plotbraids'
##' df.long <- reshape(df.wide, varying=paste("status", 1:2, sep="."), direction="long", timevar="gen")
##' df.long <- df.long[order(df.long$id),]
##' ## TODO: Generate appropriate 'class' labels for status.
##' ## TODO: Use this opportunity to demonstrate meaningful application of a colored factor.
##' ## 3. Plot the braided stream plot
##' plotbraids(status ~ gen, df.long, stratify=TRUE, steps=1:2,
##'            outside=TRUE, xlab="Generation",
##'            x.scales.labels=c("Father","Son"))
##' @return A \code{trellis} plot object
##' @author David C. Norris
##
##' @export plotbraids
plotbraids <- function(formula, data, idvar="id", stratify=FALSE, steps=1:3,
                       outside=FALSE, xlab=NULL, ...,
                       x.scales.labels=paste(formula[[2]],steps,sep="."),
                       x.scales=list(
                         alternating=FALSE,
                         relation="same",
                         labels=x.scales.labels,
                         at=seq(length(steps))-0.5),
                       scales=list(
                         x=x.scales,
                         y=list(
                           draw=FALSE,
                           relation="free")),
                       strip=TRUE,
                       lattice.options=list("axis.padding"=list(numeric=0.0, factor=0.0))){
  ## If the trt column of data is not a colored factor, then convert it to one
  statevar <- as.character(formula[[2]])
  if(!is(data[[statevar]], "colored")){
    warning("Coercing state variable '", statevar, "' to a colored factor, using arbitrary colors.")
    data[[statevar]] <- colored(data[[statevar]]) # invoked without explicit color.key, 'colored' issues Warning
  }
  dots <- list(...)
  if (!is.null(lattice.options)) {
    oopt <- lattice.options(lattice.options)
    on.exit(lattice.options(oopt), add = TRUE)
  }
  form <- latticeParseFormula(formula, data,
                              multiple = FALSE, outer = FALSE, 
                              subscripts = TRUE, drop = TRUE)
  panel <- function(x, y, subscripts, ...)
    panel.plotbraids(x, y, subscripts, ..., formula=formula, data=data, idvar=idvar, stratify=stratify, steps=steps, outside=outside)
  if (!is.function(strip)) 
    strip <- eval(strip)
  subscr <- form$subscr
  cond <- form$condition
  y <- form$left
  x <- form$right
  if (length(cond) == 0) {
    strip <- FALSE
    cond <- list(gl(1, length(x)))
  }
  if (missing(xlab)) 
    xlab <- form$right.name
  ## Copied verbatim from private package:lattice, then simplified
  ## TODO: Revisit the whole design of this function, to delegate more intelligently to lattice.
  ##       The current 'design' was likely driven by a combination of urgency with the special
  ##       demands of panelwise processing the braids data. Performing such processing steps
  ##       up-front -- even e.g. creating an object of class 'braids' -- ought to resolve this.
  trellis.skeleton <- function (formula = NULL, cond, aspect = default.args$aspect, 
                                as.table = default.args$as.table, between = default.args$between, 
                                key = NULL, legend = NULL, page = default.args$page,
                                main = default.args$main, sub = default.args$sub,
                                par.strip.text = default.args$par.strip.text, 
                                layout = default.args$layout, skip = default.args$skip,
                                strip = default.args$strip.default, 
                                strip.left = FALSE, xlab.default = NULL, ylab.default = NULL, 
                                xlab = NULL, ylab = NULL, xlab.top = NULL, ylab.right = NULL, 
                                panel, xscale.components = default.args$xscale.components, 
                                yscale.components = default.args$yscale.components, axis = default.args$axis, 
                                subscripts = TRUE, index.cond = NULL, perm.cond = NULL, ..., 
                                par.settings = NULL, plot.args = NULL, lattice.options = NULL){
      default.args <- lattice.getOption("default.args")
      if (is.null(skip)) 
          skip <- FALSE
      foo <- list(formula = formula, as.table = as.table, aspect.fill = (aspect == 
                                         "fill"), legend = NULL, #construct.legend(legend = legend, key = key), 
                  panel = panel, page = page, layout = layout, skip = skip, 
                  strip = if (is.logical(strip) && strip) "strip.default" else strip, 
                  strip.left = if (is.logical(strip.left) && strip.left) strip.custom(horizontal = FALSE) else strip.left, 
                  xscale.components = xscale.components, yscale.components = yscale.components, 
                  axis = axis, xlab = xlab, ylab = ylab, xlab.default = xlab.default, 
                  ylab.default = ylab.default, xlab.top = xlab.top, ylab.right = ylab.right, 
                  main = main, sub = sub, x.between = 0, y.between = 0, 
                  par.settings = par.settings, plot.args = plot.args, lattice.options = lattice.options, 
                  par.strip.text = par.strip.text, index.cond = index.cond, 
                  perm.cond = perm.cond)
      if (!is.null(between$x)) 
          foo$x.between <- between$x
      if (!is.null(between$y)) 
          foo$y.between <- between$y
      foo$condlevels <- lapply(cond, levels)
      list(foo = foo, dots = list(...))
  }
  aspect <- "fill"
  foo <- do.call("trellis.skeleton",
                 c(list(formula = formula, 
                        cond = cond, aspect = aspect, strip = strip, panel = panel, 
                        xlab = xlab, ylab = NULL, xlab.default = form$right.name, 
                        ylab.default = NULL, lattice.options = lattice.options), 
                   dots))
  dots <- foo$dots
  foo <- foo$foo
  foo$call <- sys.call(sys.parent())
  foo$call[[1]] <- quote(xyplot) # TODO: Change this to 'plotbraids'?
  ## TODO: Handle scale construction here in case outside=TRUE,
  ##       so as to allocate sufficient space for the labels?
  ## Copied verbatim from private lattice function, eliminating use of ':::' to pass CRAN check:
  construct.scales <- function (draw = TRUE, axs = "r", tck = 1, tick.number = 5, at = FALSE, 
                                labels = FALSE, log = FALSE, alternating = TRUE, relation = "same", 
                                abbreviate = FALSE, minlength = 4, limits = NULL, format = NULL, 
                                equispaced.log = TRUE, lty = FALSE, lwd = FALSE, cex = FALSE, 
                                rot = FALSE, col = FALSE, col.line = col, alpha = FALSE, 
                                alpha.line = alpha, font = FALSE, fontfamily = FALSE, fontface = FALSE, 
                                lineheight = FALSE, ..., x = NULL, y = NULL){
      x.scales <- y.scales <- list(draw = draw, axs = axs, tck = tck, 
                                   tick.number = tick.number, at = at, labels = labels, 
                                   log = log, alternating = alternating, relation = relation, 
                                   abbreviate = abbreviate, minlength = minlength, limits = limits, 
                                   format = format, equispaced.log = equispaced.log, lty = lty, 
                                   lwd = lwd, cex = cex, rot = rot, col = col, col.line = col.line, 
                                   alpha = alpha, alpha.line = alpha.line, font = font, 
                                   fontfamily = fontfamily, fontface = fontface, lineheight = lineheight)
      if (!is.null(x)) {
          if (is.character(x)) 
              x <- list(relation = x)
          #x <- complete_names(x, x.scales)
          x.scales[names(x)] <- x
      }
      if (!is.null(y)) {
          if (is.character(y)) 
              y <- list(relation = y)
          #y <- complete_names(y, y.scales)
          y.scales[names(y)] <- y
      }
      if (is.logical(x.scales$alternating)) 
          x.scales$alternating <- if (x.scales$alternating) 
              c(1, 2)
          else 1
      if (is.logical(y.scales$alternating)) 
          y.scales$alternating <- if (y.scales$alternating) 
              c(1, 2)
          else 1
      for (nm in c("tck", "cex", "rot")) {
          x.scales[[nm]] <- rep(x.scales[[nm]], length.out = 2)
          y.scales[[nm]] <- rep(y.scales[[nm]], length.out = 2)
      }
      if (x.scales$relation == "same" && (is.list(x.scales$at) || 
              is.list(x.scales$labels))) 
          stop("the 'at' and 'labels' components of 'scales' may not be lists when 'relation = \"same\"'")
      if (y.scales$relation == "same" && (is.list(y.scales$at) || 
              is.list(y.scales$labels))) 
          stop("the 'at' and 'labels' components of 'scales' may not be lists when 'relation = \"same\"'")
      list(x.scales = x.scales, y.scales = y.scales)
  }
  foo <- c(foo, do.call("construct.scales", scales))
  cond.max.level <- unlist(lapply(cond, nlevels))
  foo$panel.args.common <- dots
  npackets <- prod(cond.max.level)
  if (npackets != prod(sapply(foo$condlevels, length))) 
    stop("mismatch in number of packets")
  foo$panel.args <- vector(mode = "list", length = npackets)
  foo$packet.sizes <- numeric(npackets)
  if (npackets > 1) {
    dim(foo$packet.sizes) <- sapply(foo$condlevels, length)
    dimnames(foo$packet.sizes) <- lapply(foo$condlevels, 
                                         as.character)
  }
  cond.current.level <- rep(1, length(cond))
  ## Copied verbatim from lattice, to eliminate ':::' for CRAN check:
  compute.packet <- function (cond, levels) {
      id <- !(do.call("pmax", lapply(cond, is.na)))
      stopifnot(any(id))
      for (i in seq_along(cond)) {
          var <- cond[[i]]
          id <- id & (if (is.shingle(var)) 
                      ((var >= levels(var)[[levels[i]]][1]) & (var <= levels(var)[[levels[i]]][2]))
          else (as.numeric(var) == levels[i]))
      }
      id
  }
  ## Copied verbatim from lattice, to eliminate ':::' for CRAN check:
  cupdate <- function (index, maxim) {
      if (length(index) != length(maxim) || length(maxim) <= 0) 
          stop("Inappropriate arguments")
      index[1] <- index[1] + 1
      if (index[1] > maxim[1] && length(maxim) > 1) 
          c(1, cupdate(index[-1], maxim[-1]))
      else index
  }
  ## Copied verbatim from lattice:::getFunctionOrName, for use in 'limits.and.aspect' below:
  getFunctionOrName <- function (FUN) {
      if (is.function(FUN)) 
          FUN
      else if (is.character(FUN)) 
          get(FUN)
      else eval(FUN)
  }
  ## Copied verbatim from lattice:::extend.limits, for use in 'limits.and.aspect' below:
  extend.limits <- function (lim, length = 1, axs = "r",
                             prop = if (axs == "i") 0 else lattice.getOption("axis.padding")$numeric) {
      if (all(is.na(lim))) 
          NA_real_
      else if (is.character(lim)) {
          c(1, length(lim)) + c(-1, 1) * if (axs == "i") 
              0.5
          else lattice.getOption("axis.padding")$factor
      }
      else if (length(lim) == 2) {
          if (lim[1] > lim[2]) {
              ccall <- match.call()
              ccall$lim <- rev(lim)
              ans <- eval.parent(ccall)
              return(rev(ans))
          }
          if (!missing(length) && !missing(prop)) 
              stop("'length' and 'prop' cannot both be specified")
          if (length <= 0) 
              stop("'length' must be positive")
          if (!missing(length)) {
              prop <- (as.numeric(length) - as.numeric(diff(lim)))/(2*as.numeric(diff(lim)))
          }
          if (lim[1] == lim[2]) 
              lim + 0.5 * c(-length, length)
          else {
              d <- diff(as.numeric(lim))
              lim + prop * d * c(-1, 1)
          }
      }
      else {
          print(lim)
          stop("improper length of 'lim'")
      }
  }
  ## Copied verbatim from lattice:::limitsFromLimitList, for use in 'limits.and.aspect' below:
  limitsFromLimitlist <- function (have.lim, lim, relation, limitlist,
                                   used.at, numlimitlist, axs, npackets) {
      if (relation == "same") {
          all.na <- unlist(lapply(limitlist, function(x) all(is.na(x))))
          class.lim <- lapply(limitlist[!all.na], class)
          limits <- unlist(limitlist)
          if (sum(!is.na(limits)) > 0) {
              if (is.character(limits)) {
                  limits <- unique(limits[!is.na(limits)])
                  slicelen <- diff(extend.limits(limits, axs = axs))
              }
              else {
                  limits <- extend.limits(range(as.numeric(limits), 
                                                finite = TRUE), axs = axs)
                  slicelen <- diff(range(limits, finite = TRUE))
              }
              if (length(class.lim) > 0) 
                  class(limits) <- if (all(class.lim[[1]] == "integer")) 
                      "numeric"
                  else class.lim[[1]]
          }
          else {
              limits <- c(0, 1)
              slicelen <- 1
          }
          if (have.lim) {
              if (is.list(lim)) 
                  stop("limits cannot be a list when relation = same")
              old.limits <- limits
              limits <- lim
              if (!is.character(limits) && !is.character(old.limits)) {
                  limits[is.na(limits)] <- old.limits[is.na(limits)]
              }
              slicelen <- if (is.character(limits)) 
                  length(limits) + 2
              else diff(range(as.numeric(limits)))
          }
          ans <- list(limits = limits, slicelen = slicelen)
      }
      else if (relation == "sliced") {
          if (have.lim) {
              if (is.list(lim)) {
                  limits <- rep(lim, length.out = npackets)
              }
              else warning("Explicitly specified limits ignored")
          }
          slicelen <- limitlist
          for (i in seq_along(limitlist)) {
              slicelen[[i]] <- if (!is.character(limitlist[[i]])) {
                  if (any(is.finite(limitlist[[i]]))) 
                      diff(range(as.numeric(limitlist[[i]]), finite = TRUE))
                  else NA_real_
              }
              else if (!any(is.na(numlimitlist[[i]]))) 
                  diff(range(as.numeric(numlimitlist[[i]])))
              else NA_real_
          }
          slicelen <- (if (axs == "i") 
                       1
          else 1 + 2 * lattice.getOption("axis.padding")$numeric) * 
              max(unlist(slicelen), na.rm = TRUE)
          for (i in seq_along(limitlist)) {
              if (is.numeric(limitlist[[i]])) 
                  limitlist[[i]] <- extend.limits(limitlist[[i]], 
                                                  length = slicelen)
          }
          for (i in seq_along(numlimitlist)) {
              if (!all(is.na(numlimitlist[[i]]))) 
                  numlimitlist[[i]] <- extend.limits(as.numeric(numlimitlist[[i]]), 
                                                     length = slicelen)
          }
          ans <- list(limits = limitlist, used.at = used.at, numlimitlist = numlimitlist, 
                      slicelen = slicelen)
      }
      else if (relation == "free") {
          if (have.lim) {
              if (!is.list(lim)) 
                  lim <- list(lim)
              id <- which(sapply(limitlist, function(x) !all(is.na(x))))
              old.limitlist <- limitlist
              limitlist[id] <- lim
              which.null <- sapply(limitlist, is.null)
              limitlist[which.null] <- old.limitlist[which.null]
              for (i in seq_along(limitlist)) {
                  if (!is.character(limitlist[[i]]) && !is.character(old.limitlist[[i]])) {
                      isna <- is.na(limitlist[[i]])
                      limitlist[[i]][isna] <- old.limitlist[[i]][isna]
                  }
              }
          }
          for (i in seq_along(limitlist)) {
              if (!all(is.na(limitlist[[i]])) && !is.character(limitlist[[i]])) 
                  limitlist[[i]] <- extend.limits(limitlist[[i]], 
                                                  axs = axs)
          }
          slicelen <- numeric(length(limitlist))
          for (i in seq_along(limitlist)) slicelen[i] <- if (!is.character(limitlist[[i]])) 
              diff(range(as.numeric(limitlist[[i]])))
          else if (!any(is.na(numlimitlist[[i]]))) 
              diff(range(numlimitlist[[i]]))
          else NA_real_
          ans <- list(limits = limitlist, used.at = used.at, numlimitlist = numlimitlist, 
                      slicelen = slicelen)
      }
      ans
  }
  ## Copied verbatim from lattice, to eliminate ':::' for CRAN check:
  limits.and.aspect <- function (prepanel.default, prepanel = NULL, have.xlim = FALSE, 
                                 xlim = NULL, have.ylim = FALSE, ylim = NULL, x.relation, 
                                 y.relation, panel.args.common = list(), panel.args = list(), 
                                 aspect, banking = lattice.getOption("banking"),
                                 npackets = length(panel.args), x.axs = "r", y.axs = "r", ...) {
      prepanel.default.function <- getFunctionOrName(prepanel.default)
      prepanel <- getFunctionOrName(prepanel)
      if (npackets < 1) 
          stop("need at least one panel")
      x.limits <- vector("list", npackets)
      y.limits <- vector("list", npackets)
      x.used.at <- vector("list", npackets)
      y.used.at <- vector("list", npackets)
      x.num.limit <- vector("list", npackets)
      y.num.limit <- vector("list", npackets)
      dxdy <- vector("list", npackets)
      for (count in seq_len(npackets)) {
          if (is.list(panel.args[[count]])) {
              pargs <- c(panel.args.common, panel.args[[count]], 
                         list(...))
              tem <- do.call("prepanel.default.function", pargs)
              if (is.function(prepanel)) {
                  prenames <- names(formals(prepanel))
                  if (!("..." %in% prenames)) 
                      pargs <- pargs[intersect(names(pargs), prenames)]
                  pretem <- do.call("prepanel", pargs)
                  if (!is.null(pretem$xlim) && !is.character(pretem$xlim)) 
                      if (any(isna <- is.na(pretem$xlim))) 
                          pretem$xlim[isna] <- tem$xlim[isna]
                  if (!is.null(pretem$ylim) && !is.character(pretem$ylim)) 
                      if (any(isna <- is.na(pretem$ylim))) 
                          pretem$ylim[isna] <- tem$ylim[isna]
#                  tem <- updateList(tem, pretem)
                  tem <- modifyList(if (is.null(tem)) list() else tem,
                                    pretem)
              }
              x.limits[[count]] <- tem$xlim
              y.limits[[count]] <- tem$ylim
              x.used.at[[count]] <- if (is.null(tem$xat)) 
                  NA
              else tem$xat
              y.used.at[[count]] <- if (is.null(tem$yat)) 
                  NA
              else tem$yat
              x.num.limit[[count]] <- if (is.null(tem$xat)) 
                  NA
              else range(tem$xat)
              y.num.limit[[count]] <- if (is.null(tem$yat)) 
                  NA
              else range(tem$yat)
              dxdy[[count]] <- list(dx = tem$dx, dy = tem$dy)
          }
          else {
              x.limits[[count]] <- c(NA_real_, NA_real_)
              y.limits[[count]] <- c(NA_real_, NA_real_)
              x.used.at[[count]] <- NA_real_
              y.used.at[[count]] <- NA_real_
              x.num.limit[[count]] <- NA_real_
              y.num.limit[[count]] <- NA_real_
              dxdy[[count]] <- list(dx = NA_real_, dy = NA_real_)
          }
      }
      x.limits <- limitsFromLimitlist(have.lim = have.xlim, lim = xlim, 
                                      relation = x.relation, limitlist = x.limits, used.at = x.used.at, 
                                      numlimitlist = x.num.limit, axs = x.axs, npackets = npackets)
      y.limits <- limitsFromLimitlist(have.lim = have.ylim, lim = ylim, 
                                      relation = y.relation, limitlist = y.limits, used.at = y.used.at, 
                                      numlimitlist = y.num.limit, axs = y.axs, npackets = npackets)
      if (is.character(aspect)) {
          if (aspect == "xy") {
              aspect <- median(sapply(dxdy, banking) * y.limits$slicelen/x.limits$slicelen, 
                               na.rm = TRUE)
          }
          else if (aspect == "iso") {
              aspect <- median(y.limits$slicelen/x.limits$slicelen, 
                               na.rm = TRUE)
              if (y.relation == "free" || x.relation == "free") 
                  warning("'aspect=\"iso\"' approximate since 'relation=\"free\"'")
          }
          else aspect <- 1
      }
      list(x.limits = x.limits$limits, y.limits = y.limits$limits, 
           x.used.at = x.limits$used.at, y.used.at = y.limits$used.at, 
           x.num.limit = x.limits$numlimitlist, y.num.limit = y.limits$numlimitlist, 
           aspect.ratio = aspect, prepanel.default = prepanel.default, 
           prepanel = prepanel)
  }
  for (packet.number in seq_len(npackets)) {
      id <- compute.packet(cond, cond.current.level)
      foo$packet.sizes[packet.number] <- sum(id)
      foo$panel.args[[packet.number]] <- list(x = x[id], y = y[id])
      foo$panel.args[[packet.number]]$subscripts <- subscr[id]
      cond.current.level <- cupdate(cond.current.level, cond.max.level)
  }
  default.prepanel <- function(x, y, ...){ # formal args x & y are ignored, but provided to avoid R CMD check 'NOTE'
    list(xlim=c(1,length(steps))-c(0.8,0.2),
         ylim=c(0,1))
  }
  prepanel <- NULL # TODO: Can I simply specify 'prepanel', now that I no longer need to adjust panel sizes?
  ## Copied verbatim from lattice, to eliminate ':::' for CRAN check:
  cond.orders <- function (foo, ...) {
      index.cond <- vector(mode = "list", length = length(foo$condlevels))
      for (i in seq_along(foo$condlevels)) index.cond[[i]] <- seq_along(foo$condlevels[[i]])
      perm.cond <- seq_len(length(foo$condlevels))
      if (!is.null(foo$perm.cond)) {
          if (all(sort(foo$perm.cond) == perm.cond)) 
              perm.cond <- foo$perm.cond
          else stop("Invalid value of perm.cond")
      }
      if (!is.null(foo$index.cond)) {
          if (is.list(foo$index.cond) && length(foo$index.cond) == 
              length(index.cond)) {
              for (i in seq_along(foo$condlevels)) index.cond[[i]] <- index.cond[[i]][foo$index.cond[[i]]]
          }
          else if (is.function(foo$index.cond)) {
              FUN <- foo$index.cond
              nplots <- length(foo$panel.args)
              panel.order <- numeric(nplots)
              for (count in seq_len(nplots)) {
                  if (is.list(foo$panel.args[[count]])) {
                      pargs <- c(foo$panel.args.common, foo$panel.args[[count]], 
                                 list(...))
                      prenames <- names(formals(FUN))
                      if (!("..." %in% prenames)) 
                          pargs <- pargs[intersect(names(pargs), prenames)]
                      panel.order[count] <- do.call("FUN", pargs)
                  }
                  else {
                      is.na(panel.order) <- count
                  }
              }
              dim(panel.order) <- sapply(foo$condlevels, length)
              for (i in seq_along(foo$condlevels)) index.cond[[i]] <- order(apply(panel.order, 
                                                                                  i, mean, na.rm = TRUE))
          }
          else stop("Invalid value of index.cond")
      }
      list(index.cond = index.cond, perm.cond = perm.cond)
  }
  more.comp <- c(limits.and.aspect(default.prepanel, prepanel = NULL,
                                   ##have.xlim = TRUE, have.ylim = TRUE,
                                   x.relation = foo$x.scales$relation,
                                   y.relation = foo$y.scales$relation, 
                                   panel.args.common = foo$panel.args.common, panel.args = foo$panel.args, 
                                   aspect = aspect, npackets = npackets, x.axs = foo$x.scales$axs, 
                                   y.axs = foo$y.scales$axs), cond.orders(foo))
  foo[names(more.comp)] <- more.comp
  class(foo) <- "trellis"
  foo
}
