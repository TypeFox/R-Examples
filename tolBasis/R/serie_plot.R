
# usual colors when plotting time-series in TOL
.plot.colors <- function(n) {
  colors <- c("#00C600FF", "blue", "red", "yellow", "orange",
    "#D14949FF", "grey", "purple", "pink", "navyblue", "#356A6AFF",
    "#804000FF", "#05EEFAFF", "#009148FF", "#FFFFBBFF", "#29FB04FF")
  len <- length(colors)
  colors[((0:(n-1)) %% len)+1]
}

.plot.dating <- function(first, last, marks) {
  datings <- list(Yearly, HalfYearly, Quarterly, Monthly, Weekly, Daily)
  dad <- unlist(lapply(datings, function(d) {
    abs(Ddiff(Dceiling(first,d), Dfloor(last,d), d)+1-marks)
  }))
  datings[[ which(dad==min(dad))[[1]] ]]
}

.plot.date.format <- function(dating) {
  if(inherits(dating, "Yearly")) return("%Y")
  if(inherits(dating, c("Monthly", "NMonthly"))) return("%Y-%m")
  if(inherits(dating, "Weekly")) return("%Y-w%U")
  "%Y-%m-%d"
}

plot.Serie.styles = list(
  TOL = list(lab=c(10,5,7), bty="?", las=2, xaxs="i", yaxs="i", lines.lwd=2)
)

plot.Serie <- function(x, y, ...,
  from, to, ylim, dating, date.format, 
  axes=c(T,T,F,F), legend.names, style
) {
  #---------------------------------------------------------------------------
  #[series]
  series <- if(missing(y)) list(x)
    else c(list(x, y), .filter.names(list(...), ""))
  series <- .exclude.null(.flatten(series))
  if(inherits(series, "Serie")) series <- list(series)
  stopifnot(all(sapply(series, function(s) inherits(s, "Serie"))))
  #---------------------------------------------------------------------------
  #[args]
  axes <- c(axes, rep(F, 4))[1:4]
  axes.n <- sum((2^(0:3))[axes])
  bty <- substr("nnnlnnncnnnu7]no", axes.n+1, axes.n+1)
  style.args <- if(missing(style)) list()
  else if(is.null(style)) list()
  else {
    sel <- which(style==names(plot.Serie.styles))
    if(length(sel)>0) plot.Serie.styles[[ sel[[1]] ]]
    else {
      warning("'", style, "' is not a valid style")
      list()
    }
  }
  all.args <- .merge.list(.exclude.names(list(...), ""), style.args)
  if(!is.null(all.args$bty))
    if(all.args$bty=="?") all.args$bty <- bty
  #---------------------------------------------------------------------------
  #[par]
  plot.par = .exclude.prefixes(all.args, c("box", "axes", "lines", 
    "title", "legend")) 
  default.par <- par(rev(as.list(names(plot.par))))
  do.call("par", plot.par)
  #---------------------------------------------------------------------------
  # X
  # x range
  first.date <- do.call(min, lapply(series, Sfirst))
  last.date <- do.call(max, lapply(series, Slast))
  if(missing(from)) from <- first.date
  stopifnot(inherits(from, "Date"))
  stopifnot(from<last.date)  
  if(missing(to)) to <- last.date
  stopifnot(inherits(to, "Date"))
  stopifnot(first.date<to)
  stopifnot(from<to)
  # x affine transformation to [0, 1] (to avoid warnings with EPS)
  .pxmin <- as.numeric(from)
  .pxdif <- as.numeric(to) - .pxmin
  .pxtr <- function(dte) (as.numeric(dte)-.pxmin)/.pxdif
  pxrange <- c(0, 1)
  # dating and x labels
  if(missing(dating)) dating <- .plot.dating(from, to, par("lab")[1]) #see: lab
  stopifnot(inherits(dating, "Dating"))
  dates <- Dseq(from, to, dating)
  pxindices <- .pxtr(dates)
  if(missing(date.format)) date.format <- .plot.date.format(dating)
  pxlabels <- format(dates, date.format)
  # Y
  if(missing(ylim)) {
  pymin <- do.call("min", c(lapply(series, 
    function(s) { min(Ssub(s, from, to), na.rm=TRUE) }), na.rm=TRUE))
  pymax <- do.call("max", c(lapply(series, 
    function(s) { max(Ssub(s, from, to), na.rm=TRUE) }), na.rm=TRUE))
  pyrange <- c(pymin, pymax)
  } else pyrange <- ylim
  #---------------------------------------------------------------------------
  #[plot]
  plot.new()
  plot.window(pxrange, pyrange) #see: xaxs, yaxs, lab
  #---------------------------------------------------------------------------
  #[box & axes]
  #[box]
  box.args <- .merge.list(
    .extract.prefix(all.args, "box"), 
    list(bty=par("bty")),
    .filter.names(.extract.prefix(all.args, "axes"), c("col","lty","lwd")),
    list(lty=par("lty"))
  )
  do.call("box", box.args) #see: bty
  #[axes]
  axes.args <- .merge.list(
    list(at = list(pxindices, NULL), labels = list(pxlabels, TRUE)),
    .extract.prefix(all.args, "axes"),
    par(list("col","lty","lwd"))
  )
  for(i in 1:4) 
    if(axes[i]) do.call("axis", c(list(i), 
      lapply(axes.args, function(v) .extract.cycle(v, i)))) #see: las, lwd
  #[title]
  do.call("title", .extract.prefix(all.args, "title"))
  #---------------------------------------------------------------------------
  #[lines]
  lines.args <- .merge.list(
    .extract.prefix(all.args, "lines"),
    list(col=.plot.colors(length(series)), type="l", lty=par("lty"), pch=par("pch"))
  )
  for(i in length(series):1) {
    x <- .pxtr(Sdates(series[[i]]))
    y <- as.numeric(series[[i]])
    do.call("lines", c(list(x, y), lapply(lines.args, function(v) .extract.cycle(v, i))))
  }
  #---------------------------------------------------------------------------
  #[legend]
  if(missing(legend.names)) legend.names <- NULL
  if(!is.null(legend.names)) {
    # arg: legend.names
    legend.names <- c(legend.names, rep("", length(series)))[1:length(series)]
    default.names <- names(series)
    if(is.null(default.names)) default.names <- rep("", length(series))
    merge.names <- lapply(as.list(1:length(series)), function(i) {
      nme <- legend.names[[i]]
      if(nme=="") nme <- default.names[[i]]
      if(nme=="") nme <- paste0("(", i, ")")
      nme
    })
    # arg: legend.*
    legend.par.names = list("new", "mar", "oma", "fig")
    legend.lty <- unlist(sapply(1:length(series), function(i) { 
      type <- .extract.cycle(lines.args$type, i)
      lty <- if(type %in% c("p","n")) "blank" else .extract.cycle(lines.args$lty, i)
      if(inherits(lty, "numeric")) lty <- c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[lty+1]
      lty
    }))
    legend.pch <- sapply(1:length(series), function(i) { 
      type <- .extract.cycle(lines.args$type, i)
      pch <- if(type %in% c("p","b","o")) .extract.cycle(lines.args$pch, i) else -1
      if(inherits(pch, "character")) pch <- 46 # pch=="."
      pch
    })
    legend.args <- .exclude.names(.merge.list(
      .filter.names(lines.args, list("col", "lwd")),
      list(lty=legend.lty, pch=legend.pch),
      .extract.prefix(all.args, "legend"),
      list(x="bottom", legend=merge.names),
      par(list("col", "lwd", "lty")),
      list(xpd=TRUE, horiz=TRUE, inset=c(0,0), bty="n")
    ), legend.par.names)
    legend.par.args <- .exclude.null(.filter.names(.merge.list(
      list(new=TRUE),
      .extract.prefix(all.args, "legend"),
      list(fig=c(0,1,0,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
    ), legend.par.names))
    previous.par = par(legend.par.names)
    # a new plot is drawn to use all the space avoiding margins
    do.call("par", legend.par.args)
    plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
    do.call("legend", legend.args)
    do.call("par", previous.par)
  }
  #---------------------------------------------------------------------------
  do.call("par", default.par)
  invisible(NULL)
}
