"xysplom" <-
function(x, ...)
UseMethod("xysplom")

"xysplom.formula" <-
function(x, data=sys.parent(),
                            na.action=na.pass, ...) {
  dft <- do.formula.trellis.xysplom(x, data, na.action)
  other <- list(...)
  if (!("xlab" %in% names(list(...)))) other$xlab <- deparse(dft$x.formula[[2]])
  if (!("ylab" %in% names(list(...)))) other$ylab <- deparse(dft$y.formula[[2]])
  do.call("xysplom.default", c(dft[1:3], other))
}

"xysplom.default" <-
function(x, y=x, group, relation="free",
	 x.relation=relation, y.relation=relation,
         xlim.in=NULL, ylim.in=NULL,
	 corr=FALSE, beta=FALSE, abline=corr||beta, digits=3,
	 x.between=NULL, y.between=NULL,
         between.in=list(x=x.between, y=y.between),
         scales.in=list(
           x=list(relation=x.relation, alternating=FALSE),
           y=list(relation=y.relation, alternating=FALSE)),
         strip.in=strip.xysplom,  ## S-Plus requires the function here, not name
	 pch=16, cex=.75,
         panel.input="panel.xysplom", ## character name of function OK in both
         ...,
         cartesian=TRUE,
         plot=TRUE) {

  other <- list(...)
  if (!("xlab" %in% names(list(...)))) other$xlab <- deparse(substitute(x))
  if (!("ylab" %in% names(list(...)))) other$ylab <-
    if (missing(y)) other$xlab
    else deparse(substitute(y))

  if (!is.null(xlim.in)) scales.in$x$limits <- xlim.in
  if (!is.null(ylim.in)) scales.in$y$limits <- ylim.in

  if (is.matrix(x) && !is.null(dimnames(x)[[2]])) {
    dx2 <- dimnames(x)[[2]]
    dx2.done <- TRUE
  }
  else
    dx2.done <- FALSE
  x <- as.data.frame(x)
  n <- nrow(x)

  if (!dx2.done) dx2 <- dimnames(x)[[2]]
  dx2.tmp <- paste(deparse(substitute(x)), seq(length=ncol(x)), sep=".")
  if (length(dx2))
    dx2 <- ifelse(nchar(dx2), dx2, dx2.tmp)
  else
    dx2 <- dx2.tmp

  if (missing(y)) dy2 <- dx2
  else {
    if (is.matrix(y) && !is.null(dimnames(y)[[2]])) {
      dy2 <- dimnames(y)[[2]]
      dy2.done <- TRUE
    }
    else
      dy2.done <- FALSE
    y <- as.data.frame(y)
    if (!dy2.done) dy2 <- dimnames(y)[[2]]
    dy2.tmp <- paste(deparse(substitute(y)), seq(length=ncol(y)), sep=".")
    if (length(dy2))
      dy2 <- ifelse(nchar(dy2), dy2, dy2.tmp)
    else
      dy2 <- dy2.tmp
  }

  y.cn <- rep(dy2, rep(n, ncol(y)))
  x.cn <- rep(dx2, rep(n, ncol(x)))

  old.warn <- options(warn=-1)
  nxy <- n*ncol(x)*ncol(y)
  if (cartesian) ## all the y variables against all the x variables
    ccd <- data.frame(y.list        = unlist(rep(y, rep(ncol(x),ncol(y)))),
                      y             = ordered(as.vector(sapply(
                        rep(as.data.frame(matrix(y.cn,n)),
                            rep(ncol(x),ncol(y))),
                        as.matrix)), dy2),
                      x.list        = rep(unlist(x), length=nxy),
                      x             = ordered(rep(x.cn, length=nxy), dx2),
                      original.row.names = rep(dimnames(x)[[1]], length=nxy))
  else {
    ## each y against the corresponding x
    ## glitch here.  x or y MUST have the same name for all it's levels
    Lx <- length(unique(dx2))
    Ly <- length(unique(dy2))
    if (!(Lx == Ly || Lx==1 || Ly==1)) {
      stop(paste("\nWhen 'cartesian==FALSE' the left-hand side '",
                 paste(dy2, collapse=" + "),
                 "' and right-hand side of the formula '",
                 paste(dx2, collapse=" + "),
                 "' must have the same number of variables."))
    }
    ccd <- data.frame(y.list        = unlist(y),
                      y             = ordered(y.cn, unique(dy2)),
                      x.list        = unlist(x),
                      x             = ordered(x.cn, unique(dx2)),
                      original.row.names     = dimnames(x)[[1]])
    ccd$y <- ordered(paste(as.character(ccd$y), as.character(ccd$x), sep=" ~ "),
                     paste(levels(ccd$y), levels(ccd$x), sep=" ~ "))
    ccd$x <- "x"
  }
  if (missing(group) || is.null(group)) {
    if (cartesian) formula <- y.list ~ x.list | x * y
    else           formula <- y.list ~ x.list |     y
  }
  else {
    group <- interaction(group)
    if (cartesian) ccd$group <- rep(group, length=nxy)
    else           ccd$group <- rep(group, length=nxy/ncol(x))
    if (cartesian) formula <- y.list ~ x.list | x * y * group
    else           formula <- y.list ~ x.list |     y * group
  }
  options(old.warn)


  switch(paste(c("corr", "beta")[c(corr, beta)], collapse="."),
    corr={
      ccd <- cbind(ccd,
                   corr=factor(rep(digits, nrow(ccd))))
      if (missing(group) || is.null(group)) {
        if (cartesian) formula <- y.list ~ x.list | x * y * corr
        else           formula <- y.list ~ x.list |     y * corr
      }
      else {
        if (cartesian) formula <- y.list ~ x.list | x * y * group * corr
        else           formula <- y.list ~ x.list |     y * group * corr
      }
    },
    beta={
       ccd <- cbind(ccd,
                   beta=factor(rep(digits, nrow(ccd))))
      if (missing(group) || is.null(group)) {
        if (cartesian) formula <- y.list ~ x.list | x * y * beta
        else           formula <- y.list ~ x.list |     y * beta
      }
      else {
        if (cartesian) formula <- y.list ~ x.list | x * y * group * beta
        else           formula <- y.list ~ x.list |     y * group * beta
      }
    },
    corr.beta={
       ccd <- cbind(ccd,
                   corr.beta=factor(rep(digits, nrow(ccd))))
      if (missing(group) || is.null(group)) {
        if (cartesian) formula <- y.list ~ x.list | x * y * corr.beta
        else           formula <- y.list ~ x.list |     y * corr.beta
      }
      else {
        if (cartesian) formula <- y.list ~ x.list | x * y * group * corr.beta
        else           formula <- y.list ~ x.list |     y * group * corr.beta
      }
    }
  )
  panel.to.use <-
    if (missing(panel.input) && abline)
      panel=function(x,y,...) {
        panel.xyplot(x,y,...)
        panel.abline(lm(y~x, na.action=na.exclude))
      }
    else panel.input
  if (!cartesian) {
    if.R(r=formals(strip.in)$strip.names <- c(FALSE, FALSE),
         s=strip.in$strip.names <- expression(c(FALSE,FALSE))[[1]])
  }

  result <- list(formula,   ## no name: S-Plus uses "formula", R uses "x"
                 data=ccd,
                 between=between.in,
                 scales=scales.in,
                 panel=panel.to.use,
                 strip=strip.in,
                 pch=pch, cex=cex)
  result <- c(result, other)
  if (plot) do.call("xyplot", result)
  else result
}

"strip.xysplom" <-
function(which.given,
         which.panel,
         var.name,
         factor.levels,
         shingle.intervals,
         par.strip.text=trellis.par.get("add.text"),
         strip.names=c(TRUE,TRUE),
         style=1,
         ...) {
  vnwg <- var.name[which.given]
  if (match(vnwg, c("corr","beta","corr.beta"), 0)) {
###browser()
    ## if.R(r=
         {
           which.parent <- 1
           while(!(exists("rows.per.page", frame=which.parent)))
             which.parent <- which.parent + 1
           cell <- panel.number()
           xy <- get("x",pos=sys.frame(which.parent))$panel.args[[cell]]
           x <- xy$x
           y <- xy$y
         }
         ## ,
         ## s={
         ##   subs <- get("index.list",
         ##               frame=sys.parent())[[get("cell",frame=sys.parent())]]
         ##   x <- get("x",frame=sys.parent())[subs]
         ##   y <- get("y",frame=sys.parent())[subs]
         ## })
    digits <- as.numeric(factor.levels[which.panel[which.given]])
    if (vnwg != "beta") corr <- round(cor(na.exclude(cbind(x,y)))[1,2], digits)
    if (vnwg != "corr") beta <- format(coef(lm(y ~ x, na.action=na.exclude))[2], digits=4)
    strip.names <- c(TRUE,TRUE)
    factor.levels[which.panel[which.given]] <-
      switch(vnwg,
             corr=corr,
             beta=beta,
             corr.beta={
               strip.names <- c(FALSE,FALSE)
               paste("corr: ", corr, "       beta: ", beta, sep="")
             })
  }
  strip.default(which.given=which.given,
                which.panel=which.panel,
                var.name=var.name,
                factor.levels=factor.levels,
                shingle.intervals=shingle.intervals,
                par.strip.text=par.strip.text,
                strip.names=strip.names,
                style=style,
                ...)
}

"panel.xysplom" <-
function(corr, ...) panel.xyplot(...)
