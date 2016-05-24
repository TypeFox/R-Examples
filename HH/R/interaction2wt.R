interaction2wt <- function(x, ...)
  UseMethod("interaction2wt")

"interaction2wt.formula" <-
  function(x, data=sys.parent(), responselab,
           ...) {
    ## if.R(
    ##      r={
    do.formula.trellis <- NA ## make R-2.6.0dev happy
    dft <- do.formula.trellis.xysplom(x, data=data)
    y.in <- dft$y[[1]]
    x.in <- dft$x
    if (missing(responselab)) responselab <- names(dft$y)
    ## },
    ## s={
    ##   dft <- do.formula.trellis(x)
    ##   y.in <- eval(dft$expr[[1]], local=data)
    ##   x.in <- data[,dft$xlab,drop=FALSE]
    ##   if (missing(responselab)) responselab <- dft$ylab
    ## })
    if (is.null(x.in) || is.null(y.in))
      stop("both x and y are needed in formula")
    interaction2wt.default(x=x.in, response.var=y.in,
                           responselab=responselab,
                           ...)
  }

interaction2wt.default <-
  function(x, response.var,
           responselab=deparse(substitute(response.var)),
           responselab.expression = responselab,
           relation=list(x="same", y="same"),
           x.relation=relation$x, y.relation=relation$y,
           digits=3,
           x.between=1,
           y.between=1,
           between,
           cex=.75,
           rot=c(0,0),
           panel.input=panel.interaction2wt,
           strip.input=if (label.as.interaction.formula) strip.default
           else strip.interaction2wt,
           par.strip.text.input=trellis.par.get()$add.text,  ##list(cex=.7)
           scales.additional,
           main.in=paste(responselab,
             ": ", c("main", "simple")[1+simple],
             " effects and 2-way interactions", sep=""),
           xlab="", ylab="",
           simple=FALSE,
           box.ratio=if (simple) .32 else 1,
           label.as.interaction.formula=TRUE,
           ...,
           main.cex,
           key.cex.title=.8,
           key.cex.text=.7,
           factor.expressions=names.x,
           simple.pch=NULL
           ) {
    n <- nrow(x)
    k <- ncol(x)
    names.x <- names(x)
    names(names.x) <- names.x

    if (k<2) stop("interaction2wt requires at least two factors.")
    if (simple && k != 2) stop("Simple effects requires exactly two factors.")

    x.list <- x
    for (i in names(x)) {
      x[[i]] <- as.factor(x[[i]])
      x.list[[i]] <- as.numeric(x[[i]])
    }

    factor.levels <- lapply(x, levels)
    factor.position <- lapply(x, position)
    ##   scales.input <- list(x=list(
    ##                          relation=x.relation,
    ##                          alternating=FALSE,
    ##                          xaxt="n",  ## S-Plus
    ##                          draw=FALSE ## R
    ##                          ),
    ##                        y=list(relation=y.relation, alternating=2))
    ##
    ##
    ##     scales.input <- if.R(r={
    ##       list(x=list(
    ##              relation=x.relation,
    ##              alternating=FALSE,
    ##              draw=FALSE
    ##              ),
    ##            y=list(relation=y.relation, alternating=2))
    ##     }, s={
    ##       list(x=list(
    ##              relation=x.relation,
    ##              alternating=FALSE,
    ##              xaxt="n",
    ##              ),
    ##            y=list(relation=y.relation, alternating=2))
    ##     })
    ##
    xlist <- if.R(r={list(relation=x.relation, alternating=FALSE, draw=FALSE)},
                  s={list(relation=x.relation, alternating=FALSE, xaxt="n"  )})
    scales.input <- list(x=xlist,
                         y=list(relation=y.relation, alternating=2))


    if (!missing(scales.additional)) {
      scales.input$x[names(scales.additional$x)] <- scales.additional$x
      scales.input$y[names(scales.additional$y)] <- scales.additional$y
    }
    if.R(r={
      scales.input$x$at <- NULL
      scales.input$y$at <- NULL
      scales.input$x$rot <- rot[1]
      scales.input$y$rot <- rep(rot,2)[2]
    },
         s={})

    ccd <- data.frame(response.var=rep(response.var, length=n*k*k),
                      x.values    =unlist(rep(as.list(x.list), k)),
                      trace.values=unlist(rep(as.list(x.list), rep(k,k))),
                      x.factor    =factor(rep(rep(names.x, rep(n,k)), k),
                        levels=names.x),
                      trace.factor=factor(rep(    names.x, rep(n*k,   k)),
                        levels=names.x))
    if (label.as.interaction.formula) {
      ccd$x.trace <- interaction(ccd$x.factor, ccd$trace.factor)
      levels(ccd$x.trace) <- paste(responselab,
                                   outer(levels(ccd$x.factor),
                                         levels(ccd$trace.factor),
                                         FUN=paste,
                                         sep=" | "),
                                   sep=" ~ ")
      formula <- response.var ~ x.values | x.trace
    }
    else
      formula <- response.var ~ x.values | x.factor * trace.factor

    if (!missing(main.cex)) {
      main.in <- as.list(main.in)
      main.in$cex <- main.cex
    }

    if (is.null(simple.pch)) {
      simple.pch=lapply(x, function(xi) seq(along=levels(xi)))
    }

    xyplot.list <-
      list(formula,
           data=ccd,
           responselab=responselab,
           responselab.expression=responselab.expression,
           trace.values=ccd$trace.values,
           factor.levels=factor.levels,
           factor.position=factor.position,
           between=if (missing(between))
           list(x=x.between, y=y.between)
           else
           between,
           scales=scales.input,
           xaxs="e",
           prepanel=function(x,y) list(xlim=range(factor.position)+c(-.5,.5)), ##range(x)+c(-1,1)*.1*range(x)),
           panel=panel.input,
           strip=strip.input,
           par.strip.text=par.strip.text.input,
           layout=c(k, k),
           main=main.in,
           xlab="", ylab="",
           cex=cex, las=1, aspect=1,
           simple=simple,
           data.x=x,
           box.ratio=box.ratio,
           simple.pch=simple.pch,
           ...)
    if.R(r={
      cpy <- range(ccd$response.var, na.rm=TRUE)
      pcpy <- pretty(cpy)
      pcpy <- pcpy[(cpy[1] <= pcpy) & (pcpy <= cpy[2])]
      ## recover()

      lattice.options <-
        list(axis.options=list(
               bottom=list(
                 at2=factor.position,
                 labels2=factor.levels,
                 rot2=rot[1],
                 labels3=factor.expressions), ## levels(ccd$x.factor)),
               right=list(
                 at2=pcpy,
                 labels2=pcpy,
                 rot2=rot[2],
                 labels3=rep(responselab.expression, k)
                 )
               ),
             layout.heights=list(axis.xlab.padding=list(x=15, units="mm")),
             layout.widths=list(right.padding=list(x=13, units="mm")))

      if (length(xyplot.list$lattice.options) == 0)
        xyplot.list$lattice.options <- list()
      xyplot.list$lattice.options[names(lattice.options)] <- lattice.options

      keys <- vector("list")
      for (ii in seq(along=names.x)) {
        trace.id <- names(x)[ii]
        keylist <- list(title=factor.expressions[trace.id],
                        cex.title=key.cex.title,
                        border=TRUE,
                        text=list(
                          text=factor.levels[[trace.id]],
                          cex=key.cex.text),
                        lines=Rows(
                          trellis.par.get("superpose.line"),
                          seq(length(factor.levels[[trace.id]]))))
        pssl.lwd <- list(...)$par.settings$superpose.line$lwd
        if (!is.null(pssl.lwd)) keylist$lines$lwd[] <- pssl.lwd
        keys[[trace.id]] <- draw.key(keylist, draw=FALSE)

        if (simple) {
          other.id <- names(x)[3 - ii]
          keylist <- list(title=factor.expressions[other.id],
                          cex.title=key.cex.title,
                          border=TRUE,
                          text=list(
                            text=factor.levels[[other.id]],
                            cex=key.cex.text),
                          points=Rows(
                            trellis.par.get("superpose.symbol"),
                            seq(length(factor.levels[[other.id]]))))
          keylist$points$col[] <- "black"
          keylist$points$pch <- simple.pch[[other.id]]
          keylist$lines=list(col=0, size=2)
          keys[[paste(other.id, "pts", sep="")]] <- draw.key(keylist, draw=FALSE)

          ## other.id <- names(x)[3 - ii]
          ## keylist$title=factor.expressions[other.id]
          ## keylist$points <- Rows(
          ##   trellis.par.get("superpose.symbol"),
          ##   seq(length(factor.levels[[trace.id]])))
          ## keylist$points$col[] <- "black"
          ## keylist$points$pch <- simple.pch[[other.id]]
          ## keylist$lines <- NULL
          ## keylist$text=list(
          ##   text=factor.levels[[other.id]],
          ##   cex=key.cex.text)
          ## keys[[paste(other.id, "pts", sep="")]] <- draw.key(keylist, draw=FALSE)

        }
      }

      xyplot.list$legend <- list(left =
                                 list(fun = legendGrob2wt,
                                      args = keys))
      xyplot.list$axis <- axis.i2wt
      ## recover()
    },
         s={})
    do.call("xyplot", xyplot.list)
  }

"strip.interaction2wt" <-
  function(which.given,
           which.panel,
           var.name,
           factor.levels,
           shingle.intervals,
           strip.names=c(TRUE,TRUE),
           style=1,
           ...) {
    strip.default(which.given=which.given,
                  which.panel=which.panel,
                  var.name=var.name,
                  factor.levels=factor.levels,
                  shingle.intervals=shingle.intervals,
                  strip.names=strip.names,
                  style=style,
                  ...)
  }

## source("c:/HOME/rmh/HH-R.package/HH/R/interaction2wt.R")
