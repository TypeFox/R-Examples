## general.R --- File for generally useful functions.
## Author: Stephen J Eglen (2 functions from Ben Bolker's bbmisc package, GPL'ed)
## Copyright: GPL

.file.or.gz <- function(file) {
  ## Return FILE if it exists, or FILE.gz if that exists.
  ## Otherwise, return NA and generate a warning.
  if (file.exists(file)) {
    file
  } else {
    f2 <- paste(file,".gz", sep="")
    if (file.exists(f2))
      f2
    else {
      warning(paste("File", file,
                    "could not be found, nor its compressed version."))
      NA
    }
  }
}

## .plotCI() is taken from Ben Bolker's bbmisc package, which is under the GPL.
## http://www.zoo.ufl.edu/bolker/R/src/
## Update Tue 28 Nov 2006, now need .clean.args also.

## remove arguments not intended for a particular function from a string
## repeated from bbfuns/misc.R
.clean.args <- function(argstr,fn,exclude.repeats=FALSE,
                       exclude.other=NULL,dots.ok=TRUE) {
  fnargs <- names(formals(fn))
  if (length(argstr)>0 && !("..." %in% fnargs && dots.ok))  {
    badargs <- names(argstr)[!sapply(names(argstr),"%in%",c(fnargs,""))]
    for (i in badargs)
      argstr[[i]] <- NULL
  }
  if (exclude.repeats) {
    ntab <- table(names(argstr))
    badargs <- names(ntab)[ntab>1 & names(ntab)!=""]
    for (i in badargs)
      argstr[[i]] <- NULL
  }
  for (i in exclude.other)  ## additional arguments to exclude.other
    argstr[[i]] <- NULL
  argstr
}


.plotCI <- function (x, y = NULL, uiw, liw = uiw,
                    ui=NULL, li=NULL,
                    err="y",
                    sfrac = 0.01, gap=0, slty=par("lty"),
                    add=FALSE,
                    scol=NULL,
                    pt.bg=par("bg"),
                    ...)  {
  ## from Bill Venables, R-list, modified with contributions and ideas
  ## from Gregory Warnes and the list
  ## requires .clean.args()
  ## process arguments:
  arglist <- list(...)
  if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (is.null(y)) {
    if (is.null(x)) 
      stop("both x and y NULL")
    y <- as.numeric(x)
    x <- seq(along = x)
  }
  if (missing(uiw) && (is.null(ui) || is.null(li)))
    stop("must specify either relative limits or both lower and upper limits")
  if (!missing(uiw)) {  ## 
    if (err=="y") z <- y else z <- x
    ui <- z + uiw
    li <- z - liw
  }
  ## fill in default arguments
  if (is.null(arglist$"xlab"))
    arglist$"xlab" <- deparse(substitute(x))
  if (is.null(arglist$"ylab"))
    arglist$"ylab" <- deparse(substitute(y))
  if (err=="y" && is.null(arglist$"ylim"))
    arglist$"ylim" <- range(c(y, ui, li), na.rm=TRUE)
  if (err=="x" && is.null(arglist$"xlim"))
    arglist$"xlim" <- range(c(x, ui, li), na.rm=TRUE)
  if (missing(scol)) {
    if (!is.null(arglist$"col")) scol <- arglist$"col"
    else scol <- par("col")
  }
  plotpoints <- TRUE
  if (!is.null(arglist$"pch") && is.na(arglist$"pch")) {
    arglist$"pch" <- 1
    plotpoints <- FALSE
  }
  ## 
  if (!add)
    do.call("plot",c(list(x,y,type="n"),
                          .clean.args(arglist,plot)))
  if (gap==TRUE) gap <- 0.01  ## default gap size: maybe too complicated?
  ul <- c(li, ui)
  if (err=="y") {
    gap <- rep(gap,length(x))*diff(par("usr")[3:4])
    smidge <- par("fin")[1] * sfrac
    arrow.args <- c(list(lty=slty,angle=90,length=smidge,code=1,
                         col=scol),
                    .clean.args(arglist,arrows,
                               ## sje --- add "type"?
                               exclude.other=c("col","lty", "type")))
    do.call("arrows",c(list(x , li, x, pmax(y-gap,li)),
                       arrow.args))
    do.call("arrows",c(list(x , ui, x, pmin(y+gap,ui)),
                       arrow.args))
  }
  else if (err=="x") {
    gap <- rep(gap,length(x))*diff(par("usr")[1:2])
    smidge <- par("fin")[2] * sfrac
    arrow.args <- c(list(lty=slty,angle=90,length=smidge,code=1),
                    .clean.args(arglist,arrows,exclude.other=c("col","lty")))
    do.call("arrows",c(list(li, y, pmax(x-gap,li), y),
                       arrow.args))
    do.call("arrows",c(list(ui, y, pmin(x+gap,ui), y),
                       arrow.args))
  }
  ## now draw the points (in case we want to have "bg" set for points)
  if (plotpoints)
    do.call("points",c(list(x, y, bg=pt.bg),
                       .clean.args(arglist,points,
                                  exclude.other=c("xlab","ylab","xlim","ylim",
                                    "axes"))))
  invisible(list(x = x, y = y))
}

.printf <- function (...) {
  ## Helper function.
  cat(sprintf(...))
}
