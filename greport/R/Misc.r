utils::globalVariables(c('Freq', '.group.'))


#' Compute mfrow Parameter
#'
#' Compute a good \code{par("mfrow")} given the
#' number of figures to plot.
#'
#' @param n numeric. Total number of figures to place in layout.
#' @param small logical. Set to \sQuote{TRUE} if the plot area should be
#' smaller to accomodate many plots.
#' @return return numeric vector.
#' oldmfrow <- mfrowSet(8)

mfrowSuggest <- function(n, small=FALSE) {
  omf <- mf <- par('mfrow')
  if(length(mf) == 0) mf <- c(1,1)
  if(n == 1) return(mf)
  if(n > 1 & max(mf) == 1) {
    if(small) {
      mf <- if(n <= 2) {
        c(1, 2)
      } else if(n <= 4) {
        c(2,2)
      } else if(n <= 6) {
        c(2,3)
      } else if(n <= 12) {
        c(3,4)
      } else if(n <= 16) {
        c(4,4)
      } else if(n <= 20) {
        c(4,5)
      } else if(n <= 24) {
        c(4,6)
      } else if(n <= 25) {
        c(5,5)
      } else if(n <= 30) {
        c(5,6)
      } else if(n <= 36) {
        c(6,6)
      } else if(n <= 42) {
        c(6,7)
      } else {
        c(6,8)
      }
    } else {
      mf <- if(n <= 2) {
        c(1,2)
      } else if(n <= 4) {
        c(2,2)
      } else if(n <= 6) {
        c(2,3)
      } else if(n <= 9) {
        c(3,3)
      } else {
        c(4,3)
      }

      if(n > 12 & n <= 16) {
        mf <- c(4,4)
      }
    }
    }
  mf
}

#' Set greport Options
#'
#' @param \dots a series of options for which non-default values are desired:
#' \itemize{
#'  \item{\code{tx.pch}:}{symbols corresponding to treatments}
#'  \item{\code{tx.col}:}{colors corresponding to treatments}
#'  \item{\code{tx.linecol}:}{colors for lines in line plots}
#'  \item{\code{nontx.col}:}{colors for categories other than treatments}
#'  \item{\code{tx.lty}:}{line types corresponding to treatments}
#'  \item{\code{tx.lwd}:}{line widths corresponding to treatments}
#'  \item{\code{tx.var}:}{character string name of treatment variable}
#'  \item{\code{er.col}:}{2-vector with names \code{"enrolled","randomized"} containing colors to use for enrolled and randomized in needle displays}
#'  \item{\code{alpha.f}:}{single numeric specifying alpha adjustment to be applied to all colors.  Default is 0.7}
#'  \item{\code{denom}:}{named vector with overall sample sizes}
#'  \item{\code{tablelink}:}{a character string, either \code{"tooltip"} or \code{"hyperref"} (the default); use the latter to make supporting data tables be hyperlinked to tables in the appendix rather than using a pop-up tooltip}
#' \item{\code{figenv}:}{LaTeX figure environment to use, default is \code{"figure"}.  Use \code{figenv="SCfigure"} if using the LaTeX \code{sidecap} package.  \code{SCfigure} is only used for narrow images, by calling \code{putFig} with \code{sidecap=TRUE}.}
#' \item{code{figpos}:}{LaTeX figure environment position; default is \code{"htb!"}}
#'  \item{\code{gtype}:}{graphics type, \code{"pdf"} or \code{"interactive"}}
#'  \item{\code{pdfdir}:}{name of subdirectory in which to write \code{pdf} graphics}
#'  \item{\code{texdir}:}{name of subdirectory in which to write \code{LaTeX} code}
#'  \item{\code{texwhere}:}{default is \code{"texdir"} to use location specified by \code{texdir}.  Set to \code{""} to write generated non-appendix LaTeX code to the console as expected by \code{knitr}}
#'	\item{\code{defs}"}{fully qualified file name to which to write LaTeX macro definitions such as poptables}
#' }
setgreportOption <- function(...) {
  default <- getOption('greport')
  opts <- list(...)
  alpha.f <- if(length(opts$alpha.f)) opts$alpha.f else 0.7
  if(! length(default))
    default <-
      list(tx.pch = 16:17,
           tx.col = adjustcolor(c('black', '#0080ff'), alpha.f=alpha.f),
           tx.linecol = adjustcolor(c('black', '#0080ff'), alpha.f=alpha.f),
           nontx.col = adjustcolor(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
             "#66a61e", "#e6ab02", "#a6761d", "#666666"),
             alpha.f=alpha.f),  ## see colorbrewer2.org
           tx.lty = c(1, 1), tx.lwd = c(1, 2),
           tx.var = '', er.col = NULL, alpha.f = 0.7,
           denom = c(enrolled=NA, randomized=NA),
           tablelink = 'hyperref', figenv='figure', figpos='htb!',
           gtype = 'pdf', pdfdir='pdf', texdir='gentex', 
           texwhere='texdir', defs='gentex/defs.tex')

  if(length(opts)) {
    if(any(names(opts) %nin% names(default)))
      stop(paste('greport options must be one of the following:',
                 paste(names(default), collapse=' ')))
    default[names(opts)] <- opts
  }
  i <- names(opts$denom) %nin% c('enrolled', 'randomized')
  if(any(i) && sum(opts$denom[i]) != opts$denom['randomized'])
    stop('sum of # subjects randomized to each treatment must = total number randomized')
  if(! length(default$er.col))
    default$er.col <-
      adjustcolor(setdiff(c('red', 'green', "#0080ff", "#ff00ff",
                            "darkgreen", "#ff0000", "orange", "#00ff00",
                            "brown"),
                          default$tx.col)[1 : 2], alpha.f=alpha.f)
  options(greport = default)
  invisible()
}

#' Setup lattice plots using greport options
#'
#' Initializes colors and other graphical attributes based on
#' what is stored in system option \code{greport}.
latticeInit <- function() {
  opt <- getOption('greport')

  tx.col <- opt$tx.col
  nontx.col <- opt$nontx.col
  tx.pch <- opt$tx.pch
  tx.lty <- opt$tx.pty
  tx.lwd <- opt$tx.lwd
  alpha.f <- opt$alpha.f
  w <- trellis.par.get('dot.symbol')
  w$col <- tx.col[1]
  w$pch <- tx.pch[1]
  trellis.par.set(dot.symbol=w)
  w <- trellis.par.get('superpose.symbol')
  # w$col <- adjustcolor(w$col, alpha.f=alpha.f)
  # w$col <- rep(c(tx.col, setdiff(w$col, tx.col)), length=length(w$col))
  w$col <- rep(c(tx.col, nontx.col), length=length(w$col))
  w$pch <- rep(c(tx.pch,
                 setdiff(c(1, 2, 3, 4, 5, 6, 16, 17, 15, 18, 20), tx.pch)),
               length=length(w$pch))
  trellis.par.set(superpose.symbol=w)
  w <- trellis.par.get('superpose.line')
  # w$col <- adjustcolor(w$col, alpha.f=alpha.f)
  w$col <- rep(c(tx.col, nontx.col), length=length(w$col))
  w$lty <- rep(c(tx.lty, w$lty), length=length(w$lty))
  w$lwd <- rep(c(tx.lwd, w$lwd), length=length(w$lwd))
  trellis.par.set(superpose.line=w)
  w <- trellis.par.get('plot.symbol')
  w$col <- tx.col[1]
  w$pch <- tx.pch[1]
  trellis.par.set(plot.symbol=w)

  ## Some of the following seem to have no effect
  trellis.par.set(
    axis.text       = list(cex = 0.75),
    par.xlab.text   = list(cex = 0.85),
    par.ylab.text   = list(cex = 0.85),
    par.strip       = list(cex = 0.8 ),
    layout.heights  = list(strip = 0.9, main = 0, sub  = 0),
    scales = list(y = list(rot = 0)))

  lattice.options(default.args = list(as.table = TRUE))
  
  invisible(opt)
}

#'  
#' Get greport Options
#'
#' Get greport options, assigning default values of unspecified optios.
#'
#' @param opts character vector containing list of option names to retrieve.  If only one element, the result is a scalar, otherwise a list.  If \code{opts} is not specified, a list with all current option settings is returned.
#' @export

getgreportOption <- function(opts=NULL) {
  go <- getOption('greport')
  if(! length(opts)) return(go)
  go <- go[opts]
  if(length(opts) == 1) go <- go[[1]]
  go
}

#' Compute Sample Fractions
#'
#' Uses denominators stored with \code{setgreportOption} along with counts specified to \code{SampleFrac} to compute fractions of subjects in current analysis
#'
#' @param n integer vector, named with \code{"enrolled","randomized"} and optionally also including treatment levels.
#' @param nobsY a result of the the \code{nobsY} Hmisc function
#' @param table set to \code{TRUE} to return as an attribute \code{"table"} a character string containing a LaTeX tabular showing the pertinent frequencies created from \code{n} and the \code{denom} option, and if \code{nobsY} is present, adding another table with response variable-specific counts.
#' @export

sampleFrac <- function(n, nobsY=NULL, table=TRUE) {
  denom <- getgreportOption('denom')
  if(any(is.na(denom))) stop('denom must be defined with setgreportOption()')
  if(names(n)[1] != 'enrolled')
    n <- structure(c(n[1], n), names=c('enrolled', names(n)))
  if(all(names(n) %in% c('enrolled', 'randomized')))
    denom <- denom[unique(c('enrolled', names(n)))]
  if(length(n) != length(denom))
    stop('length of n does not equal length of denom')
  if(any(names(n) != names(denom)))
    stop('n does not have same names as denom in the same order')
  f <- n / denom
  if(any(f > 1.)) warning('A sample fraction > 1.0; assuming analysis is to compare randomized and non-randomized subjects\nfraction capped at 1.0')
  f <- pmin(f, 1.)
  if(! table) return(f)
  tab <- data.frame(upFirst(names(n)), denom, n)
  tab <- latexTabular(tab, align = 'lrr', translate=FALSE, hline=1,
                      headings=c('Category', '$N$', 'Used in Analysis'))
  if(length(nobsY)) {
    tab <- paste(tab, '\n\\vspace{1ex}\n\\hsepline\n\\vspace{1ex}\n')
    if(length(m <- nobsY$nobsg)) {
      m <- t(m)
      d <- cbind(Variable=rownames(m), as.data.frame(m))
      tab2 <- latexTabular(d, align=paste('l',
                                paste(rep('r', ncol(m)), collapse=''),
                                sep=''),
                           translate=FALSE, hline=1)
    }
    else {
      m <- nobsY$nobs
      d <- data.frame(Variable=names(m), N=m)
      tab2 <- latexTabular(d, align='lr', translate=FALSE, hline=1,
                           headings=c('Variable', '$N$'))
    }
    tab <- paste(tab, tab2)
  }
  attr(f, 'table') <- tab
  f
}

#' Draw Needles
#'
#' Create a LaTeX \code{picture} to draw needles for current sample sizes.  Uses colors set by call to \code{setgreportOptions}.
#'
#' @param sf output of \code{sampleFrac}
#' @param name character string name of LaTeX variable to create
#' @param file output file name (character string)
#' @param append set to \code{FALSE} to start a new \code{file}
#' @export

dNeedle <- function(sf, name, file='', append=TRUE) {
  co <- getgreportOption(c('er.col', 'tx.col'))
  co <- c(co$er.col, co$tx.col)
  latexNeedle(sf, col=co, name=name, file=file, append=append)
}


#' Put Figure
#'
#' Included a generated figure within LaTex document.  \code{tcaption} and \code{tlongcaption} only apply if \code{setgreportOption(tablelink="hyperref")}.
#'
#' @param panel character. Panel name.
#' @param name character. Name for figure.
#' @param caption character. Short caption for figure.
#' @param longcaption character. Long caption for figure.
#' @param tcaption character.  Short caption for supporting table.
#' @param tlongcaption character.  Long caption for supporting table.
#' @param poptable an optional character string containing LaTeX code that will be used as a pop-up tool tip for the figure (typically a tabular).  Set to \code{NULL} to suppress supplemental tables that back up figures.
#' @param popfull set to \code{TRUE} to make the pop-up be full-page
#' @param sidecap set to \code{TRUE} (only applies if \code{greportOption(figenv="SCfigure")}) to assume the figure is narrow and to use side captions
#' @param outtable set to \code{TRUE} to only have the caption and hyperlink to graphics in a LaTeX table environment and to leave the tabulars to free-standing LaTeX markup.  This is useful when the table is long, to prevent hyperlinking from making the table run outside the visable area.  Instead of the hyperlink area being the whole table, it will be the caption.  A \code{clearpage} is issued after the tabular.
#' @param append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @export

putFig <- function(panel, name, caption=NULL, longcaption=NULL,
                   tcaption=caption, tlongcaption=NULL,
                   poptable=NULL, popfull=FALSE, sidecap=FALSE,
                   outtable=FALSE, append=TRUE) {
  o <- getgreportOption()
  gtype     <- o$gtype
  texdir    <- o$texdir
  texwhere  <- o$texwhere
  tablelink <- o$tablelink
  figenv    <- o$figenv
  figpos    <- o$figpos
  if(! sidecap) figenv <- 'figure'

  if(length(gtype) && gtype == 'interactive') return(invisible())
  
  bcenter <- if(figenv == 'figure') '\\centerline{' else ''
  ecenter <- if(figenv == 'figure') '}' else ''
  
  panel <- gsub('\\.', '-', panel)
  name  <- gsub('\\.', '-', name)
  file  <- sprintf('%s/%s.tex', texdir, panel)
  if(texwhere == '') file <- ''

  ## if(length(caption)) caption <- latexTranslate(caption)
  ## if(length(longcaption)) longcaption <- latexTranslate(longcaption)

  sf <- function(...) paste(sprintf(...), '%\n', sep='')
  cap <- lab <- tlab <- ""
  if(length(longcaption))
    cap <- sf("\\caption[%s]{%s", caption, longcaption)
  else
    if(length(caption)) cap <- sf("\\caption{%s", caption)
  else
    cap <- '\\caption{'
  
  if(length(caption)) {
    lab  <- sf("\\label{fig:%s}", name)
    tlab <- sf("\\label{tab:%s}", name)
  }

  if(! length(poptable)) {
    cap <- paste(cap, '}', sep='')
    cat(sf("\\begin{%s}[%s]\\leavevmode%s\\includegraphics{%s.pdf}%s%s\n%s%s\\end{%s}\n", figenv, figpos, bcenter, name, ecenter, '%', cap, lab, figenv),
        file=file, append=append)
    return()
  }

  if(tablelink == 'tooltip') {
    cmd <- if(popfull) '\\tooltipw' else '\\tooltipm'
    cat(sf("\\begin{%s}[%s]\\leavevmode%s\\protect%s{\\includegraphics{%s.pdf}%s{%s}}%s%s\\end{%s}",
           figenv, figpos, bcenter, cmd, name, ecenter, poptable, cap,
           lab, figenv),
        file=file, append=append)
  } else {
    ref <- if(length(caption))
      sprintf(' {\\smaller (Figure \\ref{fig:%s})}.', name)
    else ''
    ## linkfromtab <- if(outtable)
    ##    sf('\\hyperref[fig:%s]{~\\textcolor[gray]{0.5}{$\\mapsto$}}', name)
    ##  else ''
    tcap <- if(length(tlongcaption))
      sf('\\caption[%s]{%s%s}', tcaption, tlongcaption, ref)
    else if(length(tcaption)) sf('\\caption[%s]{%s%s}', tcaption, tcaption,
                                 ref)
    else sprintf('\\caption{%s}', ref)
    cat(sf('\\begin{%s}[%s]\\hyperref[tab:%s]{\\leavevmode%s\\includegraphics{%s.pdf}%s}', figenv, figpos, name, bcenter, name, ecenter),
        file=file, append=append)

    reft <- sprintf(' {\\smaller (Table \\ref{tab:%s})}}', name)
    cap <- if(grepl('\\hfill', cap))
      gsub('\\hfill', paste(reft, '\\hfill', sep=''), cap, fixed=TRUE)
    else
      paste(cap, reft, sep='')

    cat(cap, lab, 
        sprintf('\\end{%s}\n', figenv), sep='',
        file=file, append=TRUE)
    appfile <- sprintf('%s/app.tex', texdir)

    if(outtable)
      cat(sf('\\clearpage\\begin{table}[%s]%s%s\\end{table}\n%s\\clearpage\n\n',
           figpos, tcap, tlab, poptable),
        file=appfile, append=TRUE)
    else
      cat(sf('\\begin{table}[%s]%s%s\\hyperref[fig:%s]{%s}\\end{table}\\clearpage\n\n',
             figpos, tcap, tlab, name, poptable),
          file=appfile, append=TRUE)
  }
  invisible()
}

#' Plot Initialization
#'
#' Toggle plotting.  Sets options by examining \code{setgreportOption(gtype=)}.
#'
#' @param file character.  Character string specifying file prefix.
#' @param h numeric.  Height of plot in inches, default=7.
#' @param w numeric.  Width of plot in inches, default=7.
#' @param lattice logical.  Set to \code{FALSE} to prevent \code{latticeInit} from being called.
#' @param \dots Arguments to be passed to \code{spar}.
#' @export

startPlot <- function(file, h=7, w=7, lattice=TRUE, ...) {
  gtype  <- getgreportOption('gtype')
  pdfdir <- getgreportOption('pdfdir')
  if(! length(gtype) || gtype != 'interactive') {
    file <- paste(pdfdir, '/', gsub('\\.', '-', file), '.pdf', sep='')
    pdf(file, height=h, width=w)
  }
  if(! existsFunction('spar'))
    spar <-
      function(mar=if(!axes)
               c(2.25+bot-.45*multi,2+left,.5+top+.25*multi,.5+rt) else
               c(3.25+bot-.45*multi,3.5+left,.5+top+.25*multi,.5+rt),
               lwd = if(multi)1 else 1.75,
               mgp = if(!axes) mgp=c(.75, .1, 0) else
               if(multi) c(1.5, .365, 0) else c(2.4-.4, 0.475, 0),
               tcl = if(multi)-0.25 else -0.4, xpd=FALSE,
               bot=0, left=0, top=0, rt=0, ps=if(multi) 14 else 10,
               mfrow=NULL, axes=TRUE, cex.lab=1.25, cex.axis=1.15,
               ...) {
        multi <- length(mfrow) > 0
        par(mar=mar, lwd=lwd, mgp=mgp, tcl=tcl, ps=ps, xpd=xpd,
            cex.lab=cex.lab, cex.axis=cex.axis, ...)
        if(multi) par(mfrow=mfrow)
      }
  dotlist <- list(...)
  if(length(dotlist)) {
    allow <- union(names(par()), setdiff(names(formals(spar)), '...'))
    dotlist <- dotlist[names(dotlist) %in% allow]
  }
  do.call('spar', dotlist)
  if(lattice) latticeInit()
  invisible()
}

#' @rdname startPlot
#' @export

endPlot <- function() {
  gtype <- getgreportOption('gtype')
  if(length(gtype) && gtype == 'interactive') return(invisible())
  dev.off()
  invisible()
}

#' Issue LaTeX section and/or subsection in appendix
#'
#' This is useful for copying section and subsection titles in the main body of the report to the appendix, to help in navigating supporting tables.  LaTeX backslash characters need to be doubled.
#'
#' @param section a character string that will cause a section command to be added to app.tex
#' @param subsection a character string that will cause a subsection command to be added to app.tex
#' @param main set to \code{TRUE} to also write a section or subsection command to the console to be used in building the main report body (graphical section), in which case you should also specify \code{panel} if option \code{texdir} is not an empty string
#' @param panel panel string; must be given if \code{main=TRUE} and option \code{texdir} is not \code{""}
#' @export

appsection <- function(section=NULL, subsection=NULL, main=FALSE, panel='') {
  o <- getgreportOption()
  texdir   <- o$texdir
  if(main) {
    texwhere <- o$texwhere
    file <- if(texwhere == '') '' else sprintf('%s/%s.tex', texdir, panel)
    if(length(section))    cat('\\section{', section, '}\n',
                               sep='', file=file, append=TRUE)
    if(length(subsection)) cat('\\subsection{', subsection, '}\n',
                               sep='', file=file, append=TRUE)
  }
  file  <- sprintf('%s/app.tex', texdir)
  if(length(section))    cat('\\section{', section, '}\n',
                             sep='', file=file, append=TRUE)
  if(length(subsection)) cat('\\subsection{', subsection, '}\n',
                             sep='', file=file, append=TRUE)
  invisible()
}

#' Merge Multiple Data Frames or Data Tables
#'
#' Merges an arbitrarily large series of data frames or data tables containing common \code{id} variables (keys for data tables).  Information about number of observations and number of unique \code{id}s in individual and final merged datasets is printed.  The first data frame has special meaning in that all of its observations are kept whether they match \code{id}s in other data frames or not.  For all other data frames, by default non-matching observations are dropped.  The first data frame is also the one against which counts of unique \code{id}s are compared.  Sometimes \code{merge} drops variable attributes such as \code{labels} and \code{units}.  These are restored by \code{Merge}.  If all objects are of class \code{data.table}, faster merging will be done using the \code{data.table} package's join operation.  This assumes that all objects have identical key variables and those of the variables on which to merge.
#'
#' @param \dots two or more dataframes or data tables
#' @param id a formula containing all the identification variables such that the combination of these variables uniquely identifies subjects or records of interest.  May be omitted for data tables; in that case the \code{key} function retrieves the id variables.
#' @param all set to \code{FALSE} to drop observations not found in second and later data frames (only applies if not using \code{data.table})
#' @param verbose set to \code{FALSE} to not print information about observations
#' @export
#' @examples
#' a <- data.frame(sid=1:3, age=c(20,30,40))
#' b <- data.frame(sid=c(1,2,2), bp=c(120,130,140))
#' d <- data.frame(sid=c(1,3,4), wt=c(170,180,190))
#' all <- Merge(a, b, d, id = ~ sid)
#' # For data.table, first file must be the master file and must
#' # contain all ids that ever occur.  ids not in the master will
#' # not be merged from other datasets.
#' a <- data.table(a); setkey(a, sid)
#' # data.table also does not allow duplicates without allow.cartesian=TRUE
#' b <- data.table(sid=1:2, bp=c(120,130)); setkey(b, sid)
#' d <- data.table(d); setkey(d, sid)
#' all <- Merge(a, b, d)

Merge <- function(..., id, all=TRUE, verbose=TRUE) {
  w <- list(...)
  nams <- (as.character(sys.call())[-1])[1 : length(w)]
  m <- length(nams)
  ## If argument is a function call, e.g., subset(mydata, age > 20)
  ## find name of first argument and omit any dollar sign prefix and []
  for(i in 1 : m) {
    x <-       nams[i]
    x <-       gsub('subset\\(',   '', x)
    x <-       gsub(',.*',         '', x)
    x <-       gsub('\\[.*'  ,     '', x)
    nams[i] <- gsub('(.*)\\$(.*)', '\\2', x)
  }
  d1   <- w[[1]]
  idt  <- is.data.table(d1)
  id   <- if(idt) key(d1) else all.vars(id)
  m <- length(w)
  va <- n <- nu <- integer(m)
  nin1 <- nnin1 <- rep(NA, m)
  did <- if(idt) d1[, id, with=FALSE] else d1[id]
  idc1 <- unique(as.character(interaction(did)))
  id.union <- id.intersection <- idc1
  ## Unique variables, and their labels and units
  uvar <- lab <- un <- character(0)
  for(i in 1 : m) {
    d <- w[[i]]
    nd <- names(d)
    if(any(id %nin% nd))
      stop(paste('data frame', nams[i], 'does not contain id variables',
                 paste(id, collapse=', ')))
    j <- nd %nin% uvar
    uvar <- c(uvar, nd[j])
    lab  <- c(lab,  sapply(d, label)[j])
    un   <- c(un,   sapply(d, units)[j])
    idt <- is.data.table(d)
    M <- if(i == 1) d
    else
      if(idt) d[M]
    else
      merge(M, d, by=id, all.x=TRUE, all.y=all)
    did <- if(idt) d[, id, with=FALSE] else d[id]
    idc <- unique(as.character(interaction(did)))
    di <- dim(d)
    va[i] <- di[2]
    n [i] <- di[1]
    nu[i] <- length(unique(idc))
    if(i > 1) {
      nin1 [i] <- sum(idc %in%  idc1)
      nnin1[i] <- sum(idc %nin% idc1)
      id.union <- union(id.union, idc)
      id.intersection <- intersect(id.intersection, idc)
    }
  }
  ## Restore labels and units if needed
  nm <- names(M)
  names(lab) <- uvar
  names(un ) <- uvar
  anych <- FALSE
  if(any(c(lab, un) != ''))
    for(i in 1 : ncol(M)) {
      x <- M[[i]]
      ni <- nm[i]
      changed <- FALSE
      if(lab[ni] != '' && ! length(attr(x, 'label'))) {
        label(x) <- lab[ni]
        changed <- TRUE
      }
      if(un[ni] != '' && ! length(attr(x, 'units'))) {
        units(x) <- un[ni]
        changed <- TRUE
      }
      if(changed) M[[i]] <- x
      anych <- anych | changed
    }
  
  nams  <- c(nams, 'Merged')
  va    <- c(va, ncol(M))
  n     <- c(n, nrow(M))
  did   <- if(is.data.table(M)) M[, id, with=FALSE] else M[id]
  idc   <- unique(as.character(interaction(did)))
  nu    <- c(nu, length(unique(idc)))
  nin1  <- c(nin1,  sum(idc %in%  idc1))
  nnin1 <- c(nnin1, sum(idc %nin% idc1))
  info  <- cbind(Vars=va, Obs=n, 'Unique IDs'=nu, 'IDs in #1'=nin1,
                 'IDs not in #1'=nnin1)
  rownames(info) <- nams
  if(verbose) {
    print(info)
    cat('\nNumber of unique IDs in any data frame :', length(id.union), '\n')
    cat(  'Number of unique IDs in all data frames:', length(id.intersection),
        '\n')
    if(anych) cat('\nLabels or units restored\n')
  }
  attr(M, 'info') <- info
  M
}
