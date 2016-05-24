#' Descriptive Statistics Report
#'
#' Generate graphics and LaTeX with descriptive statistics
#' 
#' \code{dReport} generates multi-panel charts, separately for categorical analysis variables and continuous ones.  The Hmisc \code{summaryP} function and its plot method are used for categorical variables, and \code{bpplotM} is used to make extended box plots for continuous ones unless \code{what='byx'}.   Stratification is by treatment or other variables.  The user must have defined a LaTeX macro \code{\\eboxpopup} (which may be defined to do nothing) with one argument.  This macro is called with argument \code{extended box plot} whenever that phrase appears in the legend, so that a \code{PDF} popup may be generated to show the prototype.  See the example in \code{report.Rnw} in the \code{tests} directory.  Similarly a popup macro \code{\\qintpopup} must be defined, which generates a tooltip for the phrase \code{quantile intervals}.
#'
#' @param formula a formula accepted by the \code{bpplotM} or \code{summaryP} functions.  \code{formula} must have an \code{id(subjectidvariable)} term if there are repeated measures, in order to get correct subject counts as \code{nobs}.
#' @param groups a superpositioning variable, usually treatment, for categorical charts.  For continuous analysis variables, \code{groups} becomes the \code{y}-axis stratification variable.  This is a single character string.
#' @param what \code{"box"} (the default) or \code{"xy"} for continuous analysis variables, or \code{"proportions"} (or shorter) for categorical ones.  Instead, specifying \code{what="byx"} results in an array of quantile intervals for continuous \code{y}, Wilson confidence intervals for proportions when \code{y} is binary, or means and parametric confidence limits when \code{y} is not continuous but is not binary.  If \code{what} is omitted or \code{what="byx"}, actions will be inferred from the most continuous variable listed in \code{formula}.  When \code{fun} is given, different behavior results (see below).
#' @param byx.type set to \code{"quantiles"} to show vertical quantile intervals of \code{y} at each \code{x} for when \code{what="byx"} and the \code{y} variable is continuous numeric, or set \code{byx.type="violin"} (the default) to plot half-violin plots at each \code{x}.
#' @param violinbox set to \code{TRUE} to add violin plots to box plots
#' @param violinbox.opts a list to pass to \code{panel.violin}
#' @param summaryPsort set to \code{TRUE} to sort categories in descending order of frequencies
#' @param stable set to \code{FALSE} to suppress creation of backup supplemental tables for graphics
#' @param fun a function that takes individual response variables (which may be matrices, as in \code{\link[survival]{Surv}} objects) and creates one or more summary statistics that will be computed while the resulting data frame is being collapsed to one row per condition.  Dot charts are drawn when \code{fun} is given.
#' @param data data frame
#' @param subset a subsetting epression for the entire analysis
#' @param na.action a NA handling function for data frames, default is \code{na.retain}
#' @param panel character string.  Name of panel, which goes into file base names and figure labels for cross-referencing
#' @param subpanel If calling \code{dReport} more than once for the same type of chart (by different values of \code{what}), specify \code{subpanel} to distinguish the multiple calls.  In that case, \code{-subpanel} will be appended to \code{panel} when creating figure labels and cross-references.
#' @param head character string.  Specifies initial text in the figure caption, otherwise a default is used
#' @param tail optional character string.  Specifies final text in the figure caption, e.g., what might have been put in a footnote in an ordinary text page.  This appears just before any needles.
#' @param continuous the minimum number of numeric values a variable must have in order to be considered continuous
#' @param h numeric.  Height of plot, in inches
#' @param w numeric.  Width of plot
#' @param outerlabels logical that if \code{TRUE}, pass \code{lattice} graphics through the \code{latticeExtra} package's \code{useOuterStrips}function if there are two conditioning (paneling) variables, to put panel labels in outer margins.
#' @param \dots. Passed to \code{summaryP} or \code{bpplotM}
#' @param append logical.  Set to \code{FALSE} to start a new panel
#' @param sopts list specifying extra arguments to pass to \code{bpplotM}, \code{summaryP}, or \code{summaryS}
#' @param popts list specifying extra arguments to pass to a plot method.  One example is \code{text.at} to specify some number beyond \code{xlim[2]} to leave extra space for numerators and denominators when using \code{summaryP} for categorical analysis variables.  Another common use is for example \code{popts=list(layout=c(columns,rows))} to be used in rendering \code{lattice} plots.  \code{key} and \code{panel} are also frequently used.
#' @param lattice set to \code{TRUE} to use \code{lattice} instead of \code{ggplot2} for proportions.  When this option is in effect, numerators and denominators are shown.
#' @export
#' @examples
#' # See test.Rnw in tests directory

dReport <-
  function(formula, groups=NULL,
           what=c('box', 'proportions', 'xy', 'byx'),
           byx.type=c('violin', 'quantiles'),
           violinbox=TRUE,
           violinbox.opts=list(col=adjustcolor('blue', alpha.f=.25),
             border=FALSE),
           summaryPsort=FALSE,
           stable=TRUE,
           fun=NULL, data=NULL, subset=NULL, na.action=na.retain,
           panel = 'desc', subpanel=NULL, head=NULL, tail=NULL,
           continuous=10, h=5.5, w=5.5, outerlabels=TRUE, append=FALSE,
           sopts=NULL, popts=NULL, lattice=FALSE)
{
  mwhat    <- missing(what)
  what     <- match.arg(what)
  byx.type <- match.arg(byx.type)
  tvar     <- getgreportOption('tx.var')

  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')
  if(length(subpanel) && grepl('[^a-zA-Z-]', subpanel))
    stop('subpanel must contain only A-Z a-z -')

#  rel          <- ggplot2::rel
#  theme        <- ggplot2::theme
#  element_text <- ggplot2::element_text
#  guides       <- ggplot2::guides
#  guide_legend <- ggplot2::guide_legend
  
  center <- 'centerline'
  legend <- NULL

  Nobs <- nobsY(formula, group=tvar,
                data=data, subset=subset, na.action=na.action)
  formula.no.id <- Nobs$formula   ## removes id()
  form <- Formula(formula)
  environment(form) <- new.env(parent = environment(form))
  en <- environment(form)
  assign(envir = en, 'id', function(x) x)

  Y <- if(length(subset)) model.frame(form, data=data, subset=subset,
                                      na.action=na.action)
   else model.frame(form, data=data, na.action=na.action)
  X <- model.part(form, data=Y, rhs=1)
  Y <- model.part(form, data=Y, lhs=1)
  rhs <- terms(form, rhs=1, specials='id')
  sr  <- attr(rhs, 'specials')
  ## specials counts from lhs variables
  wid <- sr$id
  if(length(wid)) wid <- wid - ncol(Y)

  glevels <- if(length(groups)) levels(X[[groups]])
  manygroups <- length(glevels) > 3
  nstrata <- 1
  
  if(mwhat) {
    if(length(fun)) what <- 'xy'
    else {
      y <- Y[[1]]
      if(is.character(y) || is.factor(y) || inherits(y, 'ynbind'))
        what <- 'proportions' else type <- 'box'
    }
  }
  
  labs      <- sapply(X, label)
  if(length(wid)) labs <- labs[- wid]
  stratlabs <- ifelse(labs == '',
                      if(length(wid)) names(X)[-wid] else names(X), labs)
  ylabs     <- sapply(Y, label)
  ylabs     <- ifelse(ylabs == '', names(Y), ylabs)

  past <- function(x) {
    l <- length(x)
    if(l < 2) x
    else if(l == 2) paste(x, collapse=' and ')
    else paste(paste(x[1 : (l - 1)], collapse=', '), x[l], sep=', and ')
  }

  quant <- function(y) {
    probs <- c(0.05, 0.125, 0.25, 0.375)
    probs <- sort(c(probs, 1 - probs))
    y <- y[! is.na(y)]
    if(length(y) < 3) {
      w <- c(median(y), rep(NA, 9), length(y))
      names(w) <- c('Median', format(probs), 'se', 'n')
      return(w)
    }
    w <- hdquantile(y, probs)
    m <- hdquantile(y, 0.5, se=TRUE)
    se <- as.numeric(attr(m, 'se'))
    c(Median=as.numeric(m), w, se=se, n=length(y))
  }

  meanse <- function(y) {
    y <- y[! is.na(y)]
    n <- length(y)
    se <- if(n < 2) NA else sd(y) / sqrt(n)
    if(is.logical(y) || all(y %in% c(0., 1.))) {
      p  <- mean(y)
      ci <- binconf(sum(y), n)[1, ]
      if(p == 0. || p == 1.) {
        ## Don't trust se=0 at extremes; backsolve from Wilson interval
        w  <- diff(ci[c('Lower', 'Upper')])
        se <- 0.5 * w / qnorm(0.975)
      } else se <- sqrt(p * (1. - p) / n)
    }
    else ci <- smean.cl.boot(y, na.rm=FALSE)
    z <- c(ci, se=se, n=length(y))
    z
  }

  propw <- function(y) {
    y <- y[!is.na(y)]
    n <- length(y)
    p <- mean(y)
    ci <- binconf(sum(y), n)[1, ]
    if(p == 0. || p == 1.) {
      ## Don't trust se=0 at extremes; backsolve from Wilson interval
      w  <- diff(ci[c('Lower', 'Upper')])
      se <- 0.5 * w / qnorm(0.975)
    }
    else se <- sqrt(p * (1. - p) / n)
    structure(c(ci, se=se, n=n),
              names=c('Proportion', 'Lower', 'Upper', 'se', 'n'))
  }

  latexit <- function(s, what, byx.type, file) {
    at <- attributes(s)
    xv <- at$xnames
    ## panel function did the work:
    if(what == 'byx.cont' && byx.type == 'violin') {
      g <- function(y) {
        y <- y[! is.na(y)]
        if(length(y) < 3) 
          return(c(n=length(y), Median=median(y), Q1=NA, Q3=NA))
        w <- hdquantile(y, c(0.50, 0.25, 0.75))
        r <- c(length(y), w)
        names(r) <- c('n', 'Median', '0.250', '0.750')
        r
      }
      ## Attempt to find a good number of digits to right of .
      r <- min(tapply(s$y, s$yvar, function(x) max(abs(x), na.rm=TRUE)),
               na.rm=TRUE)
      dig <- if(r == 0) 2
       else max(0, min(5, 3 - round(log10(r))))
      
      s <- with(s, summarize(y, s[c('yvar', xv)],
                             g, type='matrix', keepcolnames=TRUE))
    } else dig <- 2
    sk <- switch(what,
                 byx.cont = c(n='n', Median='Median', Q1='0.250', Q3='0.750'),
                 byx.binary   = c(n='n', Proportion='Proportion'),
                 byx.discrete = c(n='n', Mean='Mean', Lower='Lower',
                   Upper='Upper'))
    cround <- switch(what,
                     byx.cont     = 2:4,
                     byx.binary   = 2,
                     byx.discrete = 2:4)

    s$y <- s$y[, sk, drop=FALSE]
 
    s$y[, cround] <- round(s$y[, cround], dig)
    colnames(s$y) <- names(sk)
    yv <- unique(as.character(s$yvar))
    ny <- length(yv)
    ylab <- character(ny)
    names(ylab) <- yv
    for(v in yv) ylab[v] <-
      labelLatex(label=upFirst(at$ylabels[v]), units=at$yunits[v], hfill=TRUE)
    
    if(length(xv) == 2) {
      r <- reshape(s, timevar=xv[2], direction='wide', idvar=c('yvar', xv[1]))
      class(r) <- 'data.frame'
      lev <- levels(s[[xv[2]]])
      nl <- length(lev)
      yvar <- unique(as.character(r$yvar))
      w <- latex(r[colnames(r) != 'yvar'],
                 table.env=FALSE, file=file, append=TRUE, rowlabel='',
                 landscape=FALSE, size=szg,
                 rowname=rep('', nrow(r)),
                 cgroup=c('', lev),
                 n.cgroup=c(1, rep(ncol(s$y), nl)),
                 rgroup=ylab[yvar],
                 colheads=c(upFirst(xv[1]), rep(names(sk), nl)), center=center)
    }
  else {
    yvar <- unique(as.character(s$yvar))
    w <- latex(s[colnames(s) != 'yvar'],
               table.env=FALSE, file=file, append=TRUE,
               landscape=FALSE,
               rowlabel='', rowname=rep('', nrow(s)),
               rgroup=ylab[yvar], size=szg,
               colheads=c(upFirst(xv[1]), names(sk)), center=center) 
  }
    if(length(xv) == 2) 'full' else 'mini'
  }

  if(what == 'byx') {
    if(length(fun)) stop('may not specify fun= when what="byx"')
    g <- function(y) {
      if(is.logical(y)) 2
      else if(! is.numeric(y)) 0
      else length(unique(y[! is.na(y)]))
    }
    nu <- max(sapply(Y, g))
    what <- if(nu < 3) {
      fun <- propw
      'byx.binary'
    } else if(nu < continuous) {
      fun <- meanse
      'byx.discrete'
    } else {
      if(byx.type == 'quantiles') fun <- quant
      'byx.cont'
    }
  }

  file <- sprintf('%s/%s.tex', getgreportOption('texdir'), panel)
  if(getgreportOption('texwhere') == '') file <- ''
   else if(!append) cat('', file=file)

  cat('%dReport:', deparse(formula), ' what:', what, ' group levels:',
      paste(glevels, collapse=','), '\n',
      file=file, append=TRUE)

  if(what == 'box' && ! length(groups) && ncol(X) == 1)
    manygroups <- length(levels(X[[1]])) > 3

  szg <- if(manygroups) 'smaller[2]' else 'smaller'

  lb <- sprintf('%s-%s', panel, what)
  if(length(subpanel)) lb <- paste(lb, subpanel, sep='-')
  lbn <- gsub('\\.', '', gsub('-', '', lb))
  lttpop <- paste('ltt', lbn, sep='')

  ## Is first x variable on the x-axis of an x-y plot?
  fx <- (what == 'xy' && ! length(fun)) || substring(what, 1, 3) == 'byx'
  a <- if(fx) {
    if(length(ylabs) < 7)
      paste(if(what != 'xy') 'for', past(ylabs), 'vs.\\', stratlabs[1])
     else paste('for', length(ylabs), 'variables vs.\\', stratlabs[1])
  } else paste('for',
               if(length(ylabs) < 7) past(ylabs) else
               paste(length(ylabs), 'variables'))

  al <- upFirst(a, alllower=TRUE)
  al <- latexTranslate(al)
  
  if(!length(head))
    head <-
      switch(what,
       box          = paste('Extended box',
         if(violinbox) 'and violin', 'plots', al),
       proportions  = paste('Proportions', al),
       xy           =  if(length(fun)) 'Statistics' else a,
       byx.binary   = paste('Proportions and confidence limits', al),
       byx.discrete =
             paste('Means and 0.95 bootstrap percentile confidence limits', al),
       byx.cont     = paste('Medians',
         switch(byx.type, quantiles='with quantile intervals',
                violin='with violin (density) plots'),
         al)      )

  sl <- tolower(past(if((what == 'xy' && ! length(fun)) || 
                        what %in% c('byx.binary', 'byx.discrete',
                                    'byx.cont'))
                     stratlabs[-1] else stratlabs))
  cap <- if(!length(sl)) head
  else sprintf('%s stratified by %s', head, sl)

  shortcap <- cap
  tcap <- switch(what,
                 box = paste('Statistics', al),
                 proportions = paste('Proportions', al),
                 xy = if(length(fun)) 'Statistics' else a,
                 byx.binary=paste('Proportions and confidence limits', al),
                 byx.discrete=paste('Means and 0.95 bootstrap CLs', al),
                 byx.cont=paste('Medians', al))
  tcap <- if(length(sl)) sprintf('%s stratified by %s', tcap, sl)
  
  cap <- gsub('Extended box', '\\\\protect\\\\eboxpopup{Extended box}', cap)
  cap <- gsub('quantile intervals', '\\\\protect\\\\qintpopup{quantile intervals}',
              cap)
  
  startPlot(lb, h=h, w=w)
  dl <- list(formula=formula.no.id,
             data=data, subset=subset, na.action=na.action,
             outerlabels=outerlabels)
  key <- popts$key
  if(! length(key) && length(groups)) {
    klines <- list(x=.6, y=-.07, cex=.8,
                   columns=length(glevels), lines=TRUE, points=FALSE)
    key=switch(what,
      box = NULL,
      proportions = list(columns=length(glevels),
        x=.75, y=-.04, cex=.9,
        col=trellis.par.get('superpose.symbol')$col, corner=c(0,1)),
      xy = klines,
      byx.binary =,
      byx.discrete =,
      byx.cont = klines)
  }
  if(length(key)) popts$key <- key

  switch(what,
         box = {
           sopts$violin      <- violinbox
           sopts$violin.opts <- violinbox.opts
           s <- do.call('bpplotM', c(dl, sopts))
           print(s)
         },
         proportions = {
           sopts$sort <- summaryPsort
           s <- do.call('summaryP', c(dl, sopts))
           if(lattice) p <- do.call('plot', c(list(x=s, groups=groups), popts))
           else {
             popts <- if(length(groups) == 1 && groups == tvar)
               c(popts, list(col  =getgreportOption('tx.col'),
                             shape=getgreportOption('tx.pch'),
                             abblen=12))
             else list(col=getgreportOption('nontx.col'), abblen=12)
             popts$addlayer <-
               theme(axis.text.x =
                       element_text(size = rel(0.8), angle=-45,
                                    hjust=0, vjust=1),
                     strip.text.x=element_text(size=rel(0.75), color='blue'),
                     strip.text.y=element_text(size=rel(0.75), color='blue',
                       angle=0),
                     legend.position='bottom')
             p <- do.call('ggplot', c(list(data=s, groups=groups), popts))
             fnvar <- attr(p, 'fnvar')
             if(length(fnvar)) tail <- paste(tail, ' ', fnvar, '.', sep='')
             if(length(groups)) p <- p + guides(color=guide_legend(title=''),
                                                shape=guide_legend(title=''))
           }
           presult <- tryCatch(
             colorFacet(p,
                        col=adjustcolor('blue', alpha.f=0.18)),
             error=function(e) list(fail=TRUE)   )
           if(length(presult$fail) && presult$fail) print(p)
         },
         xy = {
           s <- do.call('summaryS', c(dl, list(fun=fun), sopts))
           p <- do.call('plot', c(list(x=s, groups=groups), popts))
           print(p)
         },
         byx.binary = ,
         byx.discrete =,
         byx.cont = {
           s <- do.call('summaryS', c(dl, list(fun=fun), sopts))
           ylim <- NULL
           if(what %in% c('byx.binary', 'byx.discrete') &&
              all(c('Lower', 'Upper') %in% colnames(s$y))) {
             yvl <- levels(s$yvar)
             ylim <- vector('list', length(yvl))
             names(ylim) <- yvl
             for(yv in levels(s$yvar)) {
               j <- s$yvar == yv
               ylim[[yv]] <- c(min(s$y[j, 'Lower'], na.rm=TRUE),
                               max(s$y[j, 'Upper'], na.rm=TRUE))
             }
           }
           p <- do.call('plot',
            c(list(x=s, groups=groups, ylim=ylim,
                   panel=if(byx.type == 'violin' && what == 'byx.cont')
                         medvPanel else mbarclPanel,
                   paneldoesgroups=TRUE), popts))
           print(p)
         } )

  popname <- paste('poptable', lbn, sep='')
  if(stable) cat(sprintf('\\def\\%s{\\protect\n', popname), file=file, append=TRUE)
  poptab <- NULL

  if(stable && substring(what, 1, 3) == 'byx')
    poptab <- latexit(s, what, byx.type, file=file)
  else if(stable && what == 'proportions') {
    z <- latex(s, groups=groups, size=szg, file=file, append=TRUE,
               landscape=FALSE)   ## may sometimes need landscape=manygroups
    nstrata <- attr(z, 'nstrata')
    poptab <- if(manygroups) 'full' else 'mini'
  }
  else if(what == 'box' || (what == 'xy' && length(fun))) {
    S <- summaryM(formula.no.id, data=data, subset=subset, na.action=na.action,
                  test=FALSE, groups=groups)
    if(stable) {
     z <- latex(S, table.env=FALSE, file=file, append=TRUE, prmsd=TRUE,
                npct='both', exclude1=FALSE, middle.bold=TRUE, center=center,
                round='auto', insert.bottom=FALSE, size=szg,
                landscape=manygroups)
     poptab <- if(length(S$group.freq) > 3) 'full' else 'mini'
     legend <- attr(z, 'legend')
     legend <- if(! length(legend)) ''
     else paste('. ', paste(legend, collapse='\n'), sep='')
     nstrata <- attr(z, 'nstrata')
    }
  }
  if(stable) cat('}\n', file=file, append=TRUE)

  nobs <- Nobs$nobs
  r <- range(nobs)
  nn <- if(r[1] == r[2]) r[1] else paste(r[1], 'to', r[2])
  cap <- sprintf('%s. $N$=%s', cap, nn)
  if(length(tail)) cap <- paste(cap, tail, sep='. ')
  n <- c(randomized=r[2])
  nobsg <- Nobs$nobsg
  if(length(nobsg)) n <- c(n, apply(nobsg, 1, max))
  dNeedle(sampleFrac(n, nobsY=Nobs), name=lttpop, file=file)
  cap <- sprintf('%s~\\hfill\\%s', cap, lttpop)

  endPlot()

  putFig(panel = panel, name = lb, caption = shortcap,
         longcaption = cap,  tcaption=tcap,
         tlongcaption = paste(tcap, legend, sep=''),
         poptable= if(length(poptab)) paste('\\', popname, sep=''),
         popfull = length(poptab) && poptab == 'full',
         outtable = nstrata > 1 || manygroups)
  # hyperref doesn't work with multiple tabulars (runs off page) or landscape
  invisible()
}
