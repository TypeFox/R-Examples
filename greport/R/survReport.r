#' Survival Report
#'
#' Generate a Survival Report with Kaplan-Meier Estimates
#'
#' @param formula a formula with survival (\code{Surv}) objects on the left hand side and an optional stratification factor on the right (or \code{1} if none).  The survival object component variables should be labeled; these labels are used for graph annotation.  If any of the \code{Surv} objects are competing risk objects (see \code{\link[survival]{Surv}}), event labels come from the factor levels in the variable that was the second argument to \code{Surv}, and the first factor level must correspond to right-censored observations.
#' @param data data.frame
#' @param subset optional subsetting criteria
#' @param na.action function for handling \code{NA}s while creating a data frame
#' @param ylab character. Passed to \code{\link[rms]{survplot.npsurv}} as the \code{ylab} argument.  Constructed by default.
#' @param what \code{"S"} (the default) to plot survival functions or \code{"1-S"} to plot cumulative incidence functions.  If any of the survival time objects on the left hand side are competing risk objects, the default is \code{"1-S"} and you may not change it.
#' @param conf character. See \code{\link[rms]{survplot.npsurv}}.
#' @param cause character vector or list.  If a vector, every \code{Surv} term on the left hand side of \code{formula} will have cumulative incidence plotted for all causes that appear in \code{cause}.  If a list, the list elements must correspond to the \code{Surv} terms in order, and specify which causes to display from the corresponding \code{Surv} object.  When \code{cause} is a list and one of its elements contains more than one character string, or when \code{cause} is a vector and for one \code{Surv} object it matches multiple causes, \code{survReport} produces more plots than there are \code{Surv} objects.
#' @param panel character string.  Name of panel, which goes into file base names and figure labels for cross-referencing.
#' @param subpanel character string.  If calling \code{dReport} more than once for the same type of chart (categorical or continuous), specify \code{subpanel} to distinguish the multiple calls.  In that case, \code{-subpanel} will be appended to \code{panel} when creating figure labels and cross-references.
#' @param head character string.  Specifies initial text in the figure caption, otherwise a default is used.
#' @param tail optional character string.  Specifies final text in the figure caption, e.g., what might have been put in a footnote in an ordinary text page.  This appears just before any needles.
#' @param h numeric. Height of plots.
#' @param w numeric. Width of plots in inches.
#' @param multi logical.  If \code{TRUE}, multiple figures are produced, otherwise a single figure with a matrix of survival plots is made.
#' @param markevent logical.  Applies only if \code{multi=TRUE}.  Specify \code{FALSE} to not put the event label in the extreme upper left of the plot.
#' @param mfrow numeric 2-vector, used if \code{multi=FALSE}.  If not specified, default plot matrix layout will be figured.
#' @param y.n.risk used if \code{what="1-S"}, to specify \code{y} coordinate for putting numbers at risk, typically below the \code{x}-axis label
#' @param mylim numeric 2-vector.  Used to force expansion of computed y-axis limits.  See \code{survplot}.
#' @param bot number of spaces to reserve at bottom of plot for numbers at risk, if \code{what="1-S"}
#' @param aehaz logical.  Set to \code{FALSE} to not print number of events and hazard rate on plots.
#' @param times numeric vector.  If specified, prints cumulative incidence probabilities at those times on the plots.
#' @param append logical. If \code{TRUE} output will be appended instead of overwritten.
#' @param opts list.  A list specifying arguments to \code{survReport} and \code{startPlot} that override any other arguments.  This is useful when making a long series of \code{survReport} calls with the same options, as the options can be defined up front in a list.
#' @param \dots ignored
#' @export
#' @examples
#' ## See tests directory test.Rnw for a live example
#' \dontrun{
#'   set.seed(1)
#'   n <- 400
#'   dat <- data.frame(t1=runif(n, 2, 5), t2=runif(n, 2, 5),
#'                     e1=rbinom(n, 1, .5), e2=rbinom(n, 1, .5),
#'                     treat=sample(c('a','b'), n, TRUE))
#'   dat <- upData(dat,
#'                 labels=c(t1='Time to operation',
#'                          t2='Time to rehospitalization',
#'                          e1='Operation', e2='Hospitalization',
#'                          treat='Treatment')
#'                 units=c(t1='year', t2='year'))
#'   survReport(Surv(t1, e1) + Surv(t2, e2) ~ treat, data=dat)
#'
#'   dat <- upData(dat, labels=c(t1='Follow-up Time', t2='Time'),
#'                 cause=factor(sample(c('death','MI','censor'), n, TRUE),
#'                              c('censor', 'MI', 'death')))
#'   survReport(Surv(t1, cause) ~ treat, cause='death', data=dat)
#'   survReport(Surv(t1, cause) + Surv(t2, cause) ~ treat,
#'              cause=list(c('death', 'MI'), 'death'), data=dat)
#'   # Two plots for t1, one plot for t2
#' }

survReport <- function(formula, data=NULL, subset=NULL, na.action=na.retain,
                       ylab=NULL, what=c('S', '1-S'),
                       conf=c('diffbands', 'bands', 'bars', 'none'),
                       cause=NULL,
                       panel='surv', subpanel=NULL, head=NULL, tail=NULL,
                       h=3, w=4.5, multi=FALSE, markevent=TRUE, mfrow=NULL, y.n.risk=0,
                       mylim=NULL, bot=2, aehaz=TRUE, times=NULL,
                       append=FALSE, opts=NULL, ...)
{
  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')
  if(length(subpanel) && grepl('[^a-zA-Z-]', subpanel))
    stop('subpanel must contain only A-Z a-z -')

  what <- match.arg(what)
  if(length(cause)) what <- '1-S'
  conf <- match.arg(conf)

  ## Bring arguments from opts as if they were listed outside opts
  if(length(opts) && is.list(opts))
    for(j in 1 : length(opts))
      assign(names(opts)[j], opts[[j]], immediate=TRUE)
  ## Add other startPlot and spar arguments into opts
  w <- c(list(h=h, w=w, multi=multi, mfrow=mfrow, bot=bot), list(...))
  if(! length(opts) || any(names(w) %nin% names(opts))) {
    if(! length(opts)) opts <- list(junk='junk')
    for(x in names(w)) if(x %nin% names(opts)) opts[[x]] <- w[[x]]
    opts$junk <- NULL
  }

  kmlab <- if(what == 'S') 'Kaplan-Meier estimates'
           else 'One minus Kaplan-Meier estimates'

  past <- function(x) {
    l <- length(x)
    y <- if(l < 2) x
    else if(l == 2) paste(x, collapse=' and ')
    else paste(paste(x[1 : (l - 1)], collapse=', '), x[l], sep=', and ')
    upFirst(y, alllower=TRUE)
  }
    
  stamp <- function(w) 
	  text(grconvertX(0, 'ndc', 'user'), grconvertY(1, 'ndc', 'user'),
	       w, adj=c(0, 1), xpd=NA, col='blue')
	  

  form <- Formula(formula)
  Y <- if(length(subset)) model.frame(form, data=data, subset=subset,
                                      na.action=na.action)
  else model.frame(form, data=data, na.action=na.action)
  X <- model.part(form, data=Y, rhs=1)
  Y <- model.part(form, data=Y, lhs=1)
  
  texdir <- getgreportOption('texdir')
  file   <- if(getgreportOption('texwhere') == 'gentex')
    sprintf('%s/%s.tex', texdir, panel) else ''
  if(file != '' && ! append) cat('', file=file)
  lb <- if(length(subpanel)) sprintf('%s-%s', panel, subpanel) else panel

  if(length(cause) && is.list(cause) && length(cause) != length(Y))
    stop('when cause is a list it must have length = number of Surv objects on left side of model')
  

  namx <- labx <- NULL
  if(length(X)) {
    x <- X[[1]]
    namx <- names(X)[1]
    labx  <- upFirst(ifelse(label(x) == '', namx, label(x)), alllower=TRUE)
  }

  Nobs <- nobsY(formula, group=getgreportOption('tx.var'),
                data=data, subset=subset, na.action=na.action)

  ny <- nycause <- ncol(Y)
  ## nycause is the total number of plots counting any separate plots for
  ## multiple causes

  if(length(cause)) {
    Cause <- list()
    for(i in 1 : ny) {
      y <- Y[[i]]
      states <- attr(y, 'states')
      usecause <- ''
      if(length(states)) {
        selectedCauses <- if(is.list(cause)) cause[[i]] else cause
        if(! length(selectedCauses))
          stop(paste('cause not specified for Surv object #', i))
        if(is.list(cause) && any(selectedCauses %nin% states))
          stop(paste('a selected cause is not in the list of states for Surv object #',
                     i, '\nstates:', paste(states, collapse=','),
                     '\ncause:', paste(selectedCauses, collapse=',')))
        usecause <- intersect(states, selectedCauses)
      }
      Cause[[i]] <- usecause
    }
    nycause <- length(unlist(Cause))
  }
  
  if(nycause == 1) multi <- FALSE
  if(! multi) {
    if(! length(opts$mfrow)) opts$mfrow <- mfrowSuggest(nycause)
    if(what == 'S')          opts$bot <- 0
    do.call('startPlot', c(list(file=lb, lattice=FALSE), opts))
  }

  gro <- getgreportOption()
  x.is.tx <- FALSE; ng <- 0
  if(length(X)) {
    x <- X[[1]]
    ng <- if(is.factor(x)) length(levels(x)) else
     length(unique(x[!is.na(x)]))
    if(namx == gro$tx.var) {
      x.is.tx <- TRUE
      col <- gro$tx.linecol
      lwd <- gro$tx.lwd
    }
    else {
      col <- rep(gro$nontx.col, length=ng)
      lwd <- rep(c(1, 3), length=ng)
    }
  } else {
    x <- rep('', nrow(Y))
    col <- 1
    lwd <- 2
  }

  nobs <- rep(0, 1 + x.is.tx * ng)
  evlab <- character(nycause)
  icause <- 0
  for(i in 1 : ny) {
    y <- Y[[i]]

    states <- attr(y, 'states')
    usecause <- ''
    if(length(attr(y, 'states'))) {
      if(! length(cause))
        stop('cause must be specified if any Surv objects on the left side are for competing risks')
      usecause <- Cause[[i]]
    }
    
    no <- c(randomized = nrow(y[! is.na(y)]))
    s <- npsurv(y ~ x)
    if(conf == 'diffbands' && length(s$strata) < 2) conf <- 'bands'
    if(x.is.tx) {
      no        <- c(no, s$n)
      names(no) <- c('randomized', levels(x))
    }

    for(cau in usecause) {
      icause <- icause + 1
      evlab[icause] <- if(cau == '') label(y) else cau
      if(multi) {
        lbi <- paste(lb, icause, sep='-')
        if(what == 'S') opts$bot <- 0
        do.call('startPlot', c(list(file=lbi, lattice=FALSE), opts))
      }
      yl <- ylb <- if(length(ylab)) ylab else upFirst(evlab[icause])
      yl <- if(what == 'S') paste(yl, '-Free Probability', sep='')
      else paste('Cumulative Incidence of', yl)

      cex.ylab <- par('cex.lab') * ifelse(nchar(yl) > 33, .8, 1)
      if(what == 'S')
        survplot(s, 
                 n.risk=TRUE, conf=conf, lwd=lwd,
                 lty=1, col=col, ylab=yl, mylim=mylim,
                 label.curves=list(keys='lines', key.opts=list(bty='n')),
                 levels.only=TRUE, aehaz=aehaz, times=times, 
                 cex.ylab=cex.ylab, ...)
      else
        survplot(s, state=if(length(cause)) cau,
                 fun=function(y) 1 - y,
                 n.risk=TRUE, y.n.risk=y.n.risk, conf=conf, lwd=lwd,
                 lty=1, col=col, ylab=yl, mylim=mylim,
                 label.curves=list(keys='lines', key.opts=list(bty='n')),
                 levels.only=TRUE, aehaz=aehaz, times=times,
                 cex.ylab=cex.ylab, ...)

      capconf <- if(conf == 'diffbands') ', along with half-height of 0.95 confidence limits centered at estimate midpoints. $N$=' else
    ', along with 0.95 confidence bands.  $N$='

      if(multi) {
		if(markevent) stamp(ylb)
        endPlot()
        shortcap <- if(length(head)) head
                    else if(cau == '') paste(kmlab, 'for',
                              upFirst(evlab[icause], alllower=TRUE))
                    else paste('Cumulative incidence of',
                               upFirst(cau, alllower=TRUE),
                               if(length(states) > 2) 'with competing events'
                                else 'with competing event',
                               past(setdiff(states, cau)))
        if(length(labx))
          shortcap <- paste(shortcap, 'stratified by', labx)
        cap <- paste(shortcap, capconf, no[1], '. ', tail, sep='')
        dNeedle(sampleFrac(no, Nobs), name='lttsurv', file=file)
        cap <- sprintf('%s~\\hfill\\lttsurv', cap)
        putFig(panel=panel, name=lbi, caption=shortcap, longcaption=cap)
      }
      if(! multi) for(j in 1:length(nobs)) nobs[j] <- max(nobs[j], no[j])
      names(nobs) <- names(no)
    }
  }
  
  if(! multi) {
    endPlot()
    shortcap <- if(length(head)) head else
        if(! length(cause)) kmlab else
        'Cumulative incidence estimates under competing risks'                                                                                
    shortcap <- paste(shortcap, 'for', past(evlab))
    if(length(labx))
      shortcap <- paste(shortcap, 'stratified by', labx)
    cap <- paste(shortcap, capconf, nobs[1], '. ', tail, sep='')
    dNeedle(sampleFrac(nobs, Nobs), name='lttsurv', file=file)
    cap <- sprintf('%s~\\hfill\\lttsurv', cap)
    putFig(panel=panel, name=lb, caption=shortcap, longcaption=cap)
  }
  invisible()
}
