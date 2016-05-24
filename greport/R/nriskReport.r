#' Number at Risk Report
#'
#' Graph number of subjects at risk
#' 
#' \code{nriskReport} generates multi-panel charts, separately for categorical analysis variables.  Each panel depicts the number at risk as a function of follow-up time.  The Hmisc \code{Ecdf} function is used.  Stratification is by treatment or other variables.  It is assumed that this function is only run on randomized subjects.
#' @param formula a formula with time and the left hand side, and with variables on the right side being possible stratification variables.  If no stratification put \code{1} as the right hand side.  Specify unique subject IDs by including a term \code{id()} if subjects have more than one observation.
#' @param groups a character string naming a superpositioning variable.  Must also be included in \code{formula}.
#' @param data data frame
#' @param subset a subsetting epression for the entire analysis
#' @param na.action a NA handling function for data frames, default is \code{na.retain}
#' @param ylab character string if you want to override \code{"Number Followed"}
#' @param panel character string.  Name of panel, which goes into file base names and figure labels for cross-referencing.  The default is \code{'nrisk'}.
#' @param head character string.  Specifies initial text in the figure caption, otherwise a default is used
#' @param tail optional character string.  Specifies final text in the figure caption, e.g., what might have been put in a footnote in an ordinary text page.  This appears just before any needles.
#' @param h numeric.  Height of plot, in inches
#' @param w numeric.  Width of plot
#' @param outerlabels logical that if \code{TRUE}, pass \code{lattice} graphics through the \code{latticeExtra} package's \code{useOuterStrips}function if there are two conditioning (paneling) variables, to put panel labels in outer margins.
#' @param append logical.  Set to \code{FALSE} to start a new panel
#' @param popts list specifying extra arguments to pass to \code{Ecdf}.  A common use is for example \code{popts=list(layout=c(columns,rows))} to be used in rendering \code{lattice} plots.  \code{key} and \code{panel} are also frequently used.
#' @export
#' @examples
#' # See test.Rnw in tests directory

nriskReport <-
  function(formula, groups=NULL,
           data=NULL, subset=NULL, na.action=na.retain,
           ylab='Number Followed', panel = 'nrisk', head=NULL, tail=NULL,
           h=5.5, w=5.5, outerlabels=TRUE, append=FALSE,
           popts=NULL)
{
  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')

  gro  <- getgreportOption()
  tvar <- gro$tx.var
  Nobs <- nobsY(formula, group=tvar,
                data=data, subset=subset, na.action=na.action)
  formula <- Nobs$formula   # removes id()
  
  X <- if(length(subset)) model.frame(formula, data=data, subset=subset,
                                      na.action=na.action)
   else model.frame(formula, data=data, na.action=na.action)
  xnam    <- names(X)
  tx.used <- tvar %in% xnam
  tx      <- if(tx.used) X[[tvar]]
  labs    <- sapply(X, label)
  labs    <- ifelse(labs == '', xnam, labs)
  id <- Nobs$id

  if(length(id) && anyDuplicated(id)) {
    ## Reduce data matrix to one row per subject per stratum with
    ## maximum follow-up time
    X <- data.table(X, .id.=id)
    setnames(X, xnam[1], '.y.')
    by <- if(length(xnam) > 1) paste(xnam[-1], '.id.', sep=',') else '.id.'
    mx <- function(w) as.double(if(any(! is.na(w))) max(w, na.rm=TRUE) else NA)
    X <- X[, list(.maxy.=mx(.y.)), by=by]
    X <- X[, c('.maxy.', xnam[-1]), with=FALSE]
    setnames(X, '.maxy.', xnam[1])
  }

  x1      <- X[[1]]
  xunits  <- units(x1)
  if(xunits == '') xunits <- 'days'
  sl <- if(ncol(X) > 1) upFirst(labs[-1], lower=TRUE)

  file <- sprintf('%s/%s.tex', getgreportOption('texdir'), panel)
  if(getgreportOption('texwhere') == '') file <- ''
   else if(!append) cat('', file=file)
  lb <- gsub('\\.', '', gsub('-', '', panel))
  lbt <- lb
  if(! grepl('nrisk', lb)) {
    lb  <- paste(lb,  'nrisk', sep='-')
    lbt <- paste(lbt, 'nrisk', sep='')
  }
  lttpop <- paste('ltt', lbt, sep='')

  if(! length(head))
    head <- sprintf('Number of subjects followed at least $x$ %s',
                    xunits)
  cap <- if(! length(sl)) head
  else sprintf('%s stratified by %s', head, sl)

  shortcap <- cap

  form <- paste('~', xnam[1])
  cvar <- xnam %nin% c(xnam[1], groups)
  if(any(cvar)) 
    form <- paste(form, '|', paste(xnam[cvar], collapse='*'))
  form <- as.formula(form)
  if(tx.used) {
    col <- gro$tx.linecol
    lwd <- gro$tx.lwd
  } else {
    col <- rep(c(gray(c(0, .7)), 'blue', 'red', 'green'), 10)
    lwd <- rep(c(1, 3), length=10)
  }
  dl <- list(x=form,
             data=X, na.action=na.action,
             what='1-f', col=col, lwd=lwd)
  if(length(subset)) dl$subset <- subset
  if(length(ylab))   dl$ylab   <- ylab
             
  key <- popts$key
  if(! length(key) && length(groups)) {
    glevels <- levels(X[[groups]])
    popts$key <- list(x=.6, y=-.07, cex=.8,
                      columns=length(glevels), lines=TRUE, points=FALSE)
  }

  www <- c(dl, popts)
  p <- if(! length(groups)) do.call('Ecdf', c(dl, popts))
  else {
    a <- sprintf("Ecdf(form, groups=%s, data=X, na.action=na.action, what='1-f', col=col, lwd=lwd", groups)
    if(length(subset)) a <- paste(a, ', subset=subset')
    if(length(ylab))   a <- paste(a, ', ylab=ylab')
    a <- paste(a, ')')
    p <- eval(parse(text=a))
  }
  if(outerlabels && length(dim(p)) == 2) {
#    strip <- function(which.given, which.panel, var.name,
#                      factor.levels, ...) {
#      current.var <- var.name[which.given]
#      levs <- if(current.var == 'time') lev else factor.levels
#      strip.default(which.given, which.panel, var.name, factor.levels=levs, ...#)
#    }
    p <- latticeExtra::useOuterStrips(p) #, strip=strip, strip.left=strip)
  }

  startPlot(lb, h=h, w=w)
  print(p)
  
  if(length(tail)) cap <- paste(cap, tail, sep='. ')
  no <- c(Nobs$nobs, Nobs$nobs, Nobs$nobsg)
  names(no) <- c('enrolled', 'randomized', rownames(Nobs$nobsg))
  dNeedle(sampleFrac(no, Nobs), name=lttpop, file=file)
  cap <- sprintf('%s~\\hfill\\%s', cap, lttpop)
  endPlot()
    
  putFig(panel = panel, name = lb, caption = shortcap,
         longcaption = cap)
  invisible()
}

utils::globalVariables('.y.')
