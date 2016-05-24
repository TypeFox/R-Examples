#' Exclusion Report
#'
#' Generates graphics for sequential exclusion criteria
#' 
#' With input being a series of essentially binary variables with positive indicating that a subject is excluded for a specific reason, orders the reasons so that the first excludes the highest number of subjects, the second excludes the highest number of remaining subjects, and so on.  If a randomization status variable is present, actually randomized (properly or not) subjects are excluded from counts of exclusions.  First draws a single vertical axis graph showing cumulative exclusions, then creates a 2-panel dot chart with the first panel showing that information, along with the marginal frequencies of exclusions and the second showing the number of subjects remaining in the study after the sequential exclusions.  A pop-up table is created showing those quantities plus fractions.  There is an option to not sort by descending exclusion frequencies but instead to use the original variable order.  Assumes that any factor variable exclusions that have only one level and that level indicates a positive finding, that variable has a denominator equal to the overall number of subjects.
#'
#' @param formula a formula with only a right-hand side, possibly containing a term of the form \code{pending(x)} to inform the function of which subjects have incomplete randomization ("randomization pending").  The \code{pending} variable is ignored if a subject has an exclusion marked.  A \code{randomized} variable is an optional \code{logical} vector specifying which subjects are considered to have been randomized.  The presence of this variable causes consistency checking against exclusions.  One or more \code{cond} variables provide binary/logical vectors used to define subsets of subjects for which denominators are used to compute additional fractions of exclusions that are reported in a detailed table.  The arguments of the \code{cond} function are the name of the original variable (assumed to be provided as a regular variable in \code{formula}, a single character string giving the label for the condition, and the vector of essentially binary values that specify the condition.
#' @param data input data frame
#' @param subset subsetting criteria
#' @param na.action function for handling \code{NA}s when creating analysis frame
#' @param ignoreExcl a formula with only a right-hand side, specifying the names of exclusion variable names that are to be ignored when counting exclusions (screen failures)
#' @param ignoreRand a formula with only a right-hand side, specifying the names of exclusion variable names that are to be ignored when counting randomized subjects marked as exclusions
#' @param plotExRemain set to \code{FALSE} to suppress plotting a 2-panel dot plot showing the number of subjects excluded and the fraction of enrolled subjects remaining
#' @param autoother set to \code{TRUE} to add another exclusion \code{Unspecified} that is set to \code{TRUE} for non-pending subjects that have no other exclusions
#' @param sort set to \code{FALSE} to not sort variables by descending exclusion frequency
#' @param whenapp a named character vector (with names equal to names of variables in formula).  For each variable that is only assessed (i.e., is not \code{NA}) under certain conditions, add an element to this vector naming the condition
#' @param erdata a data frame that is subsetted on the combination of \code{id} variables when \code{randomized} is present, to print auxiliary information about randomized subjects who have exclusion criteria
#' @param panel panel string
#' @param subpanel If calling \code{exReport} more than once (e.g., for different values of \code{sort}), specify \code{subpanel} to distinguish the multiple calls.  In that case, \code{-subpanel} will be appended to \code{panel} when creating figure labels and cross-references.
#' @param head character string.  Specifies initial text in the figure caption, otherwise a default is used.
#' @param tail a character string to add to end of automatic caption
#' @param apptail a character string to add to end of automatic caption for appendix table with listing of subject IDs
#' @param h height of 2-panel graph
#' @param w width of 2-panel graph
#' @param hc height of cumulative exclusion 1-panel graph
#' @param wc width of this 1-panel graph
#' @param adjustwidth used to allow wide detailed exclusion table to go into left margin in order to be centered on the physical page.  The default is \code{'-0.75in'}, which works well when using article document class with default page width.  To use the geometry package in LaTeX with margin=.45in specify \code{adjustwidth='+.90in'}.
#' @param append set to \code{TRUE} if adding to an existing sub-report
#' @param popts a list of options to pass to graphing functions
#' @param app set to \code{FALSE} to prevent writing appendix information
#' @author Frank Harrell
#' @export
#' @examples
#' # See test.Rnw in tests directory

exReport <- function(formula, data=NULL, subset=NULL, na.action=na.retain,
                     ignoreExcl=NULL, ignoreRand=NULL, plotExRemain=TRUE,
                     autoother=FALSE, sort=TRUE, whenapp=NULL, erdata=NULL,
                     panel='excl', subpanel=NULL, head=NULL, tail=NULL,
                     apptail=NULL, h=5.5, w=6.5, hc=4.5, wc=5,
                     adjustwidth='-0.75in',
                     append=FALSE, popts=NULL, app=TRUE) {

  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')
  if(length(subpanel) && grepl('[^a-zA-Z-]', subpanel))
    stop('subpanel must contain only A-Z a-z -')

  file <- sprintf('%s/%s.tex', getgreportOption('texdir'), panel)
  if(getgreportOption('texwhere') == '') file <- ''
   else if(!append) cat('', file=file)
  appfile <- sprintf('%s/app.tex', getgreportOption('texdir'))
  subp <- if(length(subpanel)) subpanel else ''

  if(length(ignoreExcl)) ignoreExcl <- all.vars(ignoreExcl)
  if(length(ignoreRand)) ignoreRand <- all.vars(ignoreRand)

#  Nobs <- nobsY(formula, group=getgreportOption('tx.var'),
#                data=data, subset=subset, na.action=na.action)
  environment(formula) <- new.env(parent = environment(formula))
  en <- environment(formula)
  assign(envir = en, 'pending',    function(x) x)
  assign(envir = en, 'randomized', function(x) x)
  assign(envir = en, 'id',         function(x) x)
  gcond <- function(x, label, condition) {
    attr(condition, 'what') <- c(variable = as.character(substitute(x)),
                                 label    = label)
    condition
  }
  assign(envir = en, 'cond', gcond)
  X <- if(length(subset)) model.frame(formula, data=data, subset=subset,
                                      na.action=na.action)
   else model.frame(formula, data=data, na.action=na.action)
  Terms <- terms(formula, specials=c('pending', 'randomized', 'cond', 'id'))
  s <- attr(Terms, 'special')
  sp  <- s$pending
  sr  <- s$randomized
  sc  <- s$cond
  si  <- s$id
  Idnames <- if(length(si)) {
    a <- names(X)[si]
    a <- gsub('id\\(', '', a)
    gsub('\\)', '', a)
  }
         
  ispos <- function(x) {
    w <- if(is.logical(x)) x
    else if(is.numeric(x)) x > 0
    else tolower(as.character(x)) %in%
           c('present', 'yes', 'y', 'true', 'positive')
    w[is.na(x)] <- FALSE
    w
  }

  mis <- function(x) if(is.factor(x) && length(levels(x)) == 1 &&
                        tolower(levels(x)) %in%
                        c('present', 'yes', 'y', 'true', 'positive'))
    rep(FALSE, length(x))
  else is.na(x) | tolower(x) %in% c('unknown','n/a','u','uncertain')

  mblue <- '#0080ff'
  
  N     <- getgreportOption('denom')[c('enrolled', 'randomized')]
  n <- norig <- nrow(X)
  if(n != N['enrolled'])
    warning(sprintf('number of observations (%s) does not equal number enrolled (%s) specified using setgreportOption(denom=)', n, N['enrolled']))

  rnd <- NULL
  if(length(sr)) rnd <- ispos(X[[sr]])

  ig <- if(length(ignoreExcl)) match(ignoreExcl, names(X))
  if(length(ig) && any(is.na(ig)))
    stop('ignoreExcl contains variables not in formula')
  margdenom <- sapply(if(length(c(sp, sr, sc, si, ig)))
                      X[, -c(sp, sr, sc, si, ig)]
                       else X,
                      function(x) sum(! mis(x)))
  Xc <- NULL
  if(length(sc)) {
    Xc <- X[, sc, drop=FALSE]
    colnames(Xc) <- sapply(Xc, function(x) attr(x, 'what')['variable'])
    Xclab        <- sapply(Xc, function(x) attr(x, 'what')['label'   ])
    names(Xclab) <- colnames(Xc)
  }

  Id <- if(length(si)) as.character(interaction(X[si]))

  npend <- 0
  if(length(sp)) {
    pending <- ispos(X[[sp]])
    ## Any observations marked as excluded should have pending ignored
    anyex <- rep(FALSE, n)
    namx  <- names(X)
    for(j in (1 : ncol(X))[- c(sp, sr, sc, si)])
      if(namx[j] %nin% ignoreExcl) anyex <- anyex | ispos(X[[j]])
    pending[anyex] <- FALSE
    ## Same with any observation marked as randomized
    if(length(rnd)) pending[rnd & ! is.na(rnd)] <- FALSE
    npend   <- sum(pending)
    n       <- n - npend
    X       <- X[! pending, - c(sp, sr, sc, si), drop=FALSE]
    if(length(Xc))  Xc  <- Xc [! pending,,  drop=FALSE]
    if(length(rnd)) rnd <- rnd[! pending]
    if(length(Id))  Id  <- Id [! pending]
  }
  else if(length(c(sr, sc, si))) X <- X[, - c(sr, sc, si), drop=FALSE]
  Xname <- names(X)
  k     <- ncol(X)
  Xlab  <- sapply(X, label)
  Xlab  <- ifelse(Xlab == '', Xname, upFirst(Xlab))

  anyre <- rep(FALSE, n)
  if(length(rnd)) {
    exclv <- character(0)
    nexr  <- integer(0)
    Ids   <- Idso <- character(0)
    for(i in 1 : k) {
      if(length(ignoreRand) && Xname[i] %in% ignoreRand) next
      x <- ispos(X[[i]])
      anyre <- anyre | (x & rnd)
      r <- sum(x & rnd, na.rm=TRUE)
      if(r > 0) {
        exclv <- c(exclv, Xlab[i])
        nexr   <- c(nexr,   r)
        if(length(Id)) {
          Ids  <- c(Ids, paste(Id[x & rnd], collapse=', '))
          Idso <- c(Idso, Id[x & rnd])
        }
      }
    }
    if(length(nexr)) {
      nnre  <- sum(anyre, na.rm=TRUE)
      exclv <- c(exclv, 'Total Subjects with Any Exclusion')
      nexr   <- c(nexr, nnre)
      if(length(Ids)) Ids <- c(Ids, '')
      E <- data.frame(Exclusion=exclv, Frequency=nexr)
    }
  }

  if(length(ignoreExcl)) {
    X <- X[names(X) %nin% ignoreExcl]
    k <- ncol(X)
    Xname <- names(X)
    Xlab  <- Xlab[Xname]
  }

  ## If randomization status provided, don't count exclusions on
  ## randomized (rightly or wrongly) subjects, otherwise count all exclusions
  use       <- if(length(rnd)) ! rnd else TRUE
  marg      <- sapply(X, function(x) sum(ispos(x) & use, na.rm=TRUE))
  
  add      <- if(sort) which.max(marg) else 1
  cadd     <- Xname[add]
  exclude  <- ispos(X[[cadd]]) & use   ## new exclusion
  nexclude <- sum(exclude, na.rm=TRUE)
  nexcludec <- if(cadd %in% names(Xc)) sum(exclude & Xc[[cadd]], na.rm=TRUE)
   else nexclude
  X[exclude, ] <- NA  ## only consider subjects not previously excl.
  cond.denom <- n
  cd         <- n - nexclude
  
  if(k > 1) for(i in 2 : k) {
    remain   <- sapply(X, function(x) sum(ispos(x) & use, na.rm=TRUE))
    add      <- if(sort) which.max(remain) else i
    xn       <- Xname[add]
    exclude  <- ispos(X[[xn]]) & use
    nex      <- sum(exclude, na.rm=TRUE)
    nexc <- if(xn %in% names(Xc)) sum(exclude & Xc[[xn]], na.rm=TRUE)
     else nex
    if(nex > 0) {
      cadd         <- c(cadd, xn)
      nexclude     <- c(nexclude,  nex)
      nexcludec    <- c(nexcludec, nexc)
      X[exclude, ] <- NA
      cond.denom   <- c(cond.denom, cd)
      cd           <- cd - sum(exclude, na.rm=TRUE)
    }
  }

  nother <- n - sum(nexclude) - N['randomized']
  othlab <- character(0)
  if(autoother && nother > 0 && 'unspecified' %nin% tolower(cadd)) {
    othlab     <- c(Unspecified = 'Unspecified')
    nexclude   <- c(nexclude,    nother      )
    nexcludec  <- c(nexcludec,   nother      )
    cadd       <- c(cadd,       'Unspecified')
    marg       <- c(marg,       NA           )
    cond.denom <- c(cond.denom, cd           )
    margdenom  <- c(margdenom,  norig        )
  }
  
  if(n - sum(nexclude) != N['randomized'])
    warning(sprintf('number enrolled (%s) minus number excluded (%s) does not equal number randomized (%s) specified to setgreportOption(denom=)',
                    n, sum(nexclude), N['randomized']))
  fracnewTotal  <- nexclude / n
  fracnewRemain <- nexclude / cond.denom
  fracremain    <- 1. - cumsum(nexclude) / n
  marg <- marg[cadd]
  excl <- cadd
  elab <- c(Xlab, othlab)
  u <- rep('', k + length(othlab))
  names(u) <- c(Xname, othlab)
  if(length(whenapp)) u[names(whenapp)] <- paste(whenapp, ', ', sep='')
  b <- ifelse(u == '', paste(' / ', margdenom, sep=''),
                       paste(' (', u, 'n=', margdenom, ')', sep=''))
  elab <- ifelse(margdenom < norig, paste(elab, b, sep=''), elab)
  swr <- function(w, ...) 
    sapply(strwrap(w, ..., simplify=FALSE),
           function(x) paste(x, collapse='\n'))
  elabl <- swr(elab, width=25)
  elabr <-  swr(elab, width=25, exdent=8)
  # When smaller font used, wrap with longer width
  elabl2 <- swr(elab, width=37)
  elabr2 <- swr(elab, width=37, exdent=8)

  names(elab) <- names(elabl) <- names(elabr) <- names(elabl2) <-
    names(elabr2) <- c(Xname, othlab)
  elab <- elab  [cadd]
  ell  <- elabl [cadd]
  elr  <- elabr [cadd]
  ell2 <- elabl2[cadd]
  elr2 <- elabr2[cadd]
  
  ## Make single axis linear graph with cumulative exclusions
  lb <- paste(panel, 'cumex', sep='-')
  if(length(subpanel)) lb <- paste(lb, subpanel, sep='-')

  startPlot(lb, h=hc, w=wc, lattice=FALSE,
            mar=c(1,0,1,0), mgp=c(1.5, .5, 0))
  plot.new()
  m <- sum(nexclude)
  cumex <- cumsum(nexclude)
  r <- c(10 * floor(cumex[1] / 10), 10 * ceiling(m / 10))
  par(usr=c(-1.04, 1.04, rev(r)))
  if(diff(r) < 100)
    axis(2, pos=0, at=seq(r[1], r[2], by= 1), tcl=-.21, labels=FALSE,
         col=gray(.8))
  major <- pretty(r, n=10)
  minor <- seq(r[1], r[2], by=10)
  axis(2, pos=0, at=major, cex.axis=.675)
  axis(2, pos=0, at=minor, tcl=-.25, labels=FALSE)
  points(rep(.0125, length(cumex)), cumex, pch=19, cex=.7, xpd=NA,
         col=mblue)
  side <- 2
  ones <- character(0)
  for(i in 1 : length(cumex)) {
    a <- nexclude[i]
    if(a == 1) {
      ones <- c(ones, ell[i])
      next
    }
    cex <- if(a / m < 0.02) .6 else if(a / m < 0.05) .8 else 1
    el <- if(cex >= .8) ell[i] else ell2[i]
    er <- if(cex >= .8) elr[i] else elr2[i]
    y <- cumex[i]
    u <- if(side == 1) el else er
    v <- paste(if(i == 1) '' else '+', a, '  ', u, sep='')
    ## See if not likely to vertically run into previous entry
    if(i < 3 || (y - cumex[i - 2]) / diff(r) > 0.01) {
      if(side == 1)
        text(-.135, y, v, adj=1, xpd=NA, col=mblue, cex=cex)
      else
        text(.07, y, v, adj=0, xpd=NA, col=mblue, cex=cex)
    } else {
      if(side == 1)
        text(-1, y + 0*.0035 * diff(r), v, adj=0, xpd=NA, col=mblue, cex=cex)
      else
        text(1, y + 0*.0035 * diff(r), v, adj=1, xpd=NA, col=mblue, cex=cex)
    }
    side <- 3 - side
  }
  if(length(ones)) text(.96, r[2] - diff(r)/15,
                        paste(c('One exclusion due to:', ones), collapse='\n'),
                        adj=1, cex=0.65)
  endPlot()
  narnd <- if(length(rnd)) ', for subjects not actually randomized' else ''
  cap <- if(length(head)) head
   else paste('Cumulative number of exclusions ($y$-axis) and number of additional exclusions after exclusions placed higher', narnd, '.', sep='')
  cap <- paste(cap,
               if(sort) 'Exclusions are sorted by descending number of incremental exclusions.'
        else 'Exclusions are in the prespecified order shown in the figure.')
  cap <- paste(cap, N['enrolled'], 'subjects were enrolled,',
               sum(pending),
               'non-excluded subjects are pending randomization, and',
               m, 'subjects were excluded.')
  if(length(rnd)) cap <- paste(cap, sum(rnd, na.rm=TRUE),
                               'subjects were randomized.')
  
  wrn1 <- wrn2 <- character(0)
  if(norig != N['enrolled'])
    wrn1 <- sprintf('\\textbf{Note}: Number of observations (%s) does not equal number officially enrolled (%s).',
                    norig, N['enrolled'])
  if(n - m != N['randomized'])
    wrn2 <- sprintf('\\textbf{Note}: Number of enrolled (%s) minus number excluded (%s) does not match official number randomized (%s).',
                    n, m, N['randomized'])

  cap <- paste(cap, tail, wrn1, wrn2)

  putFig(panel=panel, name=lb, caption='Cumulative exclusions',
         longcaption=cap)

  rf <- function(x) format(round(x, 3))
  f  <- function(x) ifelse(is.na(x), '', format(x))
  tabl <- data.frame(elab      = c(latexTranslate(elab), '\\textbf{Total}'),
                     nexclude  = c(nexclude, m),
                     marg      = c(marg, NA),
                     frac      = rf(c(nexclude / n, m / n)),
                     frace     = rf(c(nexclude / m, 1)),
                     fracremain= rf(c(fracremain, (n - m) / n)),
                     row.names = 1 : (length(elab) + 1),
                     stringsAsFactors=FALSE)

  ct <- function(...) cat(..., sep='', file=file, append=TRUE)
  
  ct('\\begin{table}[htbp]\\small\n',
     '\\caption[Exclusions]{Exclusions', narnd, '.  \\texttt{Incremental Exclusions} are those in addition to exclusions in earlier rows.  \\texttt{Marginal Exclusions} are numbers of subjects excluded for the indicated reason whether or not she was excluded for other reasons.  The three \\texttt{Fractions} are based on incremental exclusions.', if(length(tail))' ', tail,
     '\\label{tab:exclstats', subp, '}}\n',
     '\\begin{center}\\begin{adjustwidth}{', adjustwidth, '}{',
     adjustwidth, '}\n',
     '\\begin{tabular}{lrrrrr}\\hline\\hline\n',
     '\\multicolumn{1}{c}{Exclusions}&\\multicolumn{1}{c}{Incremental}&\\multicolumn{1}{c}{Marginal}&\\multicolumn{1}{c}{Fraction of}&\\multicolumn{1}{c}{Fraction of}&\\multicolumn{1}{c}{Fraction}\\tabularnewline',
     '&\\multicolumn{1}{c}{Exclusions}&\\multicolumn{1}{c}{Exclusions}&\\multicolumn{1}{c}{Enrolled}&\\multicolumn{1}{c}{Exclusions}&\\multicolumn{1}{c}{Remaining}\\tabularnewline\\hline\n')
  for(i in 1 : nrow(tabl)) {
    with(tabl, ct(as.character(elab[i]),  '&',
                  nexclude[i], '&',
                  f(marg)[i],  '&',
                  frac[i],     '&',
                  frace[i],    '&',
                  fracremain[i], '\\tabularnewline'))
    cn <- cadd[i]
    if(cn %in% names(Xc)) {
      ct('\\multicolumn{6}{l}{~~~~$\\frac{', nexcludec[i], '}{',
         sx <- sum(Xc[, cn], na.rm=TRUE), '}$ =',
         rf(nexcludec[i] / sx), ' of ', latexTranslate(Xclab[cn]),
         '}\\tabularnewline\n')
      ct('&&&&&\\tabularnewline\n')
    }
    if(i == (nrow(tabl) - 1)) ct('\\hline')
    ct('\n')
  }
  ct('\\hline\\end{tabular}\\end{adjustwidth}\\end{center}\\end{table}\n\n')

  
  ## Two-panel dot chart
  if(plotExRemain) {
    if(npend > 0) {
      elab       <- c('Pending Randomization', elab)
      nexclude   <- c(npend, nexclude)
      marg       <- c(npend,     marg)
      fracremain <- c(NA, fracremain)
    }
    excl  <- factor(elab,          levels=rev(elab))
    excl2 <- factor(c(elab, elab), levels=rev(elab))
    x <- c(nexclude, marg)
    j <- length(elab)  # was k + (npend > 0)
    hh <- c(rep('Incremental Exclusions', j),
            rep('Single Exclusions',      j))
    fracremain <- c(fracremain, rep(NA,   j))
    
    panel.ex <- 
      function (x, y, groups, ..., pch, col) {
        pn <- panel.number()
        up <- max(x, na.rm=TRUE)
        if (pn == 1) {
          ww <- 10 * floor(min(x, na.rm=TRUE) / 10)
          by <- if(up - ww < 100) c(5, 10) else c(10, 50)
          panel.abline(v=seq(ww, up, by=by[1]), lwd=.4, col=gray(.75))
          panel.abline(v=seq(ww, up, by=by[2]), lwd=.6, col=gray(.56))
          panel.abline(h = y, lwd = .4, col = gray(.7))
          panel.superpose(x, y, groups = groups, pch = pch, col = col, ...)
        } else {
          u <- .01 * floor(min(x, na.rm=TRUE) / 0.01)
          ww <- .05 * floor(min(x, na.rm=TRUE) / 0.05)
          up <- max(x, na.rm=TRUE)
          panel.abline(v = seq(u, up, by=0.01), lwd=.4, col=gray(.75))
          panel.abline(v = seq(ww, up, by=0.05), lwd=.6, col=gray(.56))
          panel.dotplot  (x, y, pch=20, col=col[1], ...)
        }
      }
    
    col <- c('black', mblue)
    r <-
      dotplot(excl2 ~ x + fracremain, panel=panel.ex,
              groups=hh,
              pch=c(20, 18), col=col,
              outer = TRUE,
              scales = list(x = list(relation = "free")),
              xlab = NULL, between = list(x = 1),
              key = list(
                points = list(col = col, pch = c(20,18)),
                text = list(c('Sequential (Incremental) Exclusions',
                  'Individual (Marginal) Exclusions'), col = col,
                  cex = 0.7),
                columns = 2, between = 0.5, space = "bottom"))
    r$condlevels[[1]] <- c("Number Excluded", "Fraction Remaining")
    
    lb <- paste(panel, 'nexfrac', sep='-')
    if(length(subpanel)) lb <- paste(lb, subpanel, sep='-')
    
    startPlot(lb, h=h, w=w)
    print(r)
    endPlot()
    
    cap <- if(length(head)) head
     else sprintf('Left panel: Incremental (sequential) and marginal (each exclusion treated separately) exclusions.  Right panel: Fraction of subjects remaining after incremental exclusions.  The denominator of the fraction is the number of subjects not pending randomization (%s).',
                  n)
    cap <- paste(cap,
     if(sort) 'Exclusions are sorted by descending number of incremental exclusions.'
     else 'Exclusions are in the prespecified order shown in the figure.')
    cap <- paste(cap, N['enrolled'], 'subjects were enrolled,',
                 sum(pending),
                 'non-excluded subjects are pending randomization, and',
                 m, 'subjects were excluded.')
    if(length(rnd)) cap <- paste(cap, sum(rnd, na.rm=TRUE),
                                 'subjects were randomized.')
    
    cap <- paste(cap, tail, wrn1, wrn2)
    
    putFig(panel=panel, name=lb,
           caption='Incremental exclusions and fraction of remaining subjects',
           longcaption=cap)
  }
    
  ## If needed, display subjects marked as randomized who are marked as
  ## meeting exclusion criteria
  if(length(rnd) && length(nexr)) {
    cat('\\clearpage\n', file=file, append=TRUE)
    cap   <- 'Frequency of exclusions for subjects marked as randomized'
    scap  <- 'Exclusions in randomized subjects'
    z     <- latex(E, file=file, append=TRUE,
                   label=sprintf('tab:exclrand%s', subp),
                   hyperref=if(app && length(Ids))
                   sprintf('tab:randsubjexcl%s', subp),
                   rowname=NULL, col.just=c('l', 'r'),
                   caption=cap, caption.lot=scap, where='htbp')
    if(app && length(Ids)) {
      if(length(apptail)) apptail <- paste('.', apptail)
      cat('\\begin{table}[htbp]\\caption[Subject IDs for randomized subjects with exclusions]{Subject IDs for randomized subjects with exclusions', apptail, '}\\label{tab:randsubjexcl', subp, '}\n\\medskip%\n',
          sep='', file=appfile, append=TRUE)
      cat(sprintf('\\hyperref[tab:exclrand%s]{$\\leftarrow$}\n\n', subp),
          file=appfile, append=TRUE)
      le <- length(nexr) - 1
      for(i in 1 : le) {
        cat('\\textbf{', latexTranslate(as.character(E$Exclusion[i])),
            '}:\\\\\n',
            file=appfile, append=TRUE, sep='')
        cat('\\parbox{5in}{', Ids[i], '}',
            if(i < le) '\\\\\n', sep='',
            file=appfile, append=TRUE)
      }
      if(length(erdata)) {
        erd <- erdata[as.character(interaction(erdata[Idnames]))
                      %in% Idso, ]
#        colnames(erd) <- latexTranslate(colnames(erd))
        z <- function(x) ifelse(is.na(x), '', as.character(x))
#        for(j in 1 : ncol(erd))
#          erd[[j]] <- latexTranslate(z(erd[[j]]))
        z <- latex(erd, file=appfile, append=TRUE, rowname=NULL,
                   table.env=FALSE, na.blank=TRUE)
      }
      cat('\\end{table}\n\n', file=appfile, append=TRUE)
    }
  }
  
  invisible()
}
