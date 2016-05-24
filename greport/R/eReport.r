#' Event Report
#'
#' Generates graphics for binary event proportions
#' 
#' Generates dot charts showing proportions on left and risk difference with confidence intervals on the right, if there is only one level of event categorization.  Input data must contain one record per event, with this record containing the event name.  If there is more than one event of a given type per subject, unique subject ID must be provided.  Denominators come from \code{greport} options and it is assumed that only randomized subjects have records.  Some of the graphics functions are modifications of those found in the HH package.  The data are expected to have one record per event, and non-events are inferred from \code{setgreportOption('denom')}.  It is also assumed that only randomized subjects are included in the dataset.  
#'
#' @param formula a formula with one or two left hand variables (the first representing major categorization and the second minor), and 1-2 right hand variables.  One of these may be enclosed in \code{id()} to indicate the presence of a unique subject ID, and the other is treatment.
#' @param data input data frame
#' @param subset subsetting criteria
#' @param na.action function for handling \code{NA}s when creating analysis frame
#' @param minincidence a number between 0 and 1 specifying the minimum incidence in any stratum that must hold before an event is included in the summary
#' @param conf.int confidence level for difference in proportions
#' @param etype a character string describing the nature of the events, for example \code{"adverse events"}, \code{"serious adverse events"}.  Used in figure captions.
#' @param panel panel string
#' @param subpanel a subpanel designation to add to \code{panel}
#' @param head character string.  Specifies initial text in the figure caption, otherwise a default is used.
#' @param tail a character string to add to end of automatic caption
#' @param h height of graph
#' @param w width of graph
#' @param append set to \code{TRUE} if adding to an existing sub-report
#' @param popts a list of options to pass to graphing functions
#' @author Frank Harrell
#' @export
#' @examples
#' # See test.Rnw in tests directory

eReport <- function(formula, data=NULL, subset=NULL, na.action=na.retain,
                    minincidence=0, conf.int=0.95,
                    etype='adverse events',
                    panel='events', subpanel=NULL, head=NULL, tail=NULL,
                    h=6, w=7, append=FALSE, popts=NULL) {

  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')
  if(length(subpanel) && grepl('[^a-zA-Z-]', subpanel))
    stop('subpanel must contain only A-Z a-z -')

  Nobs <- nobsY(formula, group=getgreportOption('tx.var'),
                data=data, subset=subset, na.action=na.action)
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

  nY <- ncol(Y)
  major <- NULL
  if(nY > 1) major <- Y[[1]]
  event <- Y[[nY]]
  id <- 1 : length(event)
  nX <- ncol(X)
  gname <- glabel <- ''
  if(nX > 1 + (length(wid) > 0))
    stop('cannot have more than one right hand variable other than id variable')
  if(length(wid)) {
    id    <- X[[wid]]
    j <- setdiff(1 : nX, wid)
  } else if(nX == 1) j <- 1
    else j <- 0
  if(j == 0) {
    group <- factor(rep('', length(event)))
    gname <- glabel <- ''
  }
  else {
    group <- X[[j]]
    gname <- names(X)[j]
    glabel <- label(group, default=gname)
  }

  event <- as.factor(event)
  levels(event) <- upFirst(levels(event))
#  event  <- as.character(event)
  uevent <- levels(event)
  nue    <- length(uevent)
  N <- getgreportOption('denom')
  n <- N[setdiff(names(N), c('enrolled', 'randomized'))]
  groups <- names(n)
  if(length(groups) != 2) stop('currently only implemented for 2 treatments')
  group <- as.character(group)

  if(! length(major)) {
    ## The following functions are in the HH package
    panel.ae.dotplot <-
      function (x, y, groups, ..., col.AB, pch.AB, lower, upper) {
        pn <- panel.number()
        if (pn == 1)
          panel.ae.leftplot(x, y, groups = groups,
                            col = col.AB, pch = pch.AB, ...)
        if (pn == 2)
          panel.ae.rightplot(x, y, ..., lwd = 6, pch = 16,
                             lower = lower, upper = upper)
      }
    
    panel.ae.leftplot <- function(x, y, groups, col.AB, ...) {
      panel.abline(h = y, lty = 2, lwd = .4, col = gray(.7))
      panel.superpose(x, y, groups = groups,
                      col = col.AB, ...)
    }
    
    panel.ae.rightplot <-
      function(x, y, ..., lwd = 6, lower, upper) {
        panel.abline(v = 0, lty = 1, lwd = .6, col = gray(.7))
        panel.abline(h = y, lty = 2, lwd = .4, col = gray(.7))
        panel.segments(lower, y, upper, y, lwd = 2)
        panel.xyplot(x, y, ..., col = 1, cex = 0.7)
        panel.points(lower, y, pch = 3, col = 1, cex = 0.4)
        panel.points(upper, y, pch = 3, col = 1, cex = 0.4)
      }
    
    ## Modification of HH's ae.dotplot that uses risk difference and
    ## uses proportions instead of percents

    zcrit <- qnorm( (1 + conf.int) / 2)
    f <- function(i) {
      idi <- id[i]
      grp <- group[i]
      n1 <- length(unique(idi[grp == groups[1]]))
      n2 <- length(unique(idi[grp == groups[2]]))
      p1 <- n1 / n[1]
      p2 <- n2 / n[2]
      se <- sqrt(p1 * (1 - p1) / n[1] + p2 * (1 - p2) / n[2])
      diff  <- p1 - p2
      lower <- diff - zcrit * se
      upper <- diff + zcrit * se
      r <- c(p1, p2, diff, lower, upper)
      names(r) <- c('p1', 'p2', 'diff', 'lower', 'upper')
      r
    }
    z  <- tapply(1 : length(event), event, f)
    z  <- t(sapply(z, function(x) x))
    if(minincidence > 0) {
      small <- pmax(z[, 'p1'], z[, 'p2'], na.rm=TRUE) < minincidence
      z <- z[! small, ]
    }
    zr <- round(z, 3)
    z  <- data.frame(event=row.names(z ),          z )
    zr <- data.frame(event=upFirst(row.names(zr)), zr)

    file <- sprintf('%s/%s.tex', getgreportOption('texdir'), panel)
    if(getgreportOption('texwhere') == '') file <- ''
     else if(! append) cat('', file=file)
    lb <- if(length(subpanel)) sprintf('%s-%s', panel, subpanel) else panel
    lbn <- gsub('\\.', '', gsub('-', '', lb))

    popname <- sprintf('\\pop%s', lbn)
    names(zr) <- c('Event', groups, 'Difference', 'Lower', 'Upper')
    pop <- latexTabular(zr, align='lrrrrr', hline=1)
    cat('\\def', popname, '{', pop, '}%\n', sep='', file=file, append=TRUE)

    ## Duplicate diff and CI so that lattice can have separate records
    ## for 2 treatments
    z <- rbind(cbind(.group.=groups[1], proportion=z$p1, z),
               cbind(.group.=groups[2], proportion=z$p2, z))

    txcol <- getgreportOption('tx.col')
    txpch <- getgreportOption('tx.pch')
    z$event <- with(z, reorder(event, diff))
    
    r <- dotplot(event ~ proportion + diff,
                 groups = .group.,
                 data = z, outer = TRUE,
                 lower = z$lower, 
                 upper = z$upper,
                 panel = panel.ae.dotplot, 
                 scales = list(x = list(relation = "free",
                                 limits = list(range(z$proportion), 
                                   range(z$lower, z$upper,
                                         na.rm=TRUE))), 
                   y = list(cex = 0.6)),
                 A.name = groups[1],  B.name = groups[2],
                 col.AB = txcol,
                 pch.AB = txpch,
                 cex.AB.points = NULL, 
                 cex.AB.y.scale = 0.6,
                 ## main = list(main.title, cex = main.cex),
                 xlab = NULL, between = list(x = 1), 
                 key = list(y = -0.2, x = 0.15,
                   points = list(col = txcol, pch = txpch),
                   text = list(groups, col = txcol, cex = 0.9),
                   columns = 2, between = 0.5, space = "bottom"))
        r$condlevels[[1]] <- c("Proportion",
                               paste("Risk Difference with",
                                     conf.int, "CI"))
  }
  else {
    stop('major event grouping not yet implemented')
  }
  
  lttpop <- paste('ltt', lbn, sep='')
  if(! length(head))
    head <- paste('Proportion of', etype,
                  'and risk differences by', upFirst(glabel, lower=TRUE),
                  'sorted by risk difference')
  if(minincidence > 0 && any(small))
    head <- paste(head, '. ', sum(small), ' events with less than ',
                  minincidence,
                  ' incidence in at least one group are not shown.',
                  sep='')
  shortcap <- paste('Proportion of', etype, 'and risk differences by',
                    upFirst(glabel, lower=TRUE))
  startPlot(lb, h=h, w=w)
  print(r)
  endPlot()
  N[1] <- N[2]   # assume only analyze randomized subjects
  dNeedle(sampleFrac(N, nobsY=Nobs), name=lttpop, file=file)
  head <- paste(head, '~\\hfill\\', lttpop, sep='')
  putFig(panel=panel, name=lb, caption=shortcap, longcaption=head,
         poptable=popname, popfull=TRUE)
  invisible()
}
