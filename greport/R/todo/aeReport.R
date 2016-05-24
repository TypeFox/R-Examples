#' Adverse Events Report
#'
#' Generate reports to summarize adverse events.
#'
#' @param data data.frame. Data used for report.
#' @param vars character vector. Variables to include in analysis.
#' @param treat factor vector. Vector of treatment group for each record.
#' @param time character. Name of time variable within dataset.
#' @param times numeric vector. Times used to subset data.
#' @param id character. Name of id variable within dataset.
#' @param plotprop logical. If \sQuote{TRUE} plot proportions of adverse events.
#' @param plotkm logical. If \sQuote{TRUE} plot Kaplan-Meier estimates of
#' cumulative probabilities of adverse events.
#' @param etimedata list. List of time data, containing treatment, event time
#' and event indicator for each variable.
#' @param tables logical. If \sQuote{TRUE} generate summary tables of adverse
#' events by times found in \code{times.tables}.
#' @param times.tables numeric vector. Times used to generate summary tables.
#' @param forceBinary logical. Set to \sQuote{TRUE} if variables are binary.
#' @param ylim numeric vector. Set y-axis limits.
#' @param h numeric. Height of plot. Default is 5in. See \code{\link[Hmisc]{setps}}.
#' @param w numeric. Width of plot. Default is 5in. See \code{\link[Hmisc]{setps}}.
#' @param digits numeric. Number of significant digits to print. Defaults to 3.
#' @export
#' @examples
#' \dontrun{
#'   load(url('http://biostat.mc.vanderbilt.edu/wiki/pub/Main/Rreport/ssafety.rda'))
#'   ae <- c('headache', 'ab.pain', 'nausea', 'dyspepsia', 'diarrhea', 'upper.resp.infect', 'coad')
#'   aeReport(ssafety, ae, 'trx', 'week', id='id', times.tables=c(4,12), forceBinary=TRUE, ylim=c(0,0.15))
#' }

aeReport <- function(data=NULL, vars, treat, time,
                     times=sort(unique(time)),
                     id=NULL, plotprop=FALSE, plotkm=TRUE, etimedata=NULL,
                     tables=TRUE, times.tables=times,
                     forceBinary=FALSE,
                     ylim=c(0,1), h=5, w=5, digits=3)
{
  vars  <- unlist(vars)
  data  <- data[data[[time]] %in% times, c(vars,treat,time,id)]
  Treat <- data[[treat]]
  Time  <- data[[time]]
  nt <- length(levels(Treat))

  ## Useful if variables represent counts and we want present/absent
  if(forceBinary) {
    for(i in 1:length(vars)) {
      x <- data[[i]]
      if(is.numeric(x)) {
        lab <- label(x)
        x <- 1*(x > 0)
        label(x) <- lab
        data[[i]] <- x
      }
    }
  }

  g1 <- function(y) {
    y <- y[!is.na(y)]
    c(Mean=mean(y), var=var(y), n=length(y))
  }

  if(plotprop) {
    startPlot('ae%d', h=h, w=w)
    oldmf <- mfrowSet(length(vars))
    mf <- par('mfrow')
    allbin <- TRUE
    for(y in data[vars]) {
      if(!is.numeric(y)) y <- 1*(y %in% c('Y','Yes','yes','YES'))

      allbin <- allbin & length(unique(y[!is.na(y)])) < 3
      s <- summarize(y, llist(Treat,Time), g1)

      curves <- vector('list',nt)
      names(curves) <- levels(Treat)
      for(w in names(curves))
        curves[[w]] <- list(x=s$Time[s$Treat==w], y=s$y[s$Treat==w])

      yl <- c(ylim[1], max(ylim[2], s$y))
      labcurve(curves, lwd=c(1,2), lty=c(1,1), col.=gray(c(0,.7)),
                xlab=label(Time), ylab=label(y),
                ylim=yl, type='b', method='none', pl=TRUE)

      ## Find out a typical sample size and the pooled variance
      ## to get a typical confidence interval for the
      ## difference in two proportions or two means

      n <- tapply(s$n, s$Treat, median)
      pv <- sum((s$n - 1) * s$var)/sum(s$n - 1)
      clhalfwidth <- 1.96*sqrt(pv*(1/n[1] + 1/n[2]))
      u <- par('usr'); par(xpd=NA)
      lines(rep(u[2]+.03*(u[2]-u[1]),2), c(u[3], u[3]+clhalfwidth),
            col=gray(0.7), lwd=0.5)
    }

    endPlot()
    for(i in 1:ceiling(length(vars) / prod(mf))) {
      cap <- paste('Proportions of adverse events by treatment over time',
                    if(i > 1) ' (continued)'
                    else '', sep='')
      st <- if(allbin) 'proportions'else 'means'

      putFig('ae', paste('ae',i,sep=''), cap,
              if(i == 1) paste(cap,
                  '. Vertical bars to the right of each plot indicate',
                  ' half-widths of typical approximate 0.95 confidence',
                  ' bars for a difference in ',paste(st,'.',sep=''),
                  ' When the distance between two ',st,
                  ' exceeds the length of this bar, differences are',
                  ' significant at approximately the 0.05 level.',
                  ' These confidence bars were computed at the median',
                  ' per-treatment sample size over time',
                  ' and used the pooled variance over treatments and time.',
                  ' \\protect\\treatkey',
                  sep=''),
              append=i > 1)
    }
  }

  if(plotkm) {
    Id <- data[[id]]
    startPlot('ae-km%d', h=h, w=w)
    oldmf <- mfrowSet(min(4, length(vars)))
    mf <- par('mfrow')
    xpd <- par(xpd=NA)
    on.exit(par(xpd=xpd, mfrow=oldmf))
  
    for(v in vars) {
      if(length(etimedata)) {
        tr <- etimedata[[1]]
        etime <- etimedata[[2]][[v]]
        event <- etimedata[[3]][[v]]
        ylab <- v
      } else {
        y <- data[[v]]
        ylab <- label(y, plot=TRUE)
        if(!is.numeric(y)) {
          y <- 1*(y %in% c('Y','Yes','yes','YES'))
        }

        ## For one subject, compute time to first AE of current type,
        ## time at which AE last assessed, and event indicator
        ## Note: using matrices is MUCH faster than using data frames
        ## with by()
        g2 <- function(z) {
          y <- z[,1]
          time <- z[,2]
          s <- is.na(y)

          if(s[1]) return(c(treat=z[1,3], etime=NA, event=NA))

          ## If any NA in y, restrict attention to y's before first NA
          if(any(s)) {
            k <- min(which(s)) - 1
            y <- y[1:k]
            time <- time[1:k]
          }
          etime <- if(any(y>0)) min(time[y>0]) else NA
          ctime <- max(time)
          event <- !is.na(etime)
          if(!event) etime <- ctime
          c(treat=z[1,3], etime=etime, event=event)
        }

        z <- cbind(y, Time, Treat)
        ## Execute g2 separately for all subjects
        r <- mApply(z, Id, g2)
        tr <- factor(r[,1], 1:length(levels(Treat)), levels(Treat))
        etime <- r[,2]
        event <- r[,3]
      }
      km <- survfit.formula(Surv(etime, event) ~ tr)
      omar <- par('mar')
      mar <- omar
      mar[1] <- mar[1]+3
      par(mar=mar)

      lows <- min(km$surv, na.rm=TRUE)
      survplot.survfit(km,
                fun = function(y) 1 - y,
                xlab = label(Time),
                ylab = ylab,
                conf = 'none', label.curves = FALSE,
                ylim = c(0, 1 - lows),
                time.inc = floor(max(times) / 5),
                n.risk = TRUE, y.n.risk = -(1 - lows) * 0.475,
                sep.n.risk = 0.045,
                lwd = c(1, 2), lty = c(1, 1), col = gray(c(0, 0.7)))

      plotKmHalfCL(km, weeks, function(y) 1 - y, offset = max(times) / 50)
      par(mar = omar)
    }
    endPlot()

    for(i in 1:ceiling(length(vars) / prod(mf))) {
      cap <- paste('Kaplan-Meier estimates of cumulative probabilities',
                    ' of adverse events by treatment over time',
                    if(i > 1) ' (continued)'
                    else '',
                    sep='')

      putFig('ae', paste('ae-km', i, sep = ''),
              cap,
              if(i == 1) paste(cap,
                  '. Dotted vertical bars indicate half-widths of',
                  ' approximate 0.95 confidence intervals for',
                  ' differences in probabilities.',
                  ' When the distance between two proportions',
                  ' exceeds the length of the bar, differences are',
                  ' significant at approximately the 0.05 level.',
                  ' \\protect\\treatkey',
                  sep = '') else
              cap,
              append = i > 1)
    }
  }

  if(!tables) return()

  if(!all(times.tables %in% times))
    stop('times.tables must be a subset of times')

  aefile <- file.path(TexDirName(), 'ae.tex')
  cat('In the following tables $N$ is the number of subjects and',
      'numbers after percents are frequencies.',
      '$P$-values are from Pearson $\\chi^2$ tests.\n',
      file = aefile, append = TRUE)

  form <- as.formula(paste(treat, paste(vars, collapse = '+'), sep = '~'))

  for(x in times.tables) {
    s <- summary(form, data = data[Time == x,], method = 'reverse',
                  test = TRUE)
    w <- latex(s, file = aefile, append = TRUE,
                middle.bold = TRUE, title = '',
                caption = paste('Adverse Events at', label(Time), x),
                prtest = 'P', digits = digits, where = 'hbp!',
                insert.bottom = FALSE, ctable = TRUE)
  }

  ## Make a new data frame unioning AEs over all times
  d <- data[c(vars, treat)]
  m <- asNumericMatrix(d)
  m <- mApply(m, Id,
              function(y) {
                nc <- ncol(y)
                v <- apply(y[, -nc, drop=FALSE], 2, any, na.rm = TRUE)
                c(v, tr = y[1,nc])
              })
  
  dm <- dimnames(m)
  dm[[2]][length(dm[[2]])] <- treat
  dimnames(m) <- dm
  v <- matrix2dataFrame(m)
  s <- summary(form, data = v, method = 'reverse', test = TRUE)
  w <- latex(s, file = aefile, append = TRUE,
             middle.bold = TRUE, title = '',
             caption='Adverse Events at Any Time',
             prtest = 'P', digits = digits, where = 'hbp!',
             insert.bottom = FALSE, ctable = TRUE)
}

#' Frequency Report
#'
#' Generate tables with the frequencies found for \code{type}.
#'
#' THIS FUNCTION IS INCOMPLETE
#'
#' @param type character vector. Type for each record.
#' @param panel character. Name for panel.
#' @param treat factor vector. Treatment group for each record.
#' @param longPanel character. Long name for panel.
#' @param typeLabel character. Label for type variable.
#' @param plotprop logical. Set to \sQuote{TRUE} to output
#' a plot of proportions.
#' @param Ntreat numeric vector. Add a column header with \sQuote{N} counts.
#' @param omitZeros logical. If \sQuote{TRUE} remove all zero counts
#' from the frequency table.
#' @param ylim numeric vector. Set y-axis limits. THIS PARAMETER IS NOT USED.
#' @param h numeric. Height of plot. Default is 5in. See \code{\link[Hmisc]{setps}}.
#' @param w numeric. Width of plot. Default is 3in. See \code{\link[Hmisc]{setps}}.
#' @param digits numeric. Number of significant digits to print. Defaults to 3.
#' @param size character. Set LaTeX table font size, see \code{\link[Hmisc]{latex}}.
#' @param longtable logical. Use LaTeX \sQuote{longtable} style.
#' See \code{\link[Hmisc]{latex}}.
#' @param lines.page numeric. Maximum number of lines in the body of a table to be
#' placed on a single page. See \code{\link[Hmisc]{latex}}.
#' @param append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @export

freqReport <- function(type,
                       panel,
                       treat,
                       longPanel = panel,
                       typeLabel = label(type),
                       plotprop = FALSE,
                       Ntreat = NULL,
                       omitZeros = TRUE,
                       ylim = c(0, 1),
                       h = 5,
                       w = 5,
                       digits = 3,
                       size = NULL,
                       longtable = FALSE,
                       lines.page = 40,
                       append = FALSE)
{

  longPanel <- paste(longPanel, 'by', typeLabel)
  nt <- length(levels(treat))
  if(any(is.na(type)))
    type <- ifelse(is.na(type), 'Unspecified', as.character(type))
  
  tab <- table(type)
  if(omitZeros) tab <- tab[tab > 0]

  tab <- as.matrix(tab)
  pan <- paste('O', panel, sep = '')

  freqfile <- file.path(TexDirName(), sprintf("%s.tex", pan))
  w <- latex(tab, file = freqfile,
             title = pan, append = append, rowlabel = typeLabel,
             caption = paste('Frequencies of', longPanel),
             ctable = ! longtable, size = size,
             longtable = longtable, lines.page = lines.page)
  
  tab <- table(type, treat)
  tab <- cbind(tab, Total = rowSums(tab))
  if(omitZeros) tab <- tab[tab[, 'Total'] > 0,]

  pan <- panel
  freqfile <- file.path(TexDirName(), sprintf("%s.tex", pan))
  w   <- latex(tab, file = freqfile,
               title = pan, append = append, rowlabel = typeLabel,
               caption = paste('Frequencies of', longPanel, 'and Treatment'),
               extracolheads = if(length(Ntreat)) {
                 paste('N', c(Ntreat, sum(Ntreat)), sep = '=')
               },
               ctable = ! longtable, size = size,
               longtable = longtable, lines.page = lines.page)

  if(plotprop) {
    startPlot(pan, h = h, w = w)
    ## add plotting code here!
    endPlot()
    cap <- paste('Frequencies and proportions of', longPanel)
    ## NEED TO DEFINE name
    name <- 'unknown figure'
    putFig(pan, name, cap,
           paste(cap,
                 '. Denominators of proportions are assumed to be constants.',
                 sep = ''))
  }
}

#' Report Summarizing All Adverse Events
#'
#' This report uses a data format in which records exist only for events,
#' and these events are classified as \code{major} and \code{minor}
#' events. Typically \code{major} refers to organ system and \code{minor}
#' to preferred term. As only events are signaled in the data, the
#' denominator (number of subjects at risk and having some likelihood of
#' having had adverse events collected had they occured) for each treatment
#' must be specified. Produces an open meeting table not stratified by
#' \code{treatment} (for open report) and one stratified by
#' \code{treatment} (for closed report) in the \code{gentex} directory.
#'
#' Values of \code{minor} are not assumed to be unique across values of \code{major}.
#'
#' @section Side Effects: creates or appends to files
#' \code{"gentex/Ox.tex"} and \code{"gentex/x.tex"} under the current
#'  working directory, where \code{x} is the value of \code{panel}.
#'
#' @param major factor or character vector specifying the major classification of events, usually organ system
#' @param minor factor or character vector specifying the minor
#' classification, usually an individual event given by the preferred
#' term. If \code{minor} is omitted, breakdowns are only by \code{major}.
#' @param treat discrete treatment variable
#' @param id subject ID, used to determine when multiple events arise from the same subject
#' @param denom a vector giving the number of subjects on each treatment
#' (in order of treatment levels), i.e., number of subjects at risk
#' @param panel base name for \code{.tex} files. Default is
#' \code{"ae"}. Table for the open report has an \code{"O"} put at the
#' front, and files are written to directory \code{gentex}.
#' @param caption text string containing the table caption. Other
#' information is automatically appended to this, especially related to
#' \code{minpct} and \code{mindif}. \code{caption} should not end with a period.
#' @param ncaption NEEDDOC
#' @param descending By default, \code{major} categories are sorted in
#' order of descending frequency (ignoring \code{minor}), and
#' \code{minor} categories are sorted in order of descending frequency
#' within each major category. Specify \code{descending='none'} to
#' have categories in order of factor levels or alphabetic for
#' character variables, or to \code{'major'} or \code{'minor'}.
#' @param sortby Specify \code{sortby='\% difference'} to have
#' \code{descending} apply to difference in percent incidence between
#' two treatments, instead of to the pooled incidence over
#' treatments. For the open report, sorting is always alphabetic by
#' category or by descending pooled (over treatments) incidence.
#' @param headerStr NEEDDOC
#' @param minpct a 2-vector specifying the minimum percents for which to
#' output results for \code{major} and \code{minor} categories,
#' respectively. The default is \code{c(0,0)} so that all results are
#' output. Specify \code{minpct=c(10,0)} for example to only print
#' \code{major} categories for which at least a ten percent incidence
#' of events in those categories are present, ignoring \code{treat},
#' and for any qualifying \code{major} category, output all
#' \code{minor} ones. The percents pertain to the per-subject
#' incidences, not counting multiple events per subject.
#' @param mindif a 2-vector specifying the minimum difference in percent
#' per-subject incidence between treatments, corresponding to
#' \code{major} and \code{minor}. The sign of differences is ignored.
#' This only applies to the case of two treatments. Output not meeting
#' the minimum difference in percents is suppressed.
#' @param size see \code{\link[Hmisc]{latex.default}}
#' @param longtable see \code{\link[Hmisc]{latex.default}}
#' @param lines.page see \code{\link[Hmisc]{latex.default}}
#' @param landscape passed to \code{latex.default}
#' @param maxcol the maximum number of columns to allow for a
#' \code{major} or \code{minor} label (note that two extra spaces are
#' added to \code{minor} labels as they are indented). Longer event
#' labels are wrapped to multiple lines. By default, LaTeX will
#' attempt to display the full labels no matter how long.
#' @param append set to \code{TRUE} to have LaTeX code appended to existing files
#' @param ref.label NEEDDOC
#' @param minfreq NEEDDOC
#' @param major.filter NEEDDOC
#' @param minor.filter NEEDDOC
#' @param appendix NEEDDOC
#' @keywords interface category
#' @concept aggregation
#' @seealso \code{\link{aeReport}}
#' @export
#' @examples
#' \dontrun{
#' major <- c('A', 'A', 'A', 'B', 'B', 'C')
#' minor <- c('a1','a1','a2','b1','b2','c1')
#' id <- c(1, 1, 1, 1, 2, 3)
#' treat <- c('drug','placebo','drug','placebo','drug','placebo')
#' aeReport2(major, minor, treat, id, denom=c(placebo=10,drug=10))
#' aeReport2(major, minor, treat, id, denom=c(placebo=10,drug=10),
#'           sortby='% difference')
#' }

aeReport2 <- function(major,
                      minor = rep('', length(major)),
                      treat,
                      id,
                      denom,
                      panel = 'ae',
                      caption = 'Summary of adverse events by system organ class and preferred term',
                      ncaption = '',
                      descending = c('both', 'major', 'minor', 'none'),
                      sortby = c('incidence', '% difference'),
                      headerStr = "AE",
                      minpct = c(0, 0),
                      mindif = c(0, 0),
                      size = NULL,
                      longtable = FALSE,
                      lines.page = 50,
                      landscape = FALSE,
                      maxcol = NULL,
                      append = FALSE,
                      ref.label = NULL,
                      minfreq = 7,
                      major.filter = TRUE,
                      minor.filter = FALSE,
                      appendix = TRUE)
{

  descending <- match.arg(descending)
  sortby <- match.arg(sortby)

  i <- is.na(major) | is.na(minor) | is.na(treat) | is.na(id)

  if(any(i))
    {
      warning(paste(sum(i),'records deleted due to NAs'))
      i <- !i
      major <- major[i]
      minor <- minor[i]
      
      treat <- treat[i]
      id <- id[i]
    }

  major <- as.character(major)
  minor <- as.character(minor)

  doit <- function(file, open.report, treat = rep('', length(major)))
    {      
      treat <- as.factor(treat)
      treats <- levels(treat)
      nt <- length(treats)
      if(any(mindif != 0) && nt > 2)
        stop('mindif only applies when there are two treatments')
    
      acap <- ''
      if(minpct[1] != 0)
        acap <- paste(' Major categories for which the per-subject incidence is',
                      ' $<$ ', minpct[1], '\\% are suppressed.', sep='')
    
    if(minpct[2] != 0)
      acap <- paste(acap,
                    ' Events for which the per-subject',
                    ' incidence is $<$ ', minpct[2], '\\% are suppressed.',
                    sep='')

    if(nt == 2)
      {
        if(mindif[1] != 0)
          acap <- paste(acap,
                        ' Major categories for which the difference in',
                        ' per-subject incidence between treatments is $<$ ',
                        mindif[1], '\\% are suppressed.', sep = '')
        
        if(mindif[2] != 0)
        acap <- paste(acap,
                      ' Events for which the difference in',
                      ' per-subject incidence between treatments is $<$ ',
                      mindif[1], '\\% are suppressed.', sep = '')
      }

      if(appendix && !is.null(minfreq))
        {
          filter.mesg <- " Table includes only those"
          if(major.filter != 0)
            filter.mesg <- paste(filter.mesg, 'major categories')

      if(major.filter != 0 && minor.filter != 0)
        filter.mesg <- paste(filter.mesg, 'and')
      
      if(minor.filter != 0)
        filter.mesg <- paste(filter.mesg, 'events')

          filter.mesg <- paste(filter.mesg,
                               ' where per-subject incidence is $>$ ',
                               minfreq, '.', sep='')
        }

    if(descending != 'none')
      acap <- paste(acap, ' Output is sorted on ',
                    if(descending == 'both') {
                      'major and minor'
                    } else {
                      descending
                    },
                    ' categories based on descending order of ',
                    if(nt==2 && sortby=='% difference') {
                      'absolute differences in incidences between treatments.'
                    } else {
                      'pooled incidences.'
                    }, sep='')
      
    if(nt > 1 && !length(names(denom)))
      {
        warning(paste('assuming that order of frequencies in denom is',
                      paste(treats, collapse = ' ')))
        names(denom) <- treats
      }

    n <- if(nt==1) sum(denom) else denom[treats]

    g <- function(x) length(unique(x))

    w <- function(z)
      {
        z[is.na(z)] <- 0
        z
      }

    sortit <- function(what, sub = 1:length(major))
      {
        x <- if(what=='major') major else minor

        if(descending %nin% c(what,'both'))
          sort(unique(x[sub]))
        else
          {
            if(sortby == '% difference' && nt == 2)
              {
                z <- w(tapply(id[sub], list(x[sub], treat[sub]), g))
                cn <- colnames(z)
                d <- z[, 1, drop = FALSE] /
                  n[cn[1]] -  z[, 2, drop = FALSE] / n[cn[2]]
                rownames(d)[order(-abs(d))]  # sort not work on matrices
              }
            else names(sort(-table(x[sub])))
          }
      }
      
      ae <- sb <- matrix(NA, nrow = 10000, ncol = nt)
      minsubset <- logical(10000)
      
      lab <- 'Any'
      i <- 1
      minsubset[i] <- TRUE
    
      majors <- sortit('major')
      for(maj in majors)
        {
          j <- which(major == maj)

          if(100 * g(id[j]) / sum(denom) < minpct[1]) next

          sbinc <- w(tapply(id[j], treat[j, drop = FALSE], g))

      if(length(sbinc) == 2 && abs(100*diff(sbinc/n)) < mindif[1]) next
        
      isminfreq <- TRUE
      if(appendix && major.filter && !is.null(minfreq) && g(id[j]) < minfreq)
        isminfreq <- FALSE

      lab <- c(lab, '')
        
      lab <- c(lab, maj)
      i <- i + 2
      ae[i,] <- table(treat[j, drop=FALSE])
      sb[i,] <- sbinc

      if(isminfreq) minsubset[c(i-1,i)]   <- TRUE

      m <- sortit('minor', j)
      if(any(m != '', na.rm = TRUE))
        {
          for(mi in m)
            {
              k <- which(major == maj & minor == mi)

              if(100 * g(id[k]) / sum(denom) < minpct[2]) {
                next
              }

              sbinc <- w(tapply(id[k], treat[k,drop = FALSE], g))
          
              if(length(sbinc) == 2 && abs(100 * diff(sbinc / n)) < mindif[2]) next
                 
              if(isminfreq && !(minor.filter && g(id[k]) < minfreq))
                minsubset[i] <- TRUE

              lab <- c(lab, paste('~~', mi, sep=''))
              i <- i + 1
              ae[i,] <- table(treat[k, drop = FALSE])
              sb[i,] <- sbinc
            }
        }
        }
      ae <- ae[1:i,, drop = FALSE]
      sb <- sb[1:i,, drop = FALSE]
      minsubset <- minsubset[1:i]
    
      x <- matrix(NA, nrow = nrow(ae), ncol = nt * 4,
                  dimnames = list(NULL,
                    rep(c(paste("\\#", headerStr, sep = ""),
                          paste("\\#", headerStr, "/N", sep = ""),
                          "\\#Sb",
                          paste("\\%", headerStr, sep = "")),
                        nt)))

      j <- 1
    for(i in 1:nt)
      {
        x[, j:(j+3)] <- cbind(ae[, i], ae[, i] / n[i],
                              sb[,i], 100 * sb[,i] / n[i])
        j <- j + 4
      }
      
    cgroup <- if(nt > 1) paste(treats, ' (N=', n, ')', sep = '')

    extern.ref <- NULL
    if(appendix && (is.null(minfreq) || any(!minsubset)))
      {
        if(is.null(ref.label))
          ref.label <- paste(panel, '.', as.integer(Sys.time()), sep='')

        latex(x, file = AppendixPath(open.report),
              append = TRUE, rowlabel = 'Event',
              cgroup = cgroup,
              n.cgroup = if(nt > 1) rep(4,nt),
              caption = paste('\\label{', ref.label, '}', caption,
                if(nt==1) paste(' (N=',n,')',sep=''),
                ' ', ncaption, '.', acap, sep = ''),
              caption.lot=caption,
              rowlabel.just = if(length(maxcol)) paste('p{',maxcol,'ex}',sep='')
            else 'l',
              where = 'hbp!', cdec = rep(c(0, 3, 0, 1), nt), size = size,
              rowname = lab, longtable = longtable, lines.page = lines.page,
              landscape = landscape)

      if(is.null(minfreq)) return()

        extern.ref <- paste(filter.mesg, ' See Appendix Table \\ref{',
                            ref.label, '} on p.\\ \\pageref{', ref.label,
                            '} for a complete table that includes low-frequency events.',
                            sep='')

        x <- x[minsubset,,drop=FALSE]
        lab <- lab[minsubset]
      }
      
      
      latex(x, file = file.path(TexDirName(), FilenameMask(file, open.report)),
            append = append, rowlabel = 'Event',
            cgroup = cgroup,
            n.cgroup = if(nt > 1){
              rep(4,nt)
            }, caption = paste(caption,
                 if(nt==1) paste(' (N=',n,')',sep=''),
                 ' ', ncaption, '.', acap, extern.ref, sep = ''),
            caption.lot=caption,
            rowlabel.just = if(length(maxcol))
            paste('p{',maxcol,'ex}',sep='') else 'l',
            where = 'hbp!', cdec = rep(c(0, 3, 0, 1), nt), size = size,
            rowname = lab, longtable = longtable, lines.page = lines.page,
            landscape = landscape)

    }

  doit(file=paste(panel, '.tex', sep = ''), open.report=TRUE)
  doit(file=paste(panel, '.tex', sep = ''), open.report=FALSE, treat)
  invisible()
}
