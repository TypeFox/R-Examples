#' Accrual Report
#'
#' Generate graphics and LaTeX to analyze subject accrual
#'
#' Typically the left-hand-side variables of the formula, in order, are date of enrollment and date of randomization, with subjects enrolled but not randomized having missing date of randomization.  Given such date variables, this function generates cumulative frequencies optionally with target enrollment/randomization numbers and with time-zooming.  Makes a variety of dot charts by right-hand-side variables:  number of subjects, number of sites, number of subjects per site, fraction of enrolled subjects randomized, number per month, number per site-month.
#'
#' @param formula formula object, with time variables on the left (separated by +) and grouping variables on the right.  Enrollment date, randomization date, region, country, and site when present must have the variables in parenthesis preceeded by the key words \code{enrollment, randomize, region, country, site}.
#' @param data data frame.
#' @param subset a subsetting epression for the entire analysis.
#' @param na.action a NA handling function for data frames, default is \code{na.retain}.
#' @param dateRange \code{Date} or character 2-vector formatted as \code{yyyy-mm-dd}.  Provides the range on the \code{x}-axis (before any zooming).
#' @param zoom \code{Date} or character 2-vector for an option zoomed-in look at accrual.
#' @param targetN integer vector with target sample sizes over time, same length as \code{targetDate}
#' @param targetDate \code{Date} or character vector corresponding to \code{targetN}
#' @param closeDate \code{Date} or characterstring.  Used for randomizations per month and per site-month - contains the dataset closing date to be able to compute the number of dates that a group (country, site, etc.) has been online since randomizating its first subject.
#' @param enrollmax numeric specifying the upper y-axis limit for cumulative enrollment when not zoomed
#' @param studynos logical.  Set to \code{FALSE} to suppress summary study numbers table.
#' @param minrand integer.  Minimum number of randomized subjects a country must have before a box plot of time to randomization is included.
#' @param panel character string.  Name of panel, which goes into file base names and figure labels for cross-referencing.
#' @param h numeric.  Height of ordinary plots, in inches.
#' @param w numeric.  Width of ordinary plots.
#' @param hb numeric.  Height of extended box plots.
#' @param wb numeric.  Weight of extended box plots.
#' @param hdot numeric.  Height of dot charts in inches.
#' @export
#' @examples
#' \dontrun{
#' # See test.Rnw in tests directory
#' }

accrualReport <-
  function(formula, data=NULL, subset=NULL, na.action=na.retain,
           dateRange=NULL, zoom=NULL, targetN=NULL, targetDate=NULL,
           closeDate=NULL, enrollmax=NULL, studynos=TRUE,
           minrand=10, panel = 'accrual',
           h=2.5, w=3.75, hb=5, wb=5, hdot=3.5)
{
  formula <- Formula(formula)
  
  if(grepl('[^a-zA-Z-]', panel))
    stop('panel must contain only A-Z a-z -')
  
  environment(formula) <- new.env(parent = environment(formula))
  en <- environment(formula)
  f <- function(x) x
  assign(envir = en, "enroll",    f)
  assign(envir = en, "randomize", f)
  assign(envir = en, "region",    f)
  assign(envir = en, "country",   f)
  assign(envir = en, "site",      f)

  file <- sprintf('%s/%s.tex', getgreportOption('texdir'), panel)
  if(getgreportOption('texwhere') == '') file <- ''
   else cat('', file=file)
  
  ltt <- function(used, name='ltt')
      dNeedle(sampleFrac(used),
              name=name, file=file, append=TRUE)

  lhs  <- terms(formula, lhs=1, specials=c('enroll', 'randomize'))
  sl   <- attr(lhs, 'specials')
  rhs  <- terms(formula, rhs=1, specials=c('region', 'country', 'site'))
  sr   <- attr(rhs, 'specials')

  Y <- if(length(subset))
    model.frame(formula, data=data, subset=subset, na.action=na.keep)
   else model.frame(formula, data=data, na.action=na.keep)
  X    <- model.part(formula, data=Y, rhs=1)
  Y    <- model.part(formula, data=Y, lhs=1)
  nY   <- NCOL(Y)
  nX   <- NCOL(X)
  namY <- all.vars(lhs)
  namX <- all.vars(rhs)
  enroll    <- sl$enroll
  randomize <- sl$randomize

  z <- function(x, nY) if(length(x)) x - nY else NULL
  ## specials counts from lhs variables
  region    <- z(sr$region,  nY)
  country   <- z(sr$country, nY)
  site      <- z(sr$site,    nY)

  penroll    <- length(enroll)    > 0
  prandomize <- length(randomize) > 0
  pregion    <- length(region)    > 0
  pcountry   <- length(country)   > 0
  psite      <- length(site)      > 0
  pclose     <- length(closeDate) > 0
  
  cr <- pcountry || pregion

  dr <- dateRange
  if(!length(dr))
    dr <- range(pretty(do.call('range', c(as.list(Y), na.rm=TRUE))))
  else dr <- as.Date(dr)
  if(length(targetN) && ! length(targetDate))
    stop('must provide targetDate if using targetN')
  if(length(targetDate)) targetDate <- as.Date(targetDate)
  if(pclose)  closeDate  <- as.Date(closeDate)

  ylabs <- namY
  for(i in 1 : nY) {
    if(penroll    && enroll == i)    ylabs[i] <- 'enrolled'
    if(prandomize && randomize == i) ylabs[i] <- 'randomized'
  }
  xlabs <- namX
  for(i in 1 : nX) {
    if(pregion  && region == i)  xlabs[i] <- 'region'
    if(pcountry && country == i) xlabs[i] <- 'country'
    if(psite    && site == i)    xlabs[i] <- 'site'
  }

  z <- k <- character(0)
  g <- function(x, digits) as.character(round(x, digits))
  if(pcountry) {
    z <- g(length(unique(X[[country]])), 0)
    k <- 'Countries'
  }
  if(psite) {
    Site <- as.character(X[[site]])
    nsites <- length(unique(Site))
    z <- c(z, g(nsites, 0))
    k <- c(k, 'Sites')
  }
  if(penroll) {
    z <- c(z, sum(! is.na(Y[[enroll]])))
    k <- c(k, 'Subjects enrolled')
  }

  if(psite && prandomize) {
    rdate <- Y[[randomize]]
    nrand <- sum(! is.na(rdate))
    persite <- nrand / nsites
    z <- c(z, c(nrand, g(persite, 1)))
    k <- c(k, c('Subjects randomized', 'Subjects per site'))
    ## maxs = for each site the # months since that site first randomized
    ##        a subject (NA if none randomized)
    ## site months is sum of maxs
    ## avg. months since first randomized = mean maxs excluding NAs
    ## rand per site per month = # rand / site months
    ## Note: # rand / # sites / avg. months != rand per site per month
    ## because some sites have not randomized any subjects.  Such sites
    ## are counted in # sites but not in site-months
    if(pclose) {
      months <- as.numeric(difftime(closeDate, rdate, units='days')) /
        (365.25 / 12)
      mx <- function(x) if(! length(x) || all(is.na(x))) NA
       else max(x, na.rm=TRUE)
      maxs       <- tapply(months, Site, mx)
      sitemonths <- sum(maxs, na.rm=TRUE)
      z <- c(z, g(max(months, na.rm=TRUE), 1),
                g(sitemonths, 1),
                g(mean(maxs, na.rm=TRUE), 1),
                g(nrand / sitemonths, 2))
      k <- c(k, paste('Months from first subject randomized (',
                      format(min(rdate, na.rm=TRUE)), ') to ',
                      format(closeDate), sep=''),
                'Site-months for sites randomizing',
                'Average months since a site first randomized',
                'Subjects randomized per site per month')
    }
  }
  if(studynos && length(z)) {
    z <- data.frame(Number=z, Category=k)
    u <- latex(z, file=file, append=TRUE, rowname=NULL,
               col.just=c('r','l'), where='!htbp',
               label=paste(panel, 'studynos', sep='-'),
               caption='Study Numbers')
  }

  ## axis.Date when given a sequence not on Jan 1 boundaries did not
  ## place axis labels at correct location
  axisDate <- function(dr) {
    cdr <- as.character(dr)
    yr  <- substring(cdr, 1, 4)
    outer <- as.Date(c(paste(yr[1], '01-01', sep='-'),
                       paste(as.numeric(yr[2]) + 1, '01-01', sep='-')))
    dseq <- seq(outer[1], outer[2], by='year')
    if(length(dseq) > 3) dseq <- dseq[- c(1, length(dseq))]
    short <- difftime(dr[2], dr[1], units='days') < 550
    axis(1, at=as.numeric(dseq),
         labels=if(! short) substring(dseq, 1, 4) else FALSE)
    dseq <- seq(dr[1], dr[2], by='month')
    mo   <- as.numeric(substring(dseq, 6, 7))
    mo   <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',
              'Nov','Dec')[mo]
    mo   <- ifelse(mo == 'Jan', substring(dseq, 1, 4), mo)
    axis(1, at=as.numeric(dseq), 
         labels=if(short) mo else FALSE,
         tcl = 0.5 * par('tcl'), cex.axis=0.6, las=3)
  }
  
  ## For each date variable in Y, make a cumulative frequency chart and
  ## optionally zoomed-in chart
  ## If target sample size is present, add that as line graph to chart
  for(j in 1 : nY) {
    y   <- Y[[j]]
    nam <- namY[j]
    lab <- ylabs[j]
    sumnna <- sum(! is.na(y))
    y <- y[! is.na(y)]
    target <- if(! length(names(targetN))) targetN else targetN[[nam]]
    dtarget <- targetDate
    if(length(target) && min(target) > 0) {
      target <- c(0, target)
      dtarget <- c(dr[1], dtarget)
    }
    lb <- sprintf('%s-cumulative-%s', panel, lab)
    shortcap <- sprintf("Subjects %s over time", lab)
    cap <- if(length(target))
             sprintf('.  The solid back line depicts the cumulative frequency.  The dotted line represent targets.', lab) else ''
    pzoom <- length(zoom) > 0
    if(pzoom) {
      zoom <- as.Date(zoom)
      cap <- paste(cap, sprintf(
        'The plot is zoomed to show %s--%s in the right panel.  The zoomed interval is depicted with vertical grayscale lines in the left panel',
        zoom[1], zoom[2]))
    }

    longcap <- paste(shortcap, cap, '~\\hfill\\lttc', sep = '')

    startPlot(lb, h=h, w=w * (1 + 0.75 * pzoom), lattice=FALSE)
    par(mfrow=c(1, 1 + pzoom), mar=c(4,3.5,2,1))
    Ecdf(as.numeric(y), what='f', xlab=sprintf('Date %s', upFirst(lab)),
         ylab='Cumulative Number',
         subtitles=FALSE, axes=FALSE,
         xlim=as.numeric(dr),
         ylim=c(0, if(length(target)) max(length(y), target)
          else if(lab == 'enrolled' && length(enrollmax)) enrollmax
          else length(y)))
    axis(2)
    axisDate(dr)
    if(length(target)) lines(dtarget, target, lty=3, lwd=1)
    box(lwd=.75, col=gray(.4))
    
    if(pzoom) {
      abline(v=as.numeric(zoom), col=gray(.85))
      Ecdf(as.numeric(y), what='f', xlab=sprintf('Date %s', upFirst(lab)),
           ylab='Cumulative Number',
           subtitles=FALSE, axes=FALSE,
           xlim=zoom,
           ylim=c(0, if(length(target)) max(sum(y <= zoom[2], na.rm=TRUE),
             max(target[dtarget <= zoom[2]])) else
             sum(y <= zoom[2], na.rm=TRUE)))
      axis(2)
      axisDate(zoom)
      if(length(target)) lines(dtarget, target, lty=3, lwd=1)
      box(lwd=.5, col=gray(.4))
    }

    endPlot()
    ltt(switch(lab, enrolled=c(enrolled=sumnna),
               randomized=c(enrolled=sumnna, randomized=sumnna)), 'lttc')
    putFig(panel = panel, name = lb, caption = shortcap,
           longcaption = longcap)
  }

  ## Extended box plots of time to randomization for randomized subjects
  if(penroll && prandomize && (pregion || pcountry)) {
    x1 <- if(pregion)  X[[region]]
    x2 <- if(pcountry) X[[country]]
    lb <- sprintf('%s-timetorand', panel)
    startPlot(lb, h=hb, w=wb, lattice=FALSE)
    y <- as.numeric(difftime(Y[[j]], Y[[enroll]], units='days'))
    use <- TRUE
    coexcl <- 0
    if(pcountry && minrand > 0) {
      ## Exclude countries randomizing fewer than minrand subject
      nrn <- tapply(y, x2, function(x) sum(! is.na(x)))
      if(any(nrn < minrand)) {
        coexcl <- sum(nrn < minrand)
        countrieskeep <- names(nrn)[nrn >= minrand]
        use <- x2 %in% countrieskeep
      }
    }
    form <- if(length(x1) && length(x2)) x2 ~ y | x1
     else if(length(x1)) x1 ~ y
     else if(length(x2)) x2 ~ y
     else x2 ~ 1
    print(bwplot(form, panel=panel.bpplot, xlab='Days to Randomization',
                 subset=use,
                 scales=list(y='free', rot=c(0,0)),
                 violin=TRUE,
                 violin.opts=list(col=adjustcolor('blue', alpha.f=.35),
                                  border=FALSE)))
    endPlot()
    Days <- y
    form <- if(length(x1)) Days ~ x1
       else if(length(x2)) Days ~ x2
       else Days ~ 1
    popname <- '\\poptabledaysrand'
    cat(sprintf('\\def%s{\\protect\n', popname), file=file, append=TRUE)
    rddata <- data.frame(Days)
    if(length(x1)) rddata$x1 <- x1
    if(length(x2)) rddata$x2 <- x2
    rddata <- subset(rddata, ! is.na(Days))
    S <- summaryM(form, data=rddata, test=FALSE)
    z <- latex(S, table.env=FALSE, file=file, append=TRUE, prmsd=TRUE,
               middle.bold=TRUE, center='none', round=1, insert.bottom=FALSE)
    cat('}\n', file=file, append=TRUE)
    popsize <- if(length(S$group.freq) > 2) 'full' else 'mini'
    legend <- attr(z, 'legend')
    legend <- if(! length(legend)) ''
     else paste('. ', paste(legend, collapse='\n'), sep='')

    excc <- if(coexcl > 0) paste('.', coexcl, 'countries with fewer than',
                                 minrand, 'randomized subjects are not shown.')
    else ''
    putFig(panel=panel, name=lb,
           longcaption=paste('\\protect\\eboxpopup{Extended box} plots and violin plots showing the distribution of days from enrollment to randomization',
             excc, '~\\hfill\\lttc', sep=''),
           caption='Days from enrollment to randomization',
           tcaption='Days from enrollment to randomization',
           tlongcaption=paste('Days from enrollment to randomization',
             legend, sep=''),
           poptable=popname, popfull=popsize == 'full')
  }
  
  ## Chart number of subjects enrolled/randomized/... and other descriptors
  ## by right-hand variables
  if(nX == 0) return(invisible())
  
  if(psite) {
    lb <- sprintf('%s-subjpersite', panel)

    startPlot(lb, h=h, w=min(7.75, nY * w), mfrow=c(1, nY),
              ps=8, lattice=FALSE)
    for(j in 1 : nY) {
      y <- X[[site]]
      y[is.na(Y[[j]])] <- NA
      lab  <- ylabs[j]
      clab <- capitalize(lab)
      nn   <- table(table(y))
      plot(as.numeric(names(nn)), as.numeric(nn),
           xlab=sprintf('Number of Subjects %s', clab),
           ylab='Number of Sites')
    }
    endPlot()
    if(nY > 1) lab <- ''
    putFig(panel=panel, name=lb,
           longcaption=sprintf('Number of sites having the given number of subjects %s~\\hfill\\lttc', lab),
           caption=sprintf('Number of sites $\\times$ number of subjects %s', lab))
  }

  ## Start with counts of subjects by non-site grouping variables
  ## Compute number of non-site right-hand variables
  ns <- setdiff(1 : nX, site)
  dat <- list()
  if(pregion)  dat$x1 <- X[[region]]
  if(pcountry) dat$x2 <- X[[country]]
  if(psite)    dat$x3 <- X[[site]]   ## new
  if(length(ns) > 2) {
    more <- setdiff(ns, c(region, country))
    k <- 2
    for(l in more) {
      k <- k + 1
      dat[[paste('x', k, sep='')]] <- X[[l]]
    }
  }
  form <- if(length(ns)) {
    xvars <- paste(paste('x', 1 : length(ns), sep=''), collapse=' + ')
    paste('y ~', xvars)
  } else if(psite) 'y ~ x3'
  else 'y ~ 1'
  form <- as.formula(form)

  by <- paste(xlabs[ns], collapse=' and ')
  types <- c('count',
             if(psite && cr) 'sites',
             if(penroll    && prandomize) 'fracrand',
             if(prandomize && pclose && cr) 'permonth',
             if(prandomize && pclose && psite && cr)
              'persitemonth')

  np <- nY * sum(c('count', 'sites') %in% types) +
             sum(c('fracrand', 'permonth', 'persitemonth') %in%
                 types)
  mf <- if(np == 1) c(1, 1) else if(np == 2) c(1, 2) else c(2, 2)
  mc <- mf[2]
  pages <- ceiling(np / prod(mf))
  width <- if(mf[2] == 1) 3.5 else 7.0
  ip    <- 0
  page  <- 0
  ended <- FALSE
  scap  <- if(psite) 'Subject and site counts'
   else 'Subject counts'
  ## if('fracrand' %in% types) scap <- paste(scap, 'and fraction randomized')
  cap  <- if(psite) 'Counts of numbers of subjects and numbers of sites'
   else 'Counts of numbers of subjects'
  ## if('fracrand' %in% types) cap <- paste(cap, 'and fraction randomized')
  cap <- paste(cap, '~\\hfill\\lttc', sep='')
  for(type in types) {
    whichy <- if(type %in% c('fracrand', 'permonth', 'persitemonth'))
      randomize else 1 : nY
    if(length(whichy)) for(j in whichy) {
      ip <- ip + 1
      if(ip == 1) {
        page <- page + 1
        lb <- if(pages == 1) sprintf('%s-count', panel)
         else sprintf('%s-count-%s', panel, page)
        ## Compute the number of rows in the current page
        r <- if(page < pages) mf[1]
         else if(np %% 4 == 0) 2
         else ceiling((np %% 4) / 2)
        height <- min(9, hdot * r)
        startPlot(lb, h=height, w=width, mfrow=c(r, mf[2]),
                  ps=if(r == 2) 10 else 8, lattice=FALSE)
        ended <- FALSE
      }
      gg <- function(x) length(unique(x[! is.na(x)]))
      
      if(type %in% c('permonth', 'persitemonth')) {
        ## Get country if there, otherwise region
        group <- if(pcountry) as.character(X[[country]])
         else    if(pregion)  as.character(X[[region]])
        ## Get more major grouping if present otherwise use above
        mgroup <- if(pregion) as.character(X[[region]]) else group

        ## Get enrollment date if present, otherwise use rand. date
        k <-  if(penroll)     enroll
         else if(prandomize)  randomize
         else 1
        months <- as.numeric(difftime(closeDate, Y[[k]], units='days')) /
          (365.25 / 12)
        ## Find maximum months on board for each group
        ## E.g. longest elapsed time within a country
        gmonths <- tapply(months, group, max, na.rm=TRUE)
        ## Find the maximum elapsed time over groups within major groups
        ## E.g. longest time for any country within that region
        mmonths <- tapply(months, mgroup, max, na.rm=TRUE)

        ## Create a major group lookup object given group
        tab <- subset(as.data.frame(table(group, mgroup)), Freq > 0)
        mg        <- as.character(tab$mgroup)
        names(mg) <- as.character(tab$group)

        ## For site-month calculation compute the maximum elapsed time
        ## per site, then sum that over all sites within a group
        ## Assume sites are unique over countries, regions
        if(type == 'persitemonth') {
          ## For each site lookup group
          tab <- subset(as.data.frame(table(Site, group)), Freq > 0)
          gr        <- as.character(tab$group)
          names(gr) <- as.character(tab$Site)
          maxs <- tapply(months, Site, max, na.rm=TRUE)

          ## Starting with only one record per site with that site's
          ## maximum time, sum the elapsed months within each group
          gsitesum <- tapply(maxs, gr[names(maxs)], sum, na.rm=TRUE)
          ## Similar over region
          msitesum <- tapply(maxs, mg[gr[names(maxs)]], sum, na.rm=TRUE)
          ## Spread to all subjects
          y <- cbind(randomized = ! is.na(Y[[randomize]]),
                     mmonths    = msitesum[mg[group]],
                     gmonths    = gsitesum[group])
          yy <- y[group == 'US',]
        }
          else y <- cbind(randomized = ! is.na(Y[[randomize]]),
                          mmonths    = mmonths[mg[group]],
                          gmonths    = gmonths[group])
        mg <- function(y) sum(y[, 1]) / y[1, 2]
        gg <- function(y) sum(y[, 1]) / y[1, 3]
      }

      switch(type,
             count = { y <- ! is.na(Y[[j]]); fun <- sum },
             sites = { y <- X[[site]]; y[is.na(Y[[j]])] <- NA; fun <- gg },
             fracrand     = { y <- ! is.na(Y[[j]]); fun <- mean },
             permonth     = { },
             persitemonth = { } )

      lab <- ylabs[j]
      clab <- capitalize(lab)
      dat$y <- y
      fmt <- function(x) format(round(x, 2))
      if(type %in% c('permonth', 'persitemonth'))
        summaryD(form, fun=gg, funm=mg, data=dat, vals=TRUE, fmtvals=fmt,
                 xlab=switch(type,
                   permonth     = 'Number Randomized Per Month',
                   persitemonth = 'Number Randomized Per Site Per Month'))
      else summaryD(form, fun=fun, data=dat, vals=TRUE,
                    fmtvals = fmt,
                    ylab = if(psite && ! length(ns)) 'Site',
                    xlab=switch(type,
                      count=sprintf('Number of Subjects %s',      clab),
                      sites=sprintf('Number of Sites That %s',    clab),
                      fracrand=sprintf('Fraction of Subjects %s', clab)))
      if(ip == 4) {
        endPlot()
        putFig(panel=panel, name=lb, longcaption=cap, caption=scap)
        ended <- TRUE
        ip <- 0
      }
    }
  }
  if(!ended) {
    endPlot()
    putFig(panel=panel, name=lb, longcaption=cap, caption=scap)
  }
  invisible()
}
