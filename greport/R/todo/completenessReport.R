#' Completeness Report
#'
#' Generate report with information about the completeness of the data.
#'
#' details
#'
#' @param data data.frame. Data to report.
#' @param vars character vector. Column names to use.
#' @param panel character. Name of panel.
#' @param Time numeric vector. Times for each record.
#' @param times numeric vector. Subset of times to use.
#' @param longPanel character. Long name of panel.
#' @param frac numeric. Ratio used to determine if a dotplot for completeness
#' should be generated for each variable.  Defaults to 0.95, that is, if the variable
#' with the most missing data has 5% more missing than the variable with the least
#' missing data, generate the dotplot.
#' @param h numeric. Height of plot, passed to \code{\link{startPlot}}.
#' @param cex numeric. Text size within plot, passed to \code{\link[Hmisc]{Dotplot}}.
#' @param fullCaption logical. THIS PARAMETER IS NOT USED.
#' @param append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @export

completenessReport <- function(data, vars, 
                               panel, Time, times,
                               longPanel=panel, frac=0.95,
                               h=4 * (nv < 15) + 5 * (nv >= 15), cex=1,
                               fullCaption=FALSE, append=TRUE)
{
  if(!exists('compFullCaptionDone'))
    storeTemp(FALSE, 'compFullCaptionDone')
  
  if(!append) cat('', file='gentex/completeness-key.tex')

  needkey <- TRUE

  vars <- unlist(vars)
  xlab <- 'Number of Non-Missing Values'
  timeUsed <- !missing(Time)
  if(timeUsed && !missing(times)) {
    s <- Time %in% times
    data <- data[s,]
    Time <- Time[s]
  }
  fn <- paste('complete', panel, sep='-')
  nv <- length(vars)
  pl <- TRUE
  
  lcap <- paste('Completeness of ', longPanel, '.', sep='')
  if(!timeUsed) {
      n <- sapply(data[vars], function(y)sum(!is.na(y)))
      ce <- combineEqual(n)
      n <- ce$x
      r <- range(n)
      if(r[1] / r[2] < frac)
        {
          variable <- factor(names(n),names(n))
          variable <- reorder(variable, n)
#          variable <- reorder.factor(variable, n)
          startPlot(fn, h=h)
          print(Dotplot(variable ~ n, xlab=xlab, cex=cex))
        }
      else pl <- FALSE
  } else {
    n  <- rowsum(1 - is.na(data[vars]), Time)
    ce <- combineEqual(n)
    n  <- ce$x
    nvarpl <- ncol(n)
    if(nvarpl < 2) {
        startPlot(fn, h=3)
        plot(as.numeric(dimnames(n)[[1]]), n[,1],
            xlab=label(Time),
            ylab='Number of Values', type='b')
        lcap <- paste(lcap,
                    ' The $y$-axis displays the number of values measured for',
                    ' \\texttt{',
                    if(length(ce$defs)) ce$defs else vars[1],
                    '}. ', sep='')
        needkey <- FALSE
    } else {
        m <- reShape(n)
        variable <- reorder(factor(m$colvar), m$n)
        times <- sort(unique(Time))
        Time <- factor(m$rowvar, times, paste(label(Time),times))
        n <- m$n
        startPlot(fn, h=h)
        print(Dotplot(variable ~ n | Time,
                    xlab=xlab, cex=cex))
        if(length(ce$codes))
        lcap <-
            paste(lcap,
                '  Letters in parentheses indicate groups of variables ',
                'having the same number of values measured, ',
                'defined elsewhere.\n', sep='')
      if(!compFullCaptionDone) {
          lcap <- paste(lcap,
                        'Variables are sorted by the average number of',
                        'non-missing values over time.')
          storeTemp(TRUE, 'compFullCaptionDone')
      }
    }
  }
  if(pl) {
      endPlot()
      putFig('completeness', fn,
             paste('Completeness of', longPanel),
             lcap, append=append)
      if(length(ce$codes) & needkey) {
        cat('\nCodes used in Figure~\\ref{fig:',fn,'} are as follows:',
            '{\\smaller[2]',
            paste(paste('\\textbf{',ce$codes,'}:',
                        '\\texttt{\\textbf{',ce$defs,'}}',sep=''),
                  collapse='; '),
            '.}\n\n', sep='',
            file=file.path(TexDirName(), 'completeness-key.tex'), append=TRUE)
      }
  } else {
    cat('For the',length(vars),longPanel,
           'variables, the number of subjects having values entered',
           if(r[1]==r[2])paste('was ', r[1],'.\n', sep='') else
           paste('ranged from ',r[1],' to ', r[2],'.\n', sep=''),
           file=file.path(TexDirName(), 'completeness.tex'), append=append)
  }
  invisible()
}
