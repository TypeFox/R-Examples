latex.fdt <- function(x,
                      columns=1:6,
                      round=2,
                      format.classes=TRUE,
                      pattern='%.2f',
                      replace.breaks=TRUE,
                      where='!tbp',
                      caption=NULL,
                      label=NULL,
                      size='', 
                      algtable=c('\\flushleft', '\\centering', '\\flushright'),
                      hline1='\\hline',
                      header=NULL,
                      hline2='\\hline',
                      algclim='l',
                      algfreq='r',
                      hline3='\\hline')
{
  if (is.null(header))
    header <- colnames(x$table)[columns]

  nfreq <- length(header) - 1

  # substitute % in header
  header <- gsub('%',
                 '\\\\%',
                 header)

  header <- paste(' & ',
                  header,
                  sep='')

  header <- c(header,
              '\\\\')

  header <- paste(header,
                  collapse=' ')

  clim <- x$table[1]

  clim <- sapply(clim,
                 as.character)

  right <- x$breaks[4]

  # format Class limits
  if (format.classes) {
    clim <- make.fdt.format.classes(clim,
                                    right,
                                    pattern)
  }

  # replace breaks of Class limits
  if (replace.breaks) {
    # remove all: [, ], ( and ) of Class limits
    clim <- gsub('[][()]',
                 '',
                 clim)

    # substitute ,
    if (right) {
      clim <- gsub(',',
                   ' $\\\\dashv$ ',
                   clim)
    } else {
      clim <- gsub(',',
                   ' $\\\\vdash$ ',
                   clim)
    }
  }

  y <- x$table[columns][-1]

  y <- cbind(clim,
             round(y, round))

  y <- sapply(y,
              as.character)

  y <- apply(y,
             2,
             function(x) paste(' & ',
                               x,
                               sep=''))

  y <- paste(apply(y,
                   1,
                   paste,
                   collapse=''),
             ' \\\\',
             sep='')

  res <- list(start=NULL,
              algtable=match.arg(algtable),
              size=size,
              caption=NULL,
              label=NULL,
              begintabular=NULL,
              hline1=hline1,
              header=header,
              hline2=hline2,
              fdt=y,
              hline3=hline3,
              endtabular=NULL,
              end=NULL)

  res$start <- sub('where',
                   where,
                   '\\begin{table}[where]')

  if(!is.null(caption))
    res$caption <- paste('\\caption{',
                         caption,
                         '}',
                         sep='')

  if(!is.null(label))
    res$label <- paste('\\label{',
                       label,
                       '}',
                       sep='')

  begintabular <- paste(c('\\begin{tabular}{r',
                          algclim,
                          rep(algfreq, nfreq),
                          '}'),
                        collapse='')

  res$begintabular <- begintabular

  res$endtabular <- '\\end{tabular}'

  res$end <- '\\end{table}'

  class(res) <- c('latex.fdt',
                  'list')

  return(res)
}
