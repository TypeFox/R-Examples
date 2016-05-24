latex.fdt_cat <- function(x,
                          columns=1:6,
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
    header <- colnames(x)[columns]

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

  y <- x[columns][-1]

  y <- cbind(x[1],
             y)

  y <- sapply(y,
              as.character)

  if(is.null(dim(y))){

   y <- sapply(y,
              function(x) paste(' & ',
                                x,
                                sep=''))

   y <- paste(paste(y,
                    collapse=''),
              ' \\\\',
              sep='')

  }
  else {

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
  
  }

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

  class(res) <- c('latex_fdt',
                  'list')

  return(res)
}
