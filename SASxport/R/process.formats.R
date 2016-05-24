## Transform SAS 'PROC CONTENTS' dataset into a form useful for
## converting raw SAS objects to/from the appropriate R objects.

process.formats <- function(finfo)
  {
    if(is.null(finfo)) return( list() )
    ## Remove leading $ from char format names
    ##  fmtname <- sub('^\\$','',as.character(finfo$FMTNAME))
    fmtname <- as.character(finfo$FMTNAME)
    finfo <- split(finfo[c('START','END','LABEL')], fmtname)
    finfo <- lapply(finfo,
                    function(f)
                    {
                      rb <- function(a)
                        {  # remove leading + trailing blanks
                          a <- sub('[[:space:]]+$', '', as.character(a))
                          sub('^[[:space:]]+', '', a)
                        }

                      st <- rb(f$START)
                      en <- rb(f$END)
                      lab <- rb(f$LABEL)
                      ##j <- is.na(st) | is.na(en)
                      ##  st %in% c('','.','NA') | en %in% c('','.','NA')
                      j <- is.na(st) | is.na(en) | st == '' | en == ''
                      if(any(j)) {
                        warning('NA in code in FORMAT definition; removed')
                        st <- st[!j]; en <- en[!j]; lab <- lab[!j]
                      }

                      if(!all(st==en))
                        stop("Format ranges are not handled.")

                      list(value = all.is.numeric(st, 'vector'),
                           label = lab)
                    })
    finfo
  }
