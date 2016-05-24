# make.labels
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/util_r/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#                            http://paidforeveryone.wordpress.com/

# convert charcter vector names into labels to show
#
# convert = how to generate labels from names
#           default=NULL, equal to names
#           TRUE to abbreviate by alphabet letter sequences
#           FALSE to draw no labels
#           character vector to give explicit labels
#           single character to use identical character
#           NA in vector can be used to remove the item

make.labels <- function(names, convert=NULL, blank='') {
  # for names=c('setosa,'versicolor','virginica')
  #
  # convert=NULL
  # will generate labels, c('setosa,'versicolor','virginica')
  #
  # convert=TRUE
  # will generate labels, c('a','b','c')
  #
  # convert=FALSE
  # will generate labels, c('','','')
  #
  # convert=c('se','ve','')
  # will generate labels, c('se','ve','')
  # convert=c('se','ve')
  # will generate labels, c('se','ve',''), too, with keeping NA
  #
  # convert='*'
  # will generate labels, c('*','*','*')

  n <- length(names)
  use.strict <- is.null(convert)
  convert <-
    if(length(convert) == 1) rep(convert, n) else
      convert[1L:n]

  how.to.convert <-
    if(use.strict) {
      rep('s', n)          # strict
    } else if(is.logical(convert)) {
      ifelse(convert, 
             'g',          # generates internally
             'n')          # blank
    } else {
      ifelse(is.na(convert), 
             NA,           # blank, or remove
             'e')          # specified explicitly
    }
  if(n > 26) {
    find.g <- how.to.convert %in% 'g'
    find.g[1:26] <- FALSE
    how.to.convert[find.g] <- 's'
  }

  labels <- 
    sapply(as.list(1L:n),
           function(i) switch(as.character(how.to.convert[i]),
                              s=as.character(names[i]),
                              g=letters[i],
                              e=as.character(convert[i]),
                              blank)
           )

  not.na <- !is.na(how.to.convert)
  levels <- sort(unique(labels[not.na]))
  single.char <- (length(levels) <= 1)
  names(labels) <- names

  list(labels=labels, names=names, not.na=not.na, 
       levels=levels, single.char=single.char,
       use.strict=use.strict, how.to.convert=how.to.convert)
}

make.legend.labels <- function(..., col=NULL, sep=':', na.rm=TRUE, 
                               single.rm=TRUE, dup.rm=TRUE, 
                               verbose=FALSE) {
  # dots accepts one or more of the output of make.labels
  # should have same numbers
  #
  # na.rm = remove by not.na
  # single.rm = remove if all labels are same
  # dup.rm = remove duplication by unique function
  # col = color vector
  # sep = separator between labels

  pars <- list(...)
  if(length(pars) == 0) return(NULL)

  par1 <- pars[[1]]
  n <- length(par1$labels)

  flipped <- lapply(as.list(1L:n), function(i) {
               list(labels=sapply(pars, function(x) x$labels[i]),
                    names=sapply(pars, function(x) x$names[i]),
                    single.char=sapply(pars, function(x) x$single.char),
                    how.to.convert=sapply(pars, function(x) x$how.to.convert[i]))
             })
  if(verbose) print(flipped)

  LABELING <- function(x) {
    paste(c(
            (if(dup.rm) unique else function(a) a)(x), 
            ''), 
          collapse=sep)
  }
  label <- 
    sapply(flipped, 
           function(x) LABELING(x$labels[!(single.rm & x$single.char) & 
                                x$how.to.convert != 's']))

  legend <- paste(label, par1$names)
  col <- rep(col, n)[1L:n]
  names(legend) <- names(col) <- par1$names

  if(na.rm) {
    legend <- legend[par1$not.na]
    col <- col[par1$not.na]
  }

  list(legend=legend, col=col, label=label, 
       name=par1$names, not.na=par1$not.na)
}

