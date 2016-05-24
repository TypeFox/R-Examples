as.likertDataFrame.listOfNamedMatrices <-
  function(x, xName=deparse(substitute(x))) {
    names.x <- names(x)
    names(x) <- NULL
    rownames.x <- unname(unlist(lapply(x, rownames)))
    xNULL <- lapply(x, `rownames<-`, NULL)
    xx <- do.call(rbind, xNULL)
    xx <- data.frame(xx,
                     Question=factor(rownames.x, levels=unique(rownames.x)),
                     Subtable=rep(
                       factor(names.x, levels=names.x),
                       sapply(xNULL, nrow)),
                     check.names=FALSE)
    attr(xx,"names.dimnames") <- names(dimnames(x[[1]]))
    xx
  }



as.likertDataFrame <- function(x, xName=deparse(substitute(x)))
  UseMethod("as.likertDataFrame")

as.likertDataFrame.array <- function(x, xName=deparse(substitute(x)))
  as.likertDataFrame(as.listOfNamedMatrices(x, xName=xName))

as.likertDataFrame.matrix <- function(x, xName=deparse(substitute(x))) {
  result <- data.frame(x, Question=rownames(x), check.names=FALSE)
  attr(result, "names.dimnames") <- names(dimnames(x))
  result
}

as.likertDataFrame.data.frame <- as.likertDataFrame.matrix



plot.likert.formula.old <- function(x, data, ..., xName=deparse(sys.call()), main) {
  force(xName)
  if (missing(data)) stop("The data argument is required.", call.=FALSE)
  if (missing(main)) {
    main <- match.call()[1:3]
    names(main)[] <- ""
    main[[1]] <- as.name("likert")
    main <- deparse(main, width.cutoff = 500L)
  }

  switch(length(x),
         stop(paste("formula must be in form\n  Question ~ . | Subtable\nor\n  Question ~ A+B+C | Subtable\nor\n  ~ A+B+C\nor\n  ~ ."), call.=FALSE),
         {
           QuestionName <- ' '
           data[, QuestionName] <- factor(rownames(data), levels=unique(rownames(data)))
           x <- as.formula(paste('` `', deparse(x))) ## length(x) is now 3
         },
         QuestionName <- as.character(x[[2]]))

  ## convert non-'.' RHS to '.', and subset data.frame to match
  ## x3 <- x[[3]]
  ## if ((x3 != ".") && ((length(x3) != 3) || (x3[[2]] != "."))) {
  ##   DOT <- as.call(~ .)[[2]]
  ##   if ((length(x3) > 1) && (x3[[1]] == "|")) {
  ##     CondNames <- as.character(x3[[3]])
  ##     LevelNamesFormula <- x3[[2]]
  ##     x3[[2]] <- DOT
  ##   }
  ##   else {
  ##     CondNames <- NULL
  ##     LevelNamesFormula <- x3
  ##     x3 <- DOT
  ##   }

  x3 <- x[[3]]
  DOT <- as.call(~ .)[[2]]

  ## the assignment of rownames can't be factored out here.
  ## The listOfNamedMatrices could have identical rownames, and then data.frame() would complain.
  if (x3 == DOT) {
    rownames(data) <- data[, QuestionName]
    names.dimnames <- attr(data, "names.dimnames")
    data[, QuestionName] <- NULL
    data <- data.matrix(data[sapply(data, is.numeric)]) ## not redundant, data.matrix converts character columns to NA
    names(dimnames(data)) <- names.dimnames
    return(
      plot.likert.default(data, xName=xName, ..., main=main)
      )
  }

  ## ((x3 != DOT)

  if (length(x3) > 1 && x3[[1]] == '|' && x3[[2]] == DOT) {
    CondNamesFormula <- x3[[3]]
    LevelNamesFormula <- NULL
    LevelNames <- names(data)[sapply(data, is.numeric)]
  }
  else {
    if (length(x3) > 1 && x3[[1]] == '|' && x3[[2]] != DOT) {
      LevelNamesFormula <- x3[[2]]
      x3[[2]] <- DOT
      CondNamesFormula <- x3[[3]]
    }
    else { ## length(x3) == 1 || x3[[1]] != '|'
      LevelNamesFormula <- x3
      x3 <- DOT
      CondNamesFormula <- NULL
    }
  }

  if (is.null(CondNamesFormula))
    CondNames <- NULL
  else {
    CondNamesRaw <- strsplit(deparse(CondNamesFormula), "+", fixed=TRUE)[[1]]
    CondNames <- sub('^[[:space:]]+', '', sub('[[:space:]]+$', '', CondNamesRaw )) ## remove leading and trailing white space, POSIX-style
    CondNames <- sub('^\"', '', sub('\"$', '', CondNames )) ## remove leading and trailing '\"' character
  }
  if (!is.null(LevelNamesFormula)) {
    LevelNamesRaw <- strsplit(deparse(LevelNamesFormula), "+", fixed=TRUE)[[1]]
    LevelNames <- sub('^[[:space:]]+', '', sub('[[:space:]]+$', '', LevelNamesRaw )) ## remove leading and trailing white space, POSIX-style
    LevelNames <- sub('^\"', '', sub('\"$', '', LevelNames )) ## remove leading and trailing '\"' character
  }
  data <- data[c(QuestionName, CondNames, LevelNames)]  ## Subset of data columns actually used


  if (is.null(CondNames)) {
    rownames(data) <- data[, QuestionName]
    data <- data[order(data[, QuestionName]),]
    data[, QuestionName] <- NULL
    return(
      plot.likert(data, xName=xName, ..., main=main)
      )
  }
  else {
    data.CondNames <- data[, CondNames]
    if (!("factor" %in% class(data.CondNames)))
        data.CondNames <- factor(data.CondNames, levels=unique(data.CondNames))
    data.rownames <- split(data[, QuestionName], data.CondNames)  ## fix this
    data[, QuestionName] <- NULL
    data <- data[,sapply(data, is.numeric), drop=FALSE] ## CondNames and maybe ` ` are removed.
    x.split <- split(data, data.CondNames)
    names.dimnames <- attr(data, "names.dimnames")
    data <- as.listOfNamedMatrices(x.split)
    `ordered.rownames<-` <- function(x, value) {
      rownames(x) <- value
      x[order(value), ]
    }
    data <- mapply(`ordered.rownames<-`, data, data.rownames, SIMPLIFY=FALSE)
    data <- lapply(data, data.matrix)
    ## data <- lapply(data, `attr<-`, "names.dimnames", names.dimnames)
    data <- lapply(data, function(z) {
      names(dimnames(z)) <- names.dimnames
      z
    })
  }
  plot.likert.list(data, xName=xName, ..., main=main)
}


if (FALSE) {
tmp <- data.frame(a=c("hundred","thousand","million"),
                  b=factor(c("hundred","thousand","million")),
                  c=factor(c("hundred","thousand","million"), levels=c("hundred","thousand","million")),
                  d=1:3,
                  stringsAsFactors=FALSE)
sapply(tmp, class)
sapply(tmp, levels)

split(tmp$d, tmp$a) ## split coerces a to an alphabetical factor
split(tmp$d, tmp$b) ## b is already an alphabetical factor
split(tmp$d, tmp$c) ## c uses the assigned factor order

likert(a ~ d, data=tmp)
likert(b ~ d, data=tmp)
likert(c ~ d, data=tmp)


tmp2 <- cbind(rbind(tmp, tmp),
              g=rep(c("e","d"), each=3),
              h=factor(rep(c("e","d"), each=3)),
              i=factor(rep(c("e","d"), each=3), levels=c("e","d")))
sapply(tmp2, levels)

likert(c ~ d | g, data=tmp2)
likert(c ~ d | h, data=tmp2)
likert(c ~ d | i, data=tmp2)
}




## if (FALSE) {
## data(NZScienceTeaching)

## likert(Question ~ ., NZScienceTeaching)

## likert(Question ~ . | Subtable, NZScienceTeaching, layout=c(1,2), scales=list(y=list(relation="free")))


## data(ProfChal)
## ## these are the twelve possible structures for formulas
## likert(Question ~ .                          , ProfChal)
## likert(Question ~ .                | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))
## likert(Question ~ Agree                      , ProfChal)
## likert(Question ~ Agree            | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))
## likert(Question ~ Disagree + Agree           , ProfChal)
## likert(Question ~ Disagree + Agree | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))
## likert(         ~ .                          , ProfChal)
## likert(         ~ .                | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))
## likert(         ~ Agree                      , ProfChal)
## likert(         ~ Agree            | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))
## likert(         ~ Disagree + Agree           , ProfChal)
## likert(         ~ Disagree + Agree | Subtable, ProfChal, layout=c(1,6), scales=list(y=list(relation="free")))

## likert(Question ~ . | Subtable, ProfChal, as.percent=TRUE, layout=c(1,6), scales=list(y=list(relation="free")))


## data(SFF8121)
## likert(Question ~ . | Subtable, as.likertDataFrame(SFF8121), layout=c(2,1))
## }
