

## stats: Freq und PercTable ====

# replaced in version 0.99.11
Freq <- function(x, breaks = hist(x, plot = FALSE)$breaks, include.lowest = TRUE,
                 ord = c("level", "desc", "asc", "name"),
                 useNA = c("no", "ifany", "always"), ...){

  # check if x is a vector (do not use is.vector()!!!)
  if(!(is.atomic(x) || is.list(x))) stop("'x' must be a vector")

  if(inherits(x, "table")){
    tab <- x

  } else {

    if(is.numeric(x) || IsDate(x)){
      x <- cut(x, breaks = breaks, include.lowest = include.lowest,
               ordered_result = TRUE, ...)
    }

    tab <- table(x, useNA = useNA)
  }

  # how should the table be sorted, by name, level or frq? (NULL means "desc")
  switch(match.arg(ord, c("level", "desc", "asc", "name")),
         level  = {  }
         , name   = { tab <- tab[rownames(tab)] }
         , asc    = { tab <- sort(tab) }
         , desc   = { tab <- -sort(-tab) }
  )

  ptab <- prop.table(tab)
  names(tab)[is.na(names(tab))] <- "<NA>"

  z <- data.frame(level = names(tab),
                  freq = as.vector(tab[]), perc = as.vector(ptab[]),
                  cumfreq = cumsum(tab[]), cumperc = cumsum(ptab[]))

  rownames(z) <- NULL # enumerate from 1:nrow(z)
  class(z) <- c("Freq", "data.frame")
  return(z)

}



print.Freq <- function(x, digits=NULL, ...) {

  # print as data.frame if something was changed
  if(ncol(x) != 5) {
    print.data.frame(x, digits=digits, ...)
  } else {

    afmt <- .fmt_abs()
    pfmt <- .fmt_per(digits=digits)

    # object x comes as list lacking an as.data.frame option...
    print(data.frame(level=x$level, freq=Format(x$freq, fmt=afmt), perc=Format(x$perc, fmt=pfmt),
                     cumfreq=Format(x$cumfreq, fmt=afmt), cumperc=Format(x$cumperc, fmt=pfmt)),
          print.gap = InDots(..., arg="print.gap", default=2))
  }
}


# ToDo: PercTable
# param cumul.count show cumulative frequencies?
# param cumul.pct show cumulative percentage?
# param total.name a string containing footer label (defaults to "Sum")
# Drop unused levels
# useNA ?
# expected values, residuals, standardized residuals



PercTable <- function (...) UseMethod("PercTable")

PercTable.default <- function (x, y = NULL, ...) {

  # all dot arguments
  dot.args <- match.call(expand.dots=FALSE)$...
  # the dot arguments which match PercTable.table
  pt.args <- dot.args[names(dot.args) %in% names(formals(PercTable.table))]
  # the dot arguments which DO NOT match PercTable.table
  tab.args <- dot.args[names(dot.args) %nin% names(formals(PercTable.table))]

  if(is.null(y)){
    tab <- do.call("table", append(list(x), tab.args) )
  } else {
    tab <- do.call("table", append(list(x, y), tab.args) )
  }
  do.call( "PercTable", append(list(tab=tab), pt.args) )

}

# PercTable.data.frame <- function(x, ...){  sapply(x, PercTable, ...) }
PercTable.matrix <- function(x, ...){  PercTable(as.table(x), ...) }


PercTable.formula <- function(formula, data, subset, na.action, ...)
{
  # this is taken basically from wilcox.test.formula

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")

  DATA <- list(table(mf))
  do.call("PercTable", c(DATA, list(...)))
}



PercTable.table <- function(tab, row.vars=NULL, col.vars = 2, justify = "right"
                            , freq=TRUE, rfrq="100",
                            expected = FALSE, residuals = FALSE, stdres = FALSE, margins = NULL
                            , ...) {

  # labels = c("Sum", "freq", "perc", "p.row", "p.col"),

  # example:
  # tab <- table(d.pizza[,c("driver","operator")])
  # PercTable(tab, rfrq="110", margins=c(1,2))

  # create prop tables and format them
  fmt.tab <- function(x, perc, w) {

    if(perc==1) {
      px <- addmargins(prop.table(addmargins(x, 1), 1), 2)
      if(1 %nin% margins) px <- px[,-ncol(px)]
      if(2 %nin% margins) px <- px[-nrow(px),]
      class(px) <- "table"
    } else if(perc==2) {
      px <- addmargins(prop.table(addmargins(x, 2), 2), 1)
      if(1 %nin% margins) px <- px[,-ncol(px)]
      if(2 %nin% margins) px <- px[-nrow(px),]
      class(px) <- "table"
    } else {
      px <- prop.table(x)
      if(!is.null(margins)) px <- addmargins(px, if(length(dim(x))==1) {1} else {3 - margins} )
    }

    # get the percent format from global option
    px[] <- Format(px, fmt=.fmt_per())

    # set 100% margins to some zero value
    # but only if main percentages are requested
    if(perc==1 & (1 %in% margins)) px[, ncol(px)] <- zero
    if(substr(rfrq, 1, 1)=="1")
      if(perc==2 & (1 %in% margins)) px[, ncol(px)] <- zero
    if(substr(rfrq, 1, 1)=="1")
      if(perc==1 & (2 %in% margins)) px[nrow(px), ] <- zero
    if(perc==2 & (2 %in% margins)) px[nrow(px), ] <- zero

    px
  }


  # set zero element
  zero <- "."

  tlst <- list(freq=tab)

  # overwrite percents if only 1-dim table
  if(length(dim(tab)) == 1) rfrq <- paste(sum(as.numeric(rfrq) > 0), "00", sep="")

  if(unlist(strsplit(rfrq, NULL))[1] == "1")
    tlst[["perc"]] <- fmt.tab(tab, perc=0)
  if(unlist(strsplit(rfrq, NULL))[2] == "1")
    tlst[["p.row"]] <- fmt.tab(tab, perc=1)
  if(unlist(strsplit(rfrq, NULL))[3] == "1")
    tlst[["p.col"]] <- fmt.tab(tab, perc=2)

  # flip 1 to 2 and 2 to 1 in margins with: 3 - margins
  if(!is.null(margins)) tlst[["freq"]] <- addmargins(tab, if(length(dim(tab))==1) {1} else {3 - margins})

  # format tab as.character
  tlst[["freq"]][] <- Format(tlst[["freq"]], fmt=.fmt_abs())
  if(freq == FALSE) tlst[["freq"]] <- NULL

  suppressWarnings(r.chisq <- chisq.test(tab))

  if(expected == TRUE) {
    tlst[["exp"]] <- Format(r.chisq$expected, fmt=.fmt_num())
    if(1 %in% margins) tlst[["exp"]] <- cbind(tlst[["exp"]], Sum=zero)
    if(2 %in% margins) tlst[["exp"]] <- rbind(tlst[["exp"]], Sum=zero)
  }
  if(residuals == TRUE){
    tlst[["res"]] <- Format(r.chisq$residuals, fmt=.fmt_num())
    if(1 %in% margins) tlst[["res"]] <- cbind(tlst[["res"]], Sum=zero)
    if(2 %in% margins) tlst[["res"]] <- rbind(tlst[["res"]], Sum=zero)
  }
  if(stdres == TRUE) {
    tlst[["stdres"]] <-  Format(r.chisq$stdres, fmt=.fmt_num())
    if(1 %in% margins) tlst[["stdres"]] <- cbind(tlst[["stdres"]], Sum=zero)
    if(2 %in% margins) tlst[["stdres"]] <- rbind(tlst[["stdres"]], Sum=zero)
  }

  if(length(tlst) == 1){
    ftab <- ftable(tlst[[1]])

  } else {
    # if(length(dim(tab)) > 1){
    #   if(is.null(vsep))
    #     vsep <- length(tlst) > 1
    #   if(vsep)         # insert a separator line
    #     tlst[[""]] <- matrix("", nrow=nrow(tlst[[1]]), ncol=ncol(tlst[[1]]))
    # }

    # build a table array, such as to be able to pass it to ftable afterwards...
    if(length(dim(tab))==1){
      ma <- do.call("cbind", tlst)
    } else {
      ma <- do.call("Abind", c(tlst, along = 3))
    }
    ftab <- ftable(ma, col.vars=col.vars, row.vars=row.vars)

  }

  justify <- match.arg(justify, c("left","right"))
  if(justify == "right")
    # align the whole stuff to the right
    ftab[] <- StrAlign(ftab, "\\r")

  names(attr(ftab, "row.vars"))[1] <- names(dimnames(tab))[1]
  if(length(names(attr(ftab, "row.vars"))) == 2)
    names(attr(ftab, "row.vars"))[2] <- ""
  names(attr(ftab, "col.vars")) <- names(dimnames(tab))[2]

  res <- list(ftab=ftab, tlst=tlst)
  vsep <- InDots(..., arg="vsep", default=ifelse(length(dim(tab)) == 1, FALSE, NA))
  if(!is.na(vsep)) res[["vsep"]] <- vsep

  class(res) <- c("PercTable")

  return(res)

}


print.PercTable <- function(x, vsep=NULL, ...){

  vsep <- Coalesce(vsep, x[["vsep"]], (length(x[["tlst"]]) > 1))

  x <- x[["ftab"]]
  cat(paste(c(rep(" ", times=max(nchar(c(names(attr(x, "row.vars")), attr(x, "row.vars")[[1]])), na.rm=TRUE) +
                    ifelse(length(attr(x, "row.vars")) == 2, max(nchar(attr(x, "row.vars")[[2]]), na.rm=TRUE) + 2, 0) + 1),
              names(attr(x, "col.vars"))), collapse=""), sep="\n")

  names(attr(x, "col.vars")) <- NULL

  txt <- capture.output(print(x, ...))
  if(vsep)
    txt[StrLeft(txt,1) != " "][-1] <- paste("\n", txt[StrLeft(txt,1) != " "][-1], sep="")

  cat(txt, sep="\n")

}



Margins <- function(tab, ...){
  lst <- lapply(1:length(dim(tab)),
                function(i) Freq(margin.table(tab, i), ...))
  names(lst) <- names(dimnames(tab))
  lst
}



ExpFreq <- function(x, freq = c("abs", "rel")) {

  # returns the expected frequencies of a table assuming independence

  # this is a copy of independence_table {vcd}
  # by David Meyer David.Meyer@R-project.org

  if (!is.array(x))
    stop("Need array of absolute frequencies!")
  frequency <- match.arg(freq)
  n <- sum(x)
  x <- x/n
  d <- dim(x)
  margins <- lapply(1:length(d), function(i) apply(x, i, sum))
  tab <- array(apply(expand.grid(margins), 1, prod), d, dimnames = dimnames(x))
  if (frequency == "rel")
    tab
  else tab * n
}


CollapseTable <- function (x, ...) {

  nargs <- length(args <- list(...))
  if (!nargs)
    return(x)

  if (inherits(x, "ftable"))
    x <- as.table(x)

  if (inherits(x, "table")) {
    tvars <- names(dimnames(x))
    x <- as.data.frame.table(x)
    freq <- x[, "Freq"]
  } else {
    stop("Argument must be a table or ftable object")
  }

  names <- names(args)

  for (i in 1:nargs) {
    vals <- args[[i]]
    nm <- names[[i]]
    if (any(nm == tvars))
      levels(x[[nm]]) <- vals
    else warning(nm, " is not among the x variables.")
  }

  xtabs(as.formula(paste("freq ~", paste(tvars, collapse = "+"))),
        data = x)

}



###

