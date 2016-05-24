# based on xtable package

##' @export
##' @method ascii ts
ascii.ts <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 0, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){

  if (inherits(x, "ts") && !is.null(ncol(x))) {
    tp.1 <- trunc(time(x))
    tp.2 <- trunc(cycle(x))
    day.abb <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
    ROWNAMES <- switch(frequency(x),
        tp.1,
        "Arg2", "Arg3",              ## Dummy arguments
        paste(tp.1, c("Q1", "Q2", "Q3", "Q4")[tp.2], sep=" "),
        "Arg5", "Arg6",
        paste("Wk.", tp.1, " ", day.abb[tp.2], sep=""),
        "Arg8", "Arg9", "Arg10", "Arg11",
        paste(tp.1, month.abb[tp.2], sep=" "))
    tmp <- data.frame(x, row.names=ROWNAMES);
  }
  else if (inherits(x, "ts") && is.null(ncol(x))) {
    COLNAMES <- switch(frequency(x),
        "Value",
        "Arg2", "Arg3",              ## Dummy arguments
        c("Q1", "Q2", "Q3", "Q4"),
        "Arg5", "Arg6",
        day.abb,
        "Arg8", "Arg9", "Arg10", "Arg11",
        month.abb)
    ROWNAMES <- seq(from=start(x)[1], to=end(x)[1])
    tmp <- data.frame(matrix(c(rep(NA, start(x)[2] - 1), x,
              rep(NA, frequency(x) - end(x)[2])),
            ncol=frequency(x), byrow=TRUE), row.names=ROWNAMES)
    names(tmp) <- COLNAMES
  }
  obj <- asciiTable$new(x = as.data.frame(tmp), include.rownames = include.rownames,
      include.colnames = include.colnames, rownames = rownames, colnames = colnames,
      format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
      caption = caption, caption.level = caption.level, width = width, frame = frame,
      grid = grid, valign = valign, header = header, footer = footer, align = align,
      col.width = col.width, style = style,
      tgroup = tgroup, n.tgroup = n.tgroup, talign = talign,
      tvalign = tvalign, tstyle = tstyle,
      bgroup = bgroup, n.bgroup = n.bgroup, balign = balign,
      bvalign = bvalign, bstyle = bstyle,
      lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign,
      lvalign = lvalign, lstyle = lstyle,
      rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
      rvalign = rvalign, rstyle = rstyle)
  return(obj)
}

##' @export
##' @method ascii zoo
ascii.zoo <- function(x, ...) {
    return(ascii(as.ts(x), ...))
}
