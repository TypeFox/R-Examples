##' @export
##' @method ascii anova
ascii.anova <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
  x <- as.data.frame(x, check.names = FALSE)
  obj <- asciiTable$new(x = x, include.rownames = include.rownames,
                        include.colnames = include.colnames,
                        rownames = rownames, colnames = colnames, format = format,
                        digits = digits, decimal.mark = decimal.mark, na.print = na.print,
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
##' @method ascii aov
ascii.aov <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
  ascii.anova(unclass(summary(x))[[1]], include.rownames = include.rownames,
              include.colnames = include.colnames, rownames = rownames, colnames = colnames,
              format = format,digits = digits, decimal.mark = decimal.mark, na.print = na.print,
              caption = caption, caption.level = caption.level, width = width, frame = frame,
              grid = grid, valign = valign, header = header, footer = footer, align = align,
              col.width = col.width, style = style,
              tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
              tstyle = tstyle,
              bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
              bstyle = bstyle,
              lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
              lstyle = lstyle,
              rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
              rvalign = rvalign, rstyle = rstyle)
}

##' @export
##' @method ascii summary.aov
ascii.summary.aov <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
  ascii.anova(unclass(x)[[1]], include.rownames = include.rownames,
         rownames = rownames, colnames = colnames,
         include.colnames = include.colnames, format = format,
         digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, caption.level = caption.level, width = width, frame = frame,
         grid = grid,valign = valign, header = header, footer = footer, align = align,
         col.width = col.width, style = style,
         tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
         tstyle = tstyle,
         bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
         bstyle = bstyle,
         lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
         lstyle = lstyle,
         rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
         rvalign = rvalign, rstyle = rstyle)
}

##' @export
##' @method ascii aovlist
ascii.aovlist <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
  y <- summary(x)
  n <- length(y)
  if (n == 1) ascii.anova(unclass(y[[1]][[1]]), include.rownames = include.rownames,
         include.colnames = include.colnames, rownames = rownames, colnames = colnames,
         format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, caption.level = caption.level, width = width, frame = frame,
         grid = grid,valign = valign, header = header, footer = footer, align = align,
         col.width = col.width, style = style,
         tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
         tstyle = tstyle,
         bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
         bstyle = bstyle,
         lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
         lstyle = lstyle,
         rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
         rvalign = rvalign, rstyle = rstyle)
  else {
    z <- y[[1]][[1]]
    for (i in 2:n) z <- rbind(z, y[[i]][[1]])
    ascii.anova(z, include.rownames = include.rownames,
         include.colnames = include.colnames, rownames = NULL, colnames = NULL, format = format,
         digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, width = width, frame = frame, grid = grid,
         valign = valign, header = header, footer = footer, align = align,
         col.width = col.width, style = style,
         tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
         tstyle = tstyle,
         bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
         bstyle = bstyle,
         lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
         lstyle = lstyle,
         rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
         rvalign = rvalign, rstyle = rstyle)
  }
}

##' @export
##' @method ascii summary.aovlist
ascii.summary.aovlist <- function (x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = "middle", header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
  n <- length(x)
  if (n == 1) ascii.anova(unclass(x[[1]][[1]]), include.rownames = include.rownames,
         include.colnames = include.colnames, rownames = rownames, colnames = colnames,
         format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, caption.level = caption.level, width = width, frame = frame,
         grid = grid, valign = valign, header = header, footer = footer, align = align,
         col.width = col.width, style = style,
         tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
         tstyle = tstyle,
         bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
         bstyle = bstyle,
         lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
         lstyle = lstyle,
         rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
         rvalign = rvalign, rstyle = rstyle)
  else {
    z <- x[[1]][[1]]
    for (i in 2:n) z <- rbind(z, x[[i]][[1]])
    ascii.anova(z, include.rownames = include.rownames,
         include.colnames = include.colnames, rownames = rownames, colnames = colnames,
         format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, caption.level = caption.level, width = width, frame = frame,
         grid = grid, valign = valign, header = header, footer = footer, align = align,
         col.width = col.width, style = style,
         tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign,
         tstyle = tstyle,
         bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign,
         bstyle = bstyle,
         lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign,
         lstyle = lstyle,
         rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign,
         rvalign = rvalign, rstyle = rstyle)
  }
}
