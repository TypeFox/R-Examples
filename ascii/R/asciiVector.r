##' @export
##' @method ascii integer
ascii.integer <- function (x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
    obj <- asciiTable$new(x = x, include.rownames = include.rownames,
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
##' @method ascii numeric
ascii.numeric <- function (x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
    obj <- asciiTable$new(x = x, include.rownames = include.rownames,
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
##' @method ascii character
ascii.character <- function (x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
    obj <- asciiTable$new(x = x, include.rownames = include.rownames,
         include.colnames = include.colnames, rownames = rownames, colnames = colnames,
         format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print,
         caption = caption, caption.level = caption.level, width = width,
         frame = frame, grid = grid, valign = valign, header = header, footer = footer,
         align = align, col.width = col.width, style = style,
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
##' @method ascii factor
ascii.factor <- function (x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = FALSE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", ...){
    obj <- asciiTable$new(x = x, include.rownames = include.rownames,
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
