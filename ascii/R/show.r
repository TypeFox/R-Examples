##' ascii table generator
##'
##' @author David Hajage
##' @export
asciiTable <- setRefClass("asciiTable",
                          fields = c("x",
                            "include.rownames",
                            "include.colnames",
                            "rownames",
                            "colnames",
                            "format",
                             "digits",
                            "decimal.mark",
                            "na.print",
                            "caption",
                            "caption.level",
                            "width",
                            "frame",
                            "grid",
                            "valign",
                            "header",
                            "footer",
                            "align",
                            "col.width",
                            "style",
                            "tgroup",
                            "n.tgroup",
                            "talign",
                            "tvalign",
                            "tstyle",
                            "bgroup",
                            "n.bgroup",
                            "balign",
                            "bvalign",
                            "bstyle",
                            "lgroup",
                            "n.lgroup",
                            "lalign",
                            "lvalign",
                            "lstyle",
                            "rgroup",
                            "n.rgroup",
                            "ralign",
                            "rvalign",
                            "rstyle"),
                          methods = list(
                            show.asciidoc = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                              'print a table with asciidoc markup'
                              show.asciidoc.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                            },

                            show.rest = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                               'print a table with restructuredText markup'
                               show.rest.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                             },

                            show.org = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                              'print a table with org-mode markup'
                              show.org.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                            },

                            show.t2t = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                              'print a table with txt2tags markup'
                              show.t2t.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                            },

                            show.textile = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                              'print a table with textile markup'
                              show.textile.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                            },

                            show.pandoc = function(x = .self$x, include.rownames = .self$include.rownames, include.colnames = .self$include.colnames, rownames = .self$rownames, colnames = .self$colnames, format = .self$format, digits = .self$digits, decimal.mark = .self$decimal.mark, na.print = .self$na.print, caption = .self$caption, caption.level = .self$caption.level, width = .self$width, frame = .self$frame, grid = .self$grid, valign = .self$valign, header = .self$header, footer = .self$footer, align = .self$align, col.width = .self$col.width, style = .self$style, lgroup = .self$lgroup, n.lgroup = .self$n.lgroup, lalign = .self$lalign, lvalign = .self$lvalign, lstyle = .self$lstyle, rgroup = .self$rgroup, n.rgroup = .self$n.rgroup, ralign = .self$ralign, rvalign = .self$rvalign, rstyle = .self$rstyle, tgroup = .self$tgroup, n.tgroup = .self$n.tgroup, talign = .self$talign, tvalign = .self$tvalign, tstyle = .self$tstyle, bgroup = .self$bgroup, n.bgroup = .self$n.bgroup, balign = .self$balign, bvalign = .self$bvalign, bstyle = .self$bstyle) {
                              'print a table with pandoc markup'
                              show.pandoc.table(x = x, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimalmark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle)
                            }
                            )
                          )

##' ascii list generator
##'
##' @author David Hajage
##' @export
asciiList <- setRefClass("asciiList",
                         fields = c("x",
                           "caption",
                           "caption.level",
                           "list.type"),
                         methods = list(
                           show.asciidoc = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with asciidoc markup'
                             show.asciidoc.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           },

                           show.rest = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with rest markup'
                             show.rest.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           },

                           show.org = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with org markup'
                             show.org.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           },

                           show.t2t = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with t2t markup'
                             show.t2t.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           },

                           show.textile = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with textile markup'
                             show.textile.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           },

                           show.pandoc = function(x = .self$x, caption = .self$caption, caption.level = .self$caption.level, list.type = .self$list.type) {
                             'print a list with pandoc markup'
                             show.pandoc.list(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
                           }
                           )
                         )

##' ascii mixed generator
##'
##' @author David Hajage
##' @export
asciiMixed <- setRefClass("asciiMixed",
                          fields =  list(args = "list"),

                          methods = list(
                            show.asciidoc = function() {
                             'print everything with asciidoc markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "asciidoc")
                                if (i != length(args)) cat("\n")
                              }
                            },

                            show.rest = function() {
                             'print everything with rest markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "rest")
                                if (i != length(args)) cat("\n")
                              }
                            },

                            show.org = function() {
                             'print everything with org markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "org")
                                if (i != length(args)) cat("\n")
                              }
                            },

                            show.t2t = function() {
                             'print everything with t2t markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "t2t")
                                if (i != length(args)) cat("\n")
                              }
                            },

                            show.textile = function() {
                             'print everything with textile markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "textile")
                                if (i != length(args)) cat("\n")
                              }
                            },

                            show.pandoc = function() {
                             'print everything with pandoc markup'
                              for (i in seq_along(args)) {
                                if (is.null(args[[i]])) next
                                print(args[[i]], type = "pandoc")
                                if (i != length(args)) cat("\n")
                              }
                            }
                            )
                          )
