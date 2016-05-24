# documentation{{{
#' Produces a latex longtable
#'
#' @param firsthead header on first page of table only; defaults to header; if you
#'   set this, you are responsible for setting any \\hline, \\toprule, or \\midrule
#'   lines
#'
#' @param lastfoot footer on last page of table only
#'
#' @param caption the caption for the table, unlinke xTab and sTab, there is no
#'   caption.top or caption.bottom option in longtable
#'
#' @param caption.firsthead,caption.head,caption.foot,caption.lastfoot paces the caption in the
#'   firsthead, head, foot, or lastfoot respectively. It is important not to set a
#'   caption in an otherwise NULL section (although an empty string is acceptable)
#'   or strange bugs can occur. It is accaptable if the section was set by default
#'   as in head and foot. Consult the longtable documentation for a more detailed
#'   explanation of these options.
#'
#' @inheritParams xTab
#'
#' @examples
#' lTab(mtcars)
#}}}
lTab <- function(x, label = NULL,
                 caption.firsthead = NULL,
                 caption.head = NULL,
                 caption.foot = NULL,
                 caption.lastfoot = NULL,
                 booktabs = .op('kLat.lTab.booktabs', 'kLat.booktabs', FALSE),
                 toprule = .book('kLat.toprule', booktabs, '\\toprule', '\\hline'),
                 bottomrule = .book('kLat.bottomrule', booktabs, '\\bottomrule', '\\hline'),
                 midrule = .book('kLat.midrule', booktabs, '\\midrule', '\\hline'),
                 align = .op('kLat.lTab.align', 'kLat.align', 'center'),
                 envir = getOption('kLat.lTab.envir', 'longtable'),
                 colsep = .op('kLat.lTab.colsep', 'kLat.colsep', ''),
                 coldef = .coldef(x, colsep),
                 rowsep = .op('kLat.lTab.rowsep', 'kLat.rowsep', ''),
                 rows = .op('kLat.lTab.rows', 'kLat.rows', FALSE),
                 head = .header(x, rows),
                 firsthead = NULL,
                 foot = bottomrule,
                 lastfoot = NULL){
  .pt(c(
       .printif(align, "\\begin{%s}"),
       sprintf('\\begin{%s}{%s}', envir, coldef),
       .printif(caption.firsthead, "\\caption{%s}\\\\"),
       .printif(.printhead(toprule, firsthead, midrule), "%s\n\\endfirsthead"),
       .printif(caption.head, "\\caption{%s}\\\\"),
       .printif(.printhead(toprule, head, midrule), "%s\n\\endhead"),
       .printif(caption.foot, "\\caption{%s}\\\\"),
       .printif(foot, "%s\n\\endfoot"),
       .printif(caption.lastfoot, "\\caption{%s}\\\\"),
       .printif(lastfoot, "%s\n\\endlastfoot"),
       .body(x, rows, rowsep),
       .printif(label, "\\label{%s}"),
       sprintf('\\end{%s}', envir),
       .printif(align, '\\end{%s}')
       ))
}
# vim:set foldmethod=marker:
