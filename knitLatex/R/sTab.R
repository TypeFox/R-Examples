# documentation{{{
#' Produces a latex supertabular environment
#'
#' @param firsthead header on first page of table only
#'
#' @param head header to appear at the top of every page of table
#'
#' @param tail footer on bottom of every page of table
#'
#' @param lasttail footer on last page of table only
#'
#' @param caption.table sets '\\tablecaption' option of supertabular, which
#'   allows for a default placement of the caption, see supertabular documentation
#'   for a more detailed explanation.
#'
#' @inheritParams xTab
#'
#' @examples
#' sTab(mtcars)
#' sTab(mtcars,
#'   caption.top = 'my super table',
#'   booktabs = TRUE,
#'   rows = TRUE)
#}}}
sTab <- function(x, label = NULL,
                 caption.top = NULL,
                 caption.bottom = NULL,
                 caption.table = NULL,
                 booktabs = .op('kLat.sTab.booktabs', 'kLat.booktabs', FALSE),
                 toprule = .book('kLat.toprule', booktabs, '\\toprule', '\\hline'),
                 bottomrule = .book('kLat.bottomrule', booktabs, '\\bottomrule', '\\hline'),
                 midrule = .book('kLat.midrule', booktabs, '\\midrule', '\\hline'),
                 align = .op('kLat.sTab.align', 'kLat.align', 'center'),
                 envir = getOption('kLat.sTab.envir', 'supertabular'),
                 colsep = .op('kLat.sTab.colsep', 'kLat.colsep', ''),
                 coldef = .coldef(x, colsep),
                 rowsep = .op('kLat.sTab.rowsep', 'kLat.rowsep', ''),
                 rows = .op('kLat.sTab.rows', 'kLat.rows', FALSE),
                 firsthead = NULL,
                 head = .header(x, rows),
                 tail = bottomrule,
                 lasttail = NULL){
  .pt(c(
       .printif(.printhead(toprule, firsthead, midrule), "\\tablefirsthead{%s}"),
       .printif(.printhead(toprule, head, midrule), "\\tablehead{%s}"),
       .printif(tail, "\\tabletail{%s}"),
       .printif(lasttail, "\\tablelasttail{%s}"),
       .printif(caption.top, "\\topcaption{%s}"),
       .printif(caption.bottom, "\\bottomcaption{%s}"),
       .printif(caption.table, "\\bottomcaption{%s}"),
       .printif(align, "\\begin{%s}"),
       sprintf("\\begin{%s}{%s}", envir, coldef),
       .printif(label, "\\label{%s}"),
       .body(x, rows, rowsep),
       sprintf("\\end{%s}", envir),
       .printif(align, "\\end{%s}")
       ))
}
# vim:set foldmethod=marker:
