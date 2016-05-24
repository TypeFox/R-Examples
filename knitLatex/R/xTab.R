# documentation{{{
#' Produces a latex table
#' 
#' @param x a data.frame or matrix to form the base of the table
#'
#' @param label set the table's label, defaults to an empty string
#'
#' @param caption.top sets the caption command placing it at the top
#'   of the table
#'
#' @param caption.bottom sets the caption command placing it at the
#'   bottom of the table
#'
#' @param position sets the position of the table i.e.
#'   \\begin\{table\}['position'], defaults to 'ht'
#'
#' @param booktabs logical value, if not set will use value of
#'   kLat.(xTab|sTab|lTab).booktabs, if not set will use value of kLat.booktabs,
#'   if not set defaults to FALSE. When TRUE toprule defaults to '\\toprule',
#'   midrule to '\\midrule', and botrule to '\\bottomrule', when FALSE those
#'   values all default to '\\hline'. Has no effect when toprule, midrule, and
#'   botrule are individually set.
#'
#' @param toprule sets the value for the top rule, if not set will be
#'   determined by the value of booktabs
#'
#' @param bottomrule sets the value for the bottom rule, if not set will
#'   be determined by the value of booktabs
#'
#' @param midrule sets the value for the mid rule, if not set will
#'   be determined by the value of booktabs
#'
#' @param align set the alignment of the environment, if not set will use value
#'   of kLat.(xTab|sTab|lTab).align, if not will use value of kLat.align, if not set
#'   defaults to 'center'
#'
#' @param envir set the environment for the table, if not set will use the value
#'   of kLat.(xTab|sTab|lTab).envir, if not set defaults to 'tabular',
#'   'supertabular', and 'longtable' for xTab, sTab, and lTab respectively
#'
#' @param colsep separator to be used between columns (i.e. '|'), if not set
#'   will use the value of kLat.(xTab|sTab|lTab).colsep, if not set will use the
#'   value of kLat.colsep, if not set defaults to an empty string. If coldef is
#'   set this value is ignored and the separators must be specificed in the coldef
#'
#' @param coldef sets column definition i.e. \\begin\{tabular\}\{'align'\},
#'   if not set defaults to numeric = right, character = left 
#'
#' @param rows logical value to determine if rownames are included in table, if
#'   not set will use the value of kLat.(xTab|sTab|lTab).rows, if not set will use
#'   the value of kLat.rows, if not set defaults to FALSE, if TRUE the column name
#'   for the rownames column defaults to an empty string
#'
#' @param rowsep the separaotr to be used between rows (i.e. '\\hline'), if not
#'   set will use the value of kLat.(xTab|sTab|lTab).rowsep, if not set will use
#'   the value of kLat.rowsep, if not set defaults to an empty string
#'
#' @param head sets the value for the table header, defaults to the column
#'   names; if you set this be sure to end with '\\\\\\\\'
#'
#' @param foot sets value of the table footer, defaults to the value of
#'   botrule
#'
#' @examples
#' xTab(mtcars)
#' xTab(mtcars, label='my table', caption.top='tab:mytable', booktabs=TRUE)
#' xTab(mtcars, head='col1 & col2 & \\eta\\\\')
#}}}
xTab <- function(x, label = NULL,
                 caption.top = NULL,
                 caption.bottom = NULL,
                 position = getOption('kLat.xTab.position', 'ht'),
                 booktabs = .op('kLat.xTab.booktabs', 'kLat.booktabs', FALSE),
                 toprule = .book('kLat.toprule', booktabs, '\\toprule', '\\hline'),
                 bottomrule = .book('kLat.botrule', booktabs, '\\bottomrule', '\\hline'),
                 midrule = .book('kLat.midrule', booktabs, '\\midrule', '\\hline'),
                 align = .op('kLat.xTab.align', 'kLat.align', 'center'),
                 envir = getOption('kLat.xTab.envir', 'tabular'),
                 colsep = .op('kLat.xTab.colsep', 'kLat.colsep', ''),
                 coldef = .coldef(x, colsep),
                 rowsep = .op('kLat.xTab.rowsep', 'kLat.rowsep', ''),
                 rows = .op('kLat.xTab.rows', 'kLat.rows', FALSE),
                 head = .header(x, rows),
                 foot = bottomrule){
  .pt(c(
       paste0('\\begin{table}', .printif(position, "[%s]")),
       .printif(align, '\\begin{%s}'),
       .printif(caption.top, "\\caption{%s}"),
       sprintf('\\begin{%s}{%s}', envir, coldef),
       .printif(head, .printhead(toprule, head, midrule)),
       .body(x, rows, rowsep),
       .printif(foot, "%s"),
       sprintf('\\end{%s}', envir),
       .printif(caption.bottom, "\\caption{%s}"),
       .printif(label, "\\label{%s}"),
       .printif(align, '\\end{%s}'),
       '\\end{table}'
       ))
}
# vim:set foldmethod=marker:
