
# Helper function used in exportDF
setOpt= function(x, defaults, colnames) {
  res= defaults
  if (!is.null(x)) {
    if (is.null(names(x)) || any(names(x) == ""))
      stop("all elements in 'x' must be named")
    if (!all(names(x) %in% colnames))
      stop(paste0("element name(s) of 'x' not in 'colnames';",
        " must be one of '",paste(colnames,collapse="', '"),"'"))
    i= match(colnames, names(x))
    res= ifelse(is.na(i), res, x[i])
  }
  return(res)
}

# Helper function to apply a list of functions to a vector of arguments of the
# same lenght
# funs: list of functions
# args: vector of arguments, each of which is passed to the corresp. element of funs
xapply= function(funs, args) {
  if (!all(sapply(funs, is.function)))
    stop("'funs' must be a list of functions")
  if (length(funs) != length(args))
    stop("'funs' and 'args' differ in length")
  return(unlist(lapply(X=1:length(funs), FUN=function(i){funs[[i]](args[i])})))
}

#' Export a Data Frame as HTML/TEX Code
#'
#' Generates code to include tabular data in a tex document or web site.
#'
#' @param x The data frame being exported.
#' @param tex Logical. Allows to switch between generation of TEX code and HTML.
#' @param colnames Displayed column names. If \code{NULL}, the original names
#'   of \code{x} are used. Otherwise it must be a named vector with element
#'   names corresponding to column names in \code{x}. It is OK to supply
#'   alternative names for selected columns only.
#' @param width Either \code{NULL} (all columns get equal width) or a named
#'   vector with element names corresponding to column names in \code{x}. If
#'   \code{tex == TRUE}, values (between 0 and 1) are needed for columns with
#'   align code 'p' only. They are interpreted as a multiplier for '\\textwidth'.
#'   If \code{tex == FALSE}, values (between 0 and 100) should be
#'   supplied for all columns of \code{x}.
#' @param align Either \code{NULL} (to use automatic alignment) or a named
#'   vector with element names corresponding to column names in \code{x}.
#'   If \code{tex == FALSE} valid alignment codes are 'left', 'right', 'center'.
#'   If \code{tex == TRUE} valid alignment codes are 'l', 'r', 'c', and 'p'. For
#'   columns with code 'p' a corresponding value of \code{width} should be set.
#'   It is OK to supply alignment codes for selected columns only.
#' @param funHead Either \code{NULL} or a list of functions whose names
#'   correspond to column names of \code{x}. The functions should have a single
#'   formal argument; the respective column names of \code{x} are used as the
#'   actual arguments. It is OK to supply functions for selected columns only
#'   (an empty function is applied to the remaining columns). See below for some
#'   typical examples.
#' @param funCell Like \code{funHead} but these functions are applied to the
#'   cells in columns rather that to the column names.
#' @param lines Logical. Switches table borders on/off.
#' @param indent Integer. Number of blanks used to indent the generated code.
#'
#' @return A character string (usually needs to be exported to a file).
#'
#' @note The functions \code{funHead} and \code{funCell} are useful to apply
#'   formatting or character replacement. For example, one could use
#'
#'   \code{function(x) {paste0("\\\\bold{",toupper(x),"}")}}
#'
#'   to generate bold, uppercase column names in a TEX table.
#'
#' @seealso The \code{xtable} packages provides similar functionality with
#'   more sophisticated options. Consider the 'pandoc' software do convert
#'   documents from one markup language to another one. Finally, consider the
#'   latex package 'datatools' for direct inclusion of delimited text files
#'   (e.g. produced by \code{write.table}) in tex documents.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # Create example table
#' df= data.frame(stringsAsFactors=FALSE, name= c("growth", "dead"),
#'   unit= c("1/d","1/d"), expression= c("r * N * (1 - N/K)"," d * N"))
#'
#' # Export as TEX: header in bold, 1st colum in italics, last column as math
#' tex= exportDF(df, tex=TRUE,
#'   colnames=c(expression="process rate expression"),
#'   width=c(expression=0.5),
#'   align=c(expression="p"),
#'   funHead=setNames(replicate(ncol(df),
#'     function(x){paste0("\\textbf{",x,"}")}),names(df)),
#'   funCell=c(name=function(x){paste0("\\textit{",x,"}")},
#'     expression=function(x){paste0("$",x,"$")})
#' )
#' cat(tex,"\n")
#'
#' # Export as HTML: non-standard colors are used for all columns
#' tf= tempfile(fileext=".html")
#' write(x= exportDF(df, tex=FALSE,
#'   funHead=setNames(replicate(ncol(df),
#'     function(x){paste0("<font color='red'>",x,"</font>")}),names(df)),
#'   funCell=setNames(replicate(ncol(df),
#'     function(x){paste0("<font color='blue'>",x,"</font>")}),names(df))
#' ), file=tf)
#' \dontrun{
#'   browseURL(tf)
#'   file.remove(tf)
#' }

exportDF= function(x,
  tex=FALSE,
  colnames=NULL,
  width= NULL,
  align= NULL,
  funHead= NULL,
  funCell= NULL,
  lines=TRUE,
  indent=2
) {
  indent= ifelse(indent <= 0, "", paste0(rep(" ",indent),collapse=""))
  # Check input
  if (is.matrix(x))
    x= as.data.frame(x, stringsAsFactors=FALSE)
  if (!is.data.frame(x))
    stop("'x' must be  data frame")
  # Set options
  left= ifelse(tex, "l", "left")
  right= ifelse(tex, "r", "right")
  none= function(x) {x}
  w= ifelse(tex, 1/ncol(x), floor(100/ncol(x)))
  colnames= setOpt(colnames, names(x), names(x))
  width=    setOpt(width, rep(w, ncol(x)), names(x))
  align=    setOpt(align, ifelse(unlist(lapply(x, FUN=is.numeric)),right,left), names(x))
  funHead=  setOpt(funHead, replicate(n=ncol(x), none), names(x))
  funCell=  setOpt(funCell, replicate(n=ncol(x), none), names(x))
  # Assemble code
  out=''

  # tex
  if (tex) {
    i= which(align == "p")
    if (length(i) > 0) {
      align[i]= paste0(align[i],"{",width[i],"\\textwidth}")
    }
    out= paste0(out,indent,'\\begin{tabular}{',paste(align,collapse=""),
      '}',ifelse(lines, '\\hline', ''),'\n')
    out= paste0(out,indent,'  ',
      paste0(paste0('',xapply(funHead,colnames),''),collapse=' & '),' \\\\',
      ifelse(lines, ' \\hline', ''),'\n')
    for (i in 1:nrow(x)) {
      out= paste0(out,indent,'  ',
        paste0(paste0('',xapply(funCell,unlist(x[i,])),''),collapse=' & '),' \\\\',
        ifelse(lines && (i == nrow(x)), ' \\hline', ''),'\n')
    }
    out= paste0(out,indent,'\\end{tabular}\n')

  # html
  } else {
    out= paste0(out,indent,'<table border=',ifelse(lines,1,0),'>\n')
    for (i in 1:length(width)) {
      out= paste0(out,indent,'  <col width="',width[i],'%">\n')
    }
    out= paste0(out,indent,'  <tr>',
      paste0(paste0('<th style="text-align:',align,'"> ',xapply(funHead,
      colnames),' </th>'),collapse=''),' </tr>\n')
    for (i in 1:nrow(x)) {
      out= paste0(out,indent,'  <tr>',
        paste0(paste0('<td style="text-align:',align,'"> ',xapply(funCell,
        unlist(x[i,])),' </td>'),collapse=''),' </tr>\n')
    }
    out= paste0(out,indent,'</table>\n')
  }

  return(out)
}

