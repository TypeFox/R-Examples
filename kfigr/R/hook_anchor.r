#' Anchor Chunk Hook Extension for package:knitr
#' 
#' \pkg{knitr} hook functions are called when the corresponding chunk options 
#'  are not \code{NULL} to do additional jobs beside the R code in chunks. 
#'  \pkg{kfigr} provides the hook "anchor" which adds an HTML anchor tag 
#'  immediately above a code chunk.
#'
#' @details the function hook_anchor is set as a hook in \pkg{knitr} when 
#'   \pkg{kfigr} is attached (and removed when \pkg{kfigr} is detached). It 
#'   writes an HMTL anchor tag directly above a code chunk in the form 
#'   \code{<a name="chunk-name"></a>} where \code{chunk-name} is the chunk 
#'   label contained in \code{options$label}.
#'
#' @references \url{http://yihui.name/knitr/hooks#chunk_hooks}
#'
#' @seealso \code{\link{figr}}, \code{\link{anchors}}
#'
#' @param before,options,envir see references
#'
#' @examples
#' \dontrun{
#' require(knitr)
#' knit_hooks$set(anchor = hook_anchor)
#' # then in code chunks, use e.g. the option anchor = "figure"
#' }
#'
#' @export
hook_anchor <- function(before, options, envir){
  if (before){  
    invisible(index(options$label, options$anchor))
    paste('<a name="', options$label, '"></a>', sep='')
  }
}
