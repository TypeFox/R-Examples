##writeTags.r
##dmontaner@cipf.es
##TiddlyWikiR library

##' @name writeTags
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords replace tag
##' @seealso \code{\link{writeTiddlers}}, \code{\link{wikify}}, \code{gsub}
##' 
##' @title Replace tags within a TiddlyWiki file.
##' 
##' @description The function replaces one or more tags within a TiddlyWiki template file.
##' It is intended to insert the results of the analysis within already existing tiddlers. 
##' 
##' @details A list with the tags to be replaced in its 'names'
##' and the replacement objects in its 'elements' is required by the function.
##' 
##' Alternatively a character vector of tags and a list of replacements objects may be provided.
##' When there is just one tag to be replaced 'x' or 'rep' may not be wrapped into the list structure;
##' 'x' would be then a character vector of length one and 'rep' may be any object in R.
##' 
##' By default the TiddlyWiki template "file" will be overwritten.
##' The options "infile" and "outfile" may be used not to overwrite the template.
##'
##' A recommendation is to use something like \code{@@this_is_my_tag@@} as a 'tag'
##' to be replaced. This character string is generally unique so there Will not be
##' any problem in overwriting it. Also the @@ indicates a highlight in the TiddlyWiki
##' syntax so you will track your tags easily in your writing.
##'
##' @param x a list 
##' @param tag a character vector containing the tags to be replaced.
##' @param rep a list of the objects that will replace each tag.
##' @param file TiddlyWiki template file.
##' @param infile template file. 
##' @param outfile output file if required to be different from the input file.
##' @param verbose verbose.
##'
##' @examples
##' \dontrun{
##'   writeTags (x, tag = names (x), rep = x, file = "myTemplate.html")
##' }
##' 
##' @export

writeTags <- function (x, tag = names (x), rep = x,
                       file, infile = file, outfile = file,
                       verbose = TRUE) {
  
  ##CHECK
  if (!is.list (rep)) {
    rep <- list (rep)
  }

  if (length (tag) != length (rep)) {
    stop ('"tag" and "replacement" lengths do not match')
  }
  
  L <- length (tag)

  if (L == 0) {
    stop ('x is an empty list')
  }

  
  ##READ FILE
  if (verbose) cat ("Reading file:", infile, fill = TRUE)
  templ <- readLines (con = infile)

  
  ##REPLACEMENTS
  for (i in 1:L) {
    replace <- wikify (rep[[i]])
    replace <- paste (replace, collapse = "\n")
    templ <- gsub (pattern = tag[i], replacement = replace, templ)
  }
  
  
  ##SAVE
  if (verbose) cat ("Writing file:", outfile, fill = TRUE)
  writeLines (text = templ, con = outfile)
}
