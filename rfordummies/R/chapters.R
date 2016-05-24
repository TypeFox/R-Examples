### Code to automatically generate ch() functions

# Read cleaned chapter code from /inst/scripts/2-clean
# 
# @param ch Chapter number
# @keywords Internal
.readExampleFile <- function(ch){
  text <- system.file(paste0("scripts/2-clean/", ch), package = "rfordummies")
  readLines(text)
}


# Generate code and example for single chapter.
# 
# @param ch Chapter number
# @keywords Internal
# @rdname genChapter
.genChapter <- function(ch){
  chapter <- "#' Print examples of chapter <<x>> of 'R for Dummies'.\n#'\n#' To print a listing of all examples of a chapter, use \\code{ch<<x>>()}. To run all the examples of \\code{ch<<x>>()}, use \\code{example(ch<<x>>)}.\n#' @export\n#' @rdname ch<<xx>>\n#' @family Chapters\n#' @seealso \\code{\\link{toc}}\n#' @example inst/scripts/2-clean/ch<<xx>>.R\nch<<xx>> <- function(){\n  text <- .readExampleFile(\"ch<<xx>>.R\")\n  cat(text, sep=\"\n\")
  invisible(text)\n}\n"
chapterExtra <- "#' @export\n#' @aliases ch<<xx>>\n#' @rdname ch<<xx>>\nch<<x>> <- ch<<xx>>"
  if(!is.numeric(ch)) stop ("ch should be numeric")
  chapter <- gsub("<<x>>", sprintf("%s", ch), chapter)
  out1 <- gsub("<<xx>>", sprintf("%02d", ch), chapter)
  
  if(ch < 10){
    chapterExtra <- gsub("<<x>>", sprintf("%s", ch), chapterExtra)
    out2 <- gsub("<<xx>>", sprintf("%02d", ch), chapterExtra)
  } else {
    out2 <- ""
  }
  sprintf("%s\n%s\n\n", out1, out2)
}

# Generate code and examples for all chapters.
# 
# @param chapters Chapter numbers
# @param path Path to package
# @param file Name of file containing ch() functions
# @keywords Internal
.generateChapters <- function(chapters=1:20, path, file="R/chapters-auto.R"){
  sink(file.path(path, file))
  on.exit(sink())  
  for(i in chapters)cat(.genChapter(i))
}

# .generateChapters(path="rfordummies")

