setEurostatTOC <-
function() {
   if (!exists(".eurostatTOC", envir = .SmarterPolandEnv)) {
   .eurostatTOC <-  read.table("http://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?sort=1&file=table_of_contents_en.txt",  sep="\t", header=T,  quote="\"", fill = TRUE, comment.char="")
    assign(".eurostatTOC", .eurostatTOC, envir = .SmarterPolandEnv)
  }
  invisible(0)
}
