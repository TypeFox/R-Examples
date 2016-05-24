#' Download and view the TrueFX(tm) Market Data Web API Developer Guide
#' @author Garrett See
#' @note the idea for this function came from the \pkg{IBrokers} function 
#' \code{IBrokersRef} written by Jeff Ryan. Thanks Jeff!
#' @param show What to show. Either \dQuote{pdf} or \dQuote{webpage} 
#'   (\dQuote{web} and \dQuote{url} will also be recognized as 
#'   \dQuote{webpage}).
#'   Alternatively, can be numeric: 1 for \dQuote{pdf}, 2 for \dQuote{webpage})
#' @return called for side-effect.  Opens the TrueFX(tm) Market Data Web API 
#'   Developer Guide in a pdf-viewer or web browser.
#' @references 
#' \url{http://www.truefx.com/dev/data/TrueFX_MarketDataWebAPI_DeveloperGuide.pdf}
#' @examples
#' \dontrun{
#' TrueFXRef()
#' TrueFXRef("web")
#' }
#' @export
TrueFXRef <- function (show = c("pdf", "webpage")) 
{
  if (is.numeric(show)) 
    show <- c("pdf", "webpage")[show]
  switch(show[[1]], pdf = {
    tmp <- tempfile()
    download.file(paste0("http://www.truefx.com/dev/data/TrueFX_",
                         "MarketDataWebAPI_DeveloperGuide.pdf"), destfile = tmp)
    if (.Platform$OS.type == "windows") {
      shell.exec(tmp)
    } else system(paste(shQuote(getOption("pdfviewer")), 
                        shQuote(tmp)), wait = FALSE)
  }, webpage = , web = , url = {
    browseURL(paste0("http://www.truefx.com/dev/data/TrueFX_",
                     "MarketDataWebAPI_DeveloperGuide.pdf"))
  })
}
