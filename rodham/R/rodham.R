#' rodham: Hillary Rodham Clinton emails
#'
#' @section Dataset:
#' \itemize{
#' \item{emails}{dataset of 29444 emails}
#' }
#'
#' @section Functions:
#' \itemize{
#' \item \code{\link{search_emails}}
#' \item \code{\link{edges_emails}}
#' \item \code{\link{get_emails}}
#' \item \code{\link{get_xpdf}}
#' }
#'
#' @examples
#' \dontrun{
#' # get emails from internal data set
#' data("emails")
#'
#' # build graph
#' edges <- edges_emails(emails)
#' g <- igraph::graph.data.frame(edges)
#' plot(g)
#'
#' # plot communities
#' cm <- igraph::walktrap.community(g)
#' plot(cm, g)
#'
#' # get extractor to extract content from emails
#' ext <- get_xpdf()
#'
#' dir.create("./emails")
#'
#' # get emails released in august
#' aug_emails <- get_emails(release = "August", save.dir = "./emails",
#'                          extractor = ext)
#'
#' # load txt files
#' files <- list.files(aug_emails)
#' content <- lapply(1:length(files), function(x){
#'    readLines(paste0(aug_emails, "/", files[[x]]))
#' })
#' }
#'
#' @importFrom methods is
#' @importFrom utils URLencode download.file setTxtProgressBar txtProgressBar
#' @importFrom utils unzip
#'
#' @docType package
#' @name rodham
NULL
