#' Plot PLOS article-level metrics data.
#'
#' This can be used in conjuction with the function \code{\link{alm_signposts}}.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' 
#' @param input A data.frame from a search from the \code{\link{alm_signposts}} function
#' @details Note that DOIs are the unit of replication of each study. When plotting,
#' if the prefix is common among all DOIs, then just the end of the DOI, the numeric
#' part is printed to make plots less ugly.
#' @seealso \code{\link{alm_signposts}}
#' @references See a tutorial/vignette for alm at
#' \url{http://ropensci.org/tutorials/alm_tutorial.html}
#' @examples \dontrun{
#' # Plot data from a single identifier gives a bar chart
#' dat <- alm_signposts(doi="10.1371/journal.pone.0029797")
#' plot_signposts(dat)
#'
#' # Plot data from many identifiers gives a line chart
#' dois <- c('10.1371/journal.pone.0001543','10.1371/journal.pone.0040117',
#'    '10.1371/journal.pone.0029797','10.1371/journal.pone.0039395')
#' dat <- alm_signposts(doi=dois)
#' plot_signposts(dat)
#' 
#' # software lagotto instance
#' urls <- c("https://github.com/najoshi/sickle","https://github.com/lh3/wgsim",
#'    "https://github.com/jstjohn/SeqPrep")
#' dat <- alm_signposts(url = urls, api_url = "http://software.lagotto.io/api/v5/articles")
#' plot_signposts(dat)
#' 
#' # scopus ids
#' ids <- c(68049122102, 14044251458, 48349097292, 28444460441)
#' dat <- alm_signposts(scp = ids)
#' plot_signposts(dat)
#' }

plot_signposts <- function(input){
  if(!is.data.frame(input))
    stop("object passed to input must be a data.frame")

  if(nrow(input)==1){
    df <- melt(input)
    ggplot(df, aes(variable, value)) +
      theme_grey(base_size=16) +
      geom_bar(stat="identity") +
      labs(x="",y="")
  } else {
    input <- stripdois(input)
    df <- melt(input)
    ggplot(df, aes(id, value)) +
      theme_grey(base_size=16) +
      geom_bar(stat="identity") +
      facet_wrap(~variable, scales="free") +
      labs(x="",y="")
  }
}

stripdois <- function(x){
  if(strsplit(x$id[1], '/')[[1]][1] == "url"){
    x$id <- sapply(strsplit(as.character(x$id), "/"), "[[", 6)
    x
  } else {
    x$id <- sapply(strsplit(as.character(x$id), "\\."), "[[", 4)
    x
  }
}
