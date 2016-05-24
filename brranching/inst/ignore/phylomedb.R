#' @title PhylomeDB
#'
#' @description Fetch phylome trees from PhylomeDB
#'
#' @export
#' @param seqid An id
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @return Newick formatted tree or nexml text.
#' @examples \dontrun{
#' # Fetch by seqid
#' id <- "Phy004LGJW_CROPO"
#' tree <- phylomedb(seqid = id)
#' plot(tree, no.margin=TRUE)
#' }

phylomedb <- function(seqid, ...) {
  args <- list(q = "search_tree", seqid = seqid)
  gzpath <- tempfile(fileext = ".tar.gz")
  tt <- GET(phydb_base, query = args, config(followlocation=1))
  stop_for_status(tt)
  out <- content(tt, as = "text")
}

phydb_base <- "http://phylomedb.org"

tar_url <- function(x) {
  txt <- content(x, 'text')
  grep("download data\\.tar\\.gz", txt)
}

do_tar <- function(x) {
  tt <- GET(url, query = args, write_disk())
}
