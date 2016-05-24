#' Download data from Pangaea.
#'
#' Grabs data as a dataframe or list of dataframes from a Pangaea data repository URI; see:
#' \url{https://www.pangaea.de/}.
#'
#' @export
#' @param doi DOI of Pangaeae single dataset, or of a collection of datasets. Expects
#' either just a DOI of the form \code{10.1594/PANGAEA.746398}, or with the URL part
#' in front, like \code{https://doi.pangaea.de/10.1594/PANGAEA.746398}
#' @param overwrite (logical) Ovewrite a file if one is found with the same name
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @param prompt (logical) Prompt before clearing all files in cache? No prompt used when DOIs
#' assed in. Default: \code{TRUE}
#' @return One or more items of class pangaea, each with a citation object, metadata object,
#' and data object. Each data object is printed as a \code{tbl_df} object, but the
#' actual object is simply a \code{data.frame}.
#' @author Naupaka Zimmerman
#' @references \url{https://www.pangaea.de}
#' @details Data files are stored in an operating system appropriate location.
#' Run \code{rappdirs::user_cache_dir("pangaear")} to get the storage location
#' on your machine.
#' @examples \dontrun{
#' # a single file
#' res <- pg_data(doi='10.1594/PANGAEA.807580')
#' res
#' res[[1]]$doi
#' res[[1]]$citation
#' res[[1]]$meta
#' res[[1]]$data
#'
#' # another single file
#' pg_data(doi='10.1594/PANGAEA.807584')
#'
#' # Many files
#' res <- pg_data(doi='10.1594/PANGAEA.761032')
#' res[[1]]
#' res[[2]]
#'
#' # Manipulating the cache
#' ## list files in the cache
#' pg_cache_list()
#'
#' ## clear all data
#' # pg_cache_clear()
#' pg_cache_list()
#'
#' ## clear a single dataset by DOI
#' pg_data(doi='10.1594/PANGAEA.812093')
#' pg_cache_list()
#' pg_cache_clear(doi='10.1594/PANGAEA.812093')
#' pg_cache_list()
#'
#' ## clear more than 1 dataset by DOI
#' lapply(c('10.1594/PANGAEA.746398','10.1594/PANGAEA.746400'), pg_data)
#' pg_cache_list()
#' pg_cache_clear(doi=c('10.1594/PANGAEA.746398','10.1594/PANGAEA.746400'))
#' pg_cache_list()
#' }

pg_data <- function(doi, overwrite = TRUE, ...) {
  dois <- check_many(doi)
  invisible(lapply(dois, function(x) {
    if ( !is_pangaea(env$path, x) ) {
      pang_GET(url = paste0(base(), x), doi = x, overwrite, ...)
    }
  }))
  out <- process_pg(dois)
  lapply(out, structure, class = "pangaea")
}

#' @export
print.pangaea <- function(x, ...) {
  cat(sprintf("<Pangaea data> %s", x$doi), sep = "\n")
  print(x$data)
}

print.meta <- function(x, ...){
  cat(x$meta, sep = "\n")
}

print.citation <- function(x, ...){
  cat(x$citation, sep = "\n")
}

pang_GET <- function(url, doi, overwrite, ...){
  dir.create(env$path, showWarnings = FALSE, recursive = TRUE)
  fname <- rdoi(doi)
  res <- httr::GET(url,
             query = list(format = "textfile", charset = "UTF-8"),
             httr::config(followlocation = TRUE),
             httr::write_disk(file.path(env$path, fname), overwrite), ...)
  httr::stop_for_status(res)
}

process_pg <- function(x){
  lapply(x, function(m){
    list(doi = m,
         citation = pg_citation(m),
         meta = get_meta(file.path(env$path, rdoi(m))),
         #data = read_csv(file.path(env$path, rdoi(m)))
         data = as_data_frame(read_csv(file.path(env$path, rdoi(m))))
    )
  })
}

pg_citation <- function(x){
  structure(list(
    citation = sprintf('See https://doi.pangaea.de/%s for the citation', x)),
            class = "citation")
}

is_pangaea <- function(x, doi){
  if ( identical(list.files(x), character(0)) ) { FALSE } else {
    if ( any(rdoi(doi) %in% list.files(x)) ) TRUE else FALSE
  }
}

get_meta <- function(x){
  lns <- readLines(x, n = 300)
  ln_no <- grep("\\*/", lns) - 1
  use <- lns[2:ln_no]
  structure(list(meta = use), class = "meta")
}

rdoi <- function(x) paste0(gsub("/|\\.", "_", x), ".txt")

check_many <- function(x){
  res <- httr::GET(fix_doi(x))
  if (!grepl("name=\"dslist\"", content(res, "text", encoding = "UTF-8"))) {
    x
  } else {
    d <- gregexpr("<div class=\"MetaHeaderItem\"><a rel=\"follow\" href=\"(https://doi.pangaea.de/.*?)\">", res)
    d <- unlist(regmatches(content(res, "text", encoding = "UTF-8"), d))
    split_d <- strsplit(d, split = "\"")
    vapply(split_d, function(x) sub("https://doi.pangaea.de/", "", x[grepl("doi",x)]), "")
  }
}

fix_doi <- function(x) {
  if (grepl("https?://doi.pangaea.de/?", x)) {
    x
  } else {
    # make sure doi is cleaned up before making a url
    if (!grepl("^10.1594", x)) stop(x, " not of right form, expecting a DOI, see pg_data help file", call. = FALSE)
    paste0(base(), x)
  }
}
