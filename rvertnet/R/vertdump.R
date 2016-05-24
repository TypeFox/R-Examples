#' Use Vertnet taxon specific dump from KNB
#'
#' @name dump
#' @param path (character) Path to a sqlite file on your machine.
#' @param group (character) One of mammals, reptiles, amphibians, fishes, or birds
#' @param table (character) sqlite table name, you can use anything you like, but
#' code defaults to the values in the \code{group} parameter
#' @param x An object of class \code{src_sqlite}
#'
#' @details \code{dump_init} creates a \code{src_sqlite} class that you can use to
#' query the data either using SQL syntax or dplyr's R syntax. \code{dump_tbl} is
#' just a wrapper around \code{\link[dplyr]{tbl}} to create a \code{tbl} class
#' object that you can use to feed directly into dplyr's verbs, like
#' \code{\link[dplyr]{select}} and \code{\link[dplyr]{filter}}
#'
#' @references
#' \url{http://blog.vertnet.org/post/115875718156/the-data-one-thing-about-vertnet-and-big-data}
#'
#' @examples \dontrun{
#' # You first need to create your SQLite databases, e.g, for amphibians:
#' ## In the terminal
#'# wget 
#'# https://knb.ecoinformatics.org/knb/d1/mn/v1/object/urn:uuid:afc58110-b9c1-4cf7-b46c-837bdc930a21
#' # mv urn\:uuid\:afc58110-b9c1-4cf7-b46c-837bdc930a21 vertnet_amphib.gz
#' # gunzip vertnet_amphib.gz
#' # sqlite3 amphibians.sqlite
#'
#' ## In SQLite
#' # sqlite> .separator ','
#' # sqlite> .import vertnet_amphib amphibians
#' 
#' # After you have a SQLite database, do
#' # library("dplyr")
#' # x <- dump_init(path = "~/github/sac/vertnetdumps/amphibians.sqlite")
#'
#' # use SQL syntax
#' # tbl(x, sql("SELECT scientificname,title FROM amphibians LIMIT 10"))
#'
#' # use R syntax
#' # tab <- x %>% dump_tbl()
#' # tab %>%
#' #  filter(year > 2010) %>%
#' #  select(scientificname, title)
#' }

#' @export
#' @rdname dump
dump_init <- function(path, group = "amphibians", table = NULL) {
  checkfourpkg("RSQLite")
  group <- match.arg(group, c("mammals", "reptiles", "amphibians", "fishes", "birds"))
  if (is.null(table)) table <- group
  x <- dplyr::src_sqlite(path = path)
  structure(x, table = table)
}

#' @export
#' @rdname dump
dump_tbl <- function(x) {
  stopifnot(is(x, "src_sqlite"))
  dplyr::tbl(x, attr(x, "table"))
}

#' @export
#' @rdname dump
dump_links <- function() {
  list(
    mammals = list(
      table = "mammals",
      doi = "10.5063/F1GQ6VPM",
      data = "1d09e64b-d25a-46e7-bdc1-0a91fb7bf8bb",
      view = "https://knb.ecoinformatics.org/#view/doi:10.5063/F1GQ6VPM",
      eml = ""
    ),
    reptiles = list(
      table = "reptiles",
      doi = "10.5063/F10P0WX6",
      data = "14a66a88-592e-4459-ae8a-b114165106c3",
      view = "https://knb.ecoinformatics.org/#view/doi:10.5063/F10P0WX6",
      eml = ""
    ),
    amphibians = list(
      table = "amphibians",
      doi = "10.5063/F1VX0DF9",
      data = "afc58110-b9c1-4cf7-b46c-837bdc930a21",
      view = "https://knb.ecoinformatics.org/#view/doi:10.5063/F1VX0DF9",
      eml = ""
    ),
    fishes = list(
      table = "fishes",
      doi = "10.5063/F1R49NQB",
      data = paste0(dump_base("data"), "74f312ce-de6f-4cae-b5f7-26d4d0ffc781"),
      view = paste0(dump_base("view"), "10.5063/F1R49NQB"),
      eml = paste0(dump_base("eml"), "10.5063/F1R49NQB")
    ),
    birds = list(
      table = "birds",
      doi = "10.5063/F1MG7MDB",
      data = paste0(dump_base("data"), "cca85e21-6647-491c-aa57-89a7c552e564"),
      view = paste0(dump_base("view"), "10.5063/F1MG7MDB"),
      eml = paste0(dump_base("eml"), "10.5063/F1MG7MDB")
    )
  )
}

# helpers ------------------
dump_base <- function(x) {
  switch(x,
         data = "https://knb.ecoinformatics.org/knb/d1/mn/v1/object/urn:uuid:",
         view = "https://knb.ecoinformatics.org/#view/doi:",
         eml = "https://knb.ecoinformatics.org/knb/d1/mn/v1/object/doi:"
  )
}
