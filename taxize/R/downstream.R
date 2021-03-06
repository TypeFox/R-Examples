#' Retrieve the downstream taxa for a given taxon name or ID.
#'
#' This function uses a while loop to continually collect children taxa down to the
#' taxonomic rank that you specify in the \code{downto} parameter. You can get data
#' from ITIS (itis) or Catalogue of Life (col). There is no method exposed by itis
#' or col for getting taxa at a specific taxonomic rank, so we do it ourselves inside
#' the function.
#'
#' @param x Vector of taxa names (character) or IDs (character or numeric) to query.
#' @param db character; database to query. One or both of \code{itis}, \code{col}, or
#' \code{gbif}. Note that each taxonomic data source has their own identifiers, so that if you
#' provide the wrong \code{db} value for the identifier you could get a result,
#' but it will likely be wrong (not what you were expecting).
#' @param downto What taxonomic rank to go down to. One of: 'superkingdom','kingdom',
#' subkingdom','infrakingdom','phylum','division','subphylum','subdivision','infradivision',
#' 'superclass','class','subclass','infraclass','superorder','order','suborder',
#' 'infraorder','superfamily','family','subfamily','tribe','subtribe','genus','subgenus',
#' 'section','subsection','species','subspecies','variety','form','subvariety','race',
#' 'stirp','morph','aberration','subform', or 'unspecified'
#' @param intermediate (logical) If TRUE, return a list of length two with target
#' taxon rank names, with additional list of data.frame's of intermediate
#' taxonomic groups. Default: FALSE
#' @param rows (numeric) Any number from 1 to inifity. If the default NA, all rows are
#' considered. Note that this parameter is ignored if you pass in a taxonomic id of any of the
#' acceptable classes: tsn, colid.
#' @param ... Further args passed on to \code{itis_downstream} or \code{col_downstream}
#'
#' @return A named list of data.frames with the downstream names of every supplied taxa.
#' You get an NA if there was no match in the database.
#'
#' @export
#' @examples \dontrun{
#' # Plug in taxon IDs
#' ## col Ids have to be character, as they are alphanumeric IDs
#' downstream("015be25f6b061ba517f495394b80f108", db = "col", downto = "species")
#' ## ITIS tsn ids can be numeric or character
#' downstream("154395", db = "itis", downto = "species")
#' downstream(154395, db = "itis", downto = "species")
#'
#' # Plug in taxon names
#' downstream("Insecta", db = 'col', downto = 'order')
#' downstream("Apis", db = 'col', downto = 'species')
#' downstream("Apis", db = 'itis', downto = 'species')
#' downstream(c("Apis","Epeoloides"), db = 'itis', downto = 'species')
#' downstream(c("Apis","Epeoloides"), db = 'col', downto = 'species')
#' downstream("Ursus", db = 'gbif', downto = 'species')
#' downstream(get_gbifid("Ursus"), db = 'gbif', downto = 'species')
#'
#' # Plug in IDs
#' id <- get_colid("Apis")
#' downstream(id, downto = 'species')
#'
#' ## Equivalently, plug in the call to get the id via e.g., get_colid into downstream
#' identical(downstream(id, downto = 'species'),
#'          downstream(get_colid("Apis"), downto = 'species'))
#'
#' id <- get_colid("Apis")
#' downstream(id, downto = 'species')
#' downstream(get_colid("Apis"), downto = 'species')
#'
#' # Many taxa
#' sp <- names_list("genus", 3)
#' downstream(sp, db = 'col', downto = 'species')
#' downstream(sp, db = 'itis', downto = 'species')
#' downstream(sp, db = 'gbif', downto = 'species')
#'
#' # Both data sources
#' ids <- get_ids("Apis", db = c('col','itis'))
#' downstream(ids, downto = 'species')
#' ## same result
#' downstream(get_ids("Apis", db = c('col','itis')), downto = 'species')
#'
#' # Collect intermediate names
#' ## itis
#' downstream('Bangiophyceae', db="itis", downto="genus")
#' downstream('Bangiophyceae', db="itis", downto="genus", intermediate=TRUE)
#' downstream(get_tsn('Bangiophyceae'), downto="genus")
#' downstream(get_tsn('Bangiophyceae'), downto="genus", intermediate=TRUE)
#' ## col
#' downstream(get_colid("Animalia"), downto="class")
#' downstream(get_colid("Animalia"), downto="class", intermediate=TRUE)
#'
#' # Use the rows parameter
#' ## note how in the second function call you don't get the prompt
#' downstream("Poa", db = 'col', downto="species")
#' downstream("Poa", db = 'col', downto="species", rows=1)
#'
#' # use curl options
#' res <- downstream("Apis", db = 'col', downto = 'species', config=verbose())
#' res <- downstream("Apis", db = 'itis', downto = 'species', config=verbose())
#' res <- downstream("Ursus", db = 'gbif', downto = 'species', config=verbose())
#' }
downstream <- function(...){
  UseMethod("downstream")
}

#' @export
#' @rdname downstream
downstream.default <- function(x, db = NULL, downto = NULL, intermediate = FALSE, rows=NA, ...){
  nstop(downto, "downto")
  nstop(db)
  switch(db,
         itis = {
           id <- process_stream_ids(x, db, get_tsn, rows = rows, ...)
           setNames(downstream(id, downto = tolower(downto), intermediate = intermediate, ...), x)
         },
         col = {
           id <- process_stream_ids(x, db, get_colid, rows = rows, ...)
           setNames(downstream(id, downto = tolower(downto), intermediate = intermediate, ...), x)
         },
         gbif = {
           id <- process_stream_ids(x, db, get_gbifid, rows = rows, ...)
           setNames(downstream(id, downto = tolower(downto), intermediate = intermediate, ...), x)
         },
         stop("the provided db value was not recognised", call. = FALSE)
  )
}

process_stream_ids <- function(input, db, fxn, ...){
  g <- tryCatch(as.numeric(as.character(input)), warning = function(e) e)
  if (is(g, "numeric") || is.character(input) && grepl("[[:digit:]]", input)) {
    as_fxn <- switch(db, itis = as.tsn, col = as.colid, gbif = as.gbifid)
    as_fxn(input, check = FALSE)
  } else {
    eval(fxn)(input, ...)
  }
}

#' @export
#' @rdname downstream
downstream.tsn <- function(x, db = NULL, downto = NULL, intermediate = FALSE, ...) {
  fun <- function(y, downto, intermediate, ...){
    # return NA if NA is supplied
    if (is.na(y)) {
      NA
    } else {
		  itis_downstream(tsns = y, downto = downto, intermediate = intermediate, ...)
    }
  }
  out <- lapply(x, fun, downto = downto, intermediate = intermediate, ...)
  structure(out, class = 'downstream', db = 'itis', .Names = x)
}

#' @export
#' @rdname downstream
downstream.colid <- function(x, db = NULL, downto = NULL, intermediate = FALSE, ...) {
  fun <- function(y, downto, intermediate, ...){
    # return NA if NA is supplied
    if (is.na(y)) {
      NA
    } else {
      col_downstream(id = y, downto = downto, intermediate = intermediate, ...)
    }
  }
  out <- lapply(x, fun, downto = downto, intermediate = intermediate, ...)
  structure(simp(out), class = 'downstream', db = 'col')
}

#' @export
#' @rdname downstream
downstream.gbifid <- function(x, db = NULL, downto = NULL, intermediate = FALSE, ...) {
  fun <- function(y, downto, intermediate, ...){
    # return NA if NA is supplied
    if (is.na(y)) {
      NA
    } else {
      gbif_downstream(key = y, downto = downto, intermediate = intermediate, ...)
    }
  }
  out <- lapply(x, fun, downto = downto, intermediate = intermediate, ...)
  structure(out, class = 'downstream', db = 'gbif')
}

#' @export
#' @rdname downstream
downstream.ids <- function(x, db = NULL, downto = NULL, intermediate = FALSE, ...) {
  fun <- function(y, downto, intermediate, ...){
    # return NA if NA is supplied
    if (is.na(y)) {
      NA
    } else {
      downstream(y, downto = downto, intermediate = intermediate, ...)
    }
  }
  structure(lapply(x, fun, downto = downto, intermediate = intermediate, ...), class = 'downstream_ids')
}

simp <- function(x) if (length(x) == 1) x[[1]] else x
