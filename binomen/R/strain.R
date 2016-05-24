#' Parse taxon or taxondf objects by a range of names
#'
#' @export
#' @param .data Input, object of class taxon
#' @param ... Logical predicates. Multiple conditions are combined with &.
#' See Details.
#' @return A single or list of \code{taxon} class objects
#' @details Example predicates:
#' \itemize{
#'  \item . > family = Get all taxa greater than family
#'  \item . < family = Get all taxa less than family
#'  \item . == family = Get all taxa equal to family
#'  \item . != family = Get all taxa not equal to family
#'  \item genus < order = Get all taxa between genus and order
#'  \item genus .. family = Get all taxa between genus and order
#' }
#' @examples
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#' out %>% strain(. < family)
#' out %>% strain(. < genus)
#' out %>% strain(. > family)
#' out %>% strain(. < family)
strain <- function(.data, ...) {
  UseMethod("strain")
}

#' @export
strain.taxon <- function(.data, ...) {
  tmp <- va_rs(.data, ...)
  strain_parse(w = .data, vars = tmp)
}

va_rs <- function(.data, ...) {
  va_rs_(.data, .dots = lazyeval::lazy_dots(...))
}

va_rs_ <- function(.data, ..., .dots) {
  lazyeval::all_dots(.dots, ...)
}

# strain.taxa <- function(.data, ...) {
#   lapply(.data, span, ...)
# }
#
# strain.taxondf <- function(.data, ...) {
#   var <- vars(...)
#   check_vars(var, names(.data))
#   matches <- sapply(var, grep, x = names(.data))
#   .data[fill_nums(matches)]
# }

strain_parse <- function(w, vars){
  grps <- w$grouping
  vars2 <- make_vars(vars)
  id1 <- rankid_get(vars2[[1]]$lower)
  id2 <- rankid_get(vars2[[1]]$upper)
  ids <- unlist(tc(list(id1, id2)))
  switch(vars2[[1]]$operator,
    `<` = {
      start <- which(rank_table$rankid == ids) + 1
      rks <- splitem(rank_table[seq(start, NROW(rank_table)), "ranks"])
      matches <- Filter(function(x) length(x) > 0, sapply(rks, grep, x = names(grps)))
      w$grouping <- do.call("grouping", grps[names(matches)])
      return(w)
    },
    `>` = {
      start <- which(rank_table$rankid == ids)
      rks <- splitem(rank_table[seq(1, start), "ranks"])
      matches <- Filter(function(x) length(x) > 0, sapply(rks, grep, x = names(grps)))
      w$grouping <- do.call("grouping", grps[names(matches)])
      return(w)
    }
  )
}

make_vars <- function(x) {
  lapply(x, function(z) {
    as.list(setNames(tolower(strsplit(deparse(z$expr), "\\s+")[[1]]),
             c('lower', 'operator', 'upper')))
  })
}

rankid_get <- function(x) {
  if (x == ".") {
    NULL
  } else {
    rank_table[rank_table$ranks %in% x, "rankid"]
  }
}

ranks_poss <- function() {
  sort(
    unique(
      do.call(c, sapply(rank_table$ranks, strsplit, split = ",", USE.NAMES = FALSE))
      )
  )
}

splitem <- function(x) {
  unlist(sapply(x, function(z) strsplit(z, ",")[[1]], USE.NAMES = FALSE))
}
