"%&&%" <- function(a, b) {
  if (is.null(a)) a else force(b)
}

"%||%" <- function(a, b) {
  if (is.null(a)) force(b) else a
}

"%|na|%" <- function(a, b) {
  if (is.null(a) || all(is.na(a))) force(b) else a
}

"%|char|%" <- function(a, b) {
  if (all(!nzchar(a))) force(b) else a
}

"%|empty|%" <- function(a, b) {
  if (length(a) == 0L) force(b) else a
}

"%ni%" <- Negate("%in%")

.escape <- function(x) {
  x <- gsub("\\s+", " ", x)
  gsub(" (and) | (or) | (not) "," \\U\\1\\U\\2\\U\\3 ", x, perl = TRUE)
}

.collapse <- function(id) {
  if (is.null(id)) NULL else paste0(id, collapse = ",")
}

merge_list <- function(x, y) {
  if (length(x) == 0) return(y)
  if (length(y) == 0) return(x) 
  i <- is.na(match(names(y), names(x)))
  if (any(i)) {
    x[names(y)[which(i)]] <- y[which(i)]
  }
  x
}

merge_linkset <- function(x) {
  uid <- unique(unlist(x, use.names = FALSE))
  db <- unique(vapply(x, attr, "database", FUN.VALUE = "", USE.NAMES = FALSE))
  assertthat::assert_that(assertthat::is.scalar(db))
  attr(uid, "database") <- db
  uid
}

compact <- function(x) {
  x[!vapply(x, is.null, FALSE, USE.NAMES = FALSE)]
}

compactNA <- function(x) {
  x[!vapply(x, function(x) all(is.na(x)), FALSE, USE.NAMES = FALSE)]
}

count_char <- function(char, s) {
  s2 <- gsub(char, "", s)
  nchar(s) - nchar(s2)
}

#' Extract the content of XML leaf nodes
#' 
#' @param doc An object of class \code{XMLInternalDocument}.
#' @param path An XPath expression.
#' @param as Mode of return value (\code{character}, \code{integer}, \code{numeric}).
#' @param default Default return value.
#' @param ... Arguments passed to \code{\link[XML]{xpathApply}}.
#' @keywords internal
#' @export
xvalue <- function(doc, path, as = "character", default = NA_character_, ...) {
  AS <- match.fun(paste0('as.', as))
  res <- unlist(XML::xpathApply(doc, path, XML::xmlValue, ...)) %||% default
  res %&&% AS(res)
}

#' Extract the tag name of XML nodes
#' 
#' @inheritParams xvalue
#' @keywords internal
#' @export
xname <- function(doc, path, as = 'character', default = NA_character_, ...) {
  AS <- match.fun(paste0('as.', as))
  res <- unlist(XML::xpathApply(doc, path, XML::xmlName, ...)) %||% default
  res %&&% AS(res)
}

#' Extract the attributes of XML nodes
#' 
#' @inheritParams xvalue
#' @param name Name of the attribute to be extracted.
#' @keywords internal
#' @export
xattr <- function(doc, path, name, as = 'character', default = NA_character_, ...) {
  AS <- match.fun(paste0('as.', as))
  res <- unlist(XML::xpathApply(doc, path, XML::xmlGetAttr, name = name, ...)) %||% default
  res %&&% AS(res)
}

#' Extract a node set from an XML document
#' 
#' @param doc An object of class \code{XMLInternalDocument}.
#' @param path An XPath expression.
#' @param ... Arguments passed to \code{\link[XML]{xpathApply}}.
#' 
#' @keywords internal
#' @export
xset <- function(doc, path, ...) {
  XML::xpathApply(doc, path, fun = NULL, ...)
}

has_attr <- function(x, which) {
  a <- attr(x, which, exact = TRUE)
  !is.null(a) && !is.na(a) 
}

ellipsize <- function(obj, offset = 0, width = getOption("width"), ellipsis = "...") {
  str <- encodeString(obj)
  ifelse(nchar(str) > width - nchar(ellipsis) - offset,
         paste0(substring(str, 1, width - nchar(ellipsis) - offset), ellipsis),
         str)
}

#' Flatten nested lists
#'
#' Generate a function that accepts an arbitrarily deeply nested \code{list}.
#' Use \code{flatten.at} to set the level of nestedness at which the flatterner
#' will start to flatten.
#'
#' @param flatten.at An \code{integer} specifying the layer after which to 
#' start the flattening. \code{1} means to start at the very top.
#' @return A function that can be used to flatten a nested list
#' @keywords internal
#' @importFrom stats setNames
#' @export
make_flattener <- function(flatten.at = 1) {
  level <- 1
  .flatten <- function(x) {
    nm <- names(x)
    out <- list()
    for (i in seq_along(x)) {
      if (is.recursive(x[[i]])) {
        if (!is.null(nm[i])) {
          names(x[[i]]) <- paste0(nm[i], ".", names(x[[i]]))
        }
        level <<- level + 1
        out <- c(out, Recall(x[[i]]))
      } else {
        out <- c(out, x[i])
      }
    }
    level <<- level - 1
    out
  }
  
  function(x) {
    nm <- names(x)
    out <- list()
    for (i in seq_along(x)) {
      if (is.recursive(x[[i]]) && level >= flatten.at) {
        if (!is.null(nm[i])) {
          names(x[[i]]) <- paste0(nm[i], ".", names(x[[i]]))
        }
        level <<- level + 1
        out <- c(out, .flatten(x[[i]]))
      } else if (is.recursive(x[[i]]) && level < flatten.at) {
        level <<- level + 1
        out <- c(out, setNames(list(Recall(x[[i]])), nm[i]))
      } else {
        out <- c(out, x[i])
      }
    }
    if (level > 1) {
      level <<- level - 1
    }
    out
  }
}

flatten2 <- make_flattener(flatten.at = 2)

#' Set the NCBI rettype
#' 
#' @param db A valid NCBI database.
#' @param rettype Optional.
#' @param retmode Optional.
#' @keywords internal
#' @export
ncbi_retrieval_type <- function(db, rettype = NULL, retmode = NULL) {
  if (is.null(rettype) && !is.null(retmode)) {
    stop("No retrieval type specified", call. = FALSE)
  }
  rt <- set_rettype(db, rettype)
  rm <- set_retmode(db, rt, retmode)
  list(rettype = rt %||% "", retmode = rm)
}

set_rettype <- function(db, rt = NULL) {
  db <- switch(db, nucleotide = 'nuccore', db)
  rt  <- rt %|char|% NULL
  switch(db,
         bioproject = match.arg(rt, c("xml", "docsum", "uilist")),
         biosample  = match.arg(rt, c("full", "docsum", "uilist")),
         biosystems = match.arg(rt, c("xml", "docsum", "uilist")),
         gds        = match.arg(rt, c("summary", "docsum", "uilist")),
         gene       = rt %&&% match.arg(rt, c("gene_table", "docsum", "uilist")),
         homologene = rt %&&% match.arg(rt, c("alignmentscores", "fasta", "homologene", "docsum", "uilist")),
         mesh       = match.arg(rt, c("full", "docsum", "uilist")),
         nlmcatalog = rt %&&% match.arg(rt, c("docsum", "uilist")),
         nuccore    = rt %&&% match.arg(rt, c("fasta", "acc", "seqid", "native", "gb", "gbc",
                                              "gbwithparts", "fasta_cds_na", "fasta_cds_aa",
                                              "ft", "docsum", "uilist")),
         nucest     = rt %&&% match.arg(rt, c("fasta", "acc", "seqid", "native", "gb", "gbc",
                                              "est", "docsum", "uilist")),
         nucgss     = rt %&&% match.arg(rt, c("fasta", "acc", "seqid", "native", "gb", "gbc",
                                              "gss", "docsum", "uilist")),
         popset     = rt %&&% match.arg(rt, c("fasta", "acc", "seqid", "native", "gb", "gbc",
                                              "docsum", "uilist")),
         protein    = rt %&&% match.arg(rt, c("fasta", "acc", "seqid", "native", "gp", "gpc",
                                              "ipg", "ft", "docsum", "uilist")),
         pmc        = rt %&&% match.arg(rt, c("medline", "docsum", "uilist")),
         pubmed     = rt %&&% match.arg(rt, c("medline", "uilist", "abstract", "docsum")),
         sequences  = rt  %&&% match.arg(rt, c("acc", "fasta", "seqid", "uilist", "docsum")),
         snp        = rt %&&% match.arg(rt, c("flt", "fasta", "rsr", "ssexemplar", "chr",
                                              "docset", "uilist", "docsum")),
         sra        = match.arg(rt, c("full", "uilist", "docsum")),
         taxonomy   = rt %&&% match.arg(rt, c("uilist", "docsum")),
         stop('Database ', sQuote(db), ' not supported', call.=FALSE))
}

set_retmode <- function(db, rt, rm = NULL) {
  if (!is.null(rt) && rt == "docsum") {
    return("xml")
  }
  if (!is.null(rt) && rt == "uilist") {
    return(match.arg(rm, c("xml", "text")))
  }
  db <- switch(db, nucleotide = 'nuccore', db)
  switch(db,
         bioproject = switch(rt,
                             xml = match.arg(rm, c("xml"))
         ),
         biosample = switch(rt,
                            full = match.arg(rm, c("xml", "text"))
         ),
         biosystems = switch(rt,
                             xml = match.arg(rm, c("xml"))
         ),
         gds = switch(rt,
                      summary = match.arg(rm, c("text"))
         ),
         gene = switch(rt %||% 'null',
                       null = match.arg(rm, c("xml", "asn.1")),
                       gene_table = match.arg(rm, c("text"))
         ),
         homologene = switch(rt %||% 'null',
                             null = match.arg(rm, c("xml", "asn.1")),
                             alignmentscores = match.arg(rm, c("text")),
                             fasta = match.arg(rm, c("text")),
                             homologen = match.arg(rm, c("text"))
         ),
         mesh = switch(rt,
                       full = match.arg(rm, c("text"))
         ),
         nlmcatalog = switch(rt %||% 'null',
                             null = match.arg(rm, c("xml", "text"))
         ),
         nuccore = switch(rt %||% 'null',
                          null = match.arg(rm, c("text", "asn.1")),
                          native = match.arg(rm, c("xml")),
                          acc = match.arg(rm, c("text")),
                          fasta = match.arg(rm, c("xml", "text")),
                          seqid = match.arg(rm, c("text")),
                          gb = match.arg(rm, c("text", "xml")),
                          gbc = match.arg(rm, c("xml")),
                          ft = match.arg(rm, c("text")),
                          gbwithparts = match.arg(rm, c("text")),
                          fasta_cds_na = match.arg(rm, c("text")),
                          fasta_cds_aa = match.arg(rm, c("text"))
                          
         ),
         nucest = switch(rt %||% 'null',
                         null = match.arg(rm, c("text", "asn.1")),
                         native = match.arg(rm, c("xml")),
                         acc = match.arg(rm, c("text")),
                         fasta = match.arg(rm, c("xml", "text")),
                         seqid = match.arg(rm, c("text")),
                         gb = match.arg(rm, c("text", "xml")),
                         gbc = match.arg(rm, c("xml")),
                         est = match.arg(rm, c("xml"))
                         
         ),
         nucgss = switch(rt %||% 'null',
                         null = match.arg(rm, c("text", "asn.1")),
                         native = match.arg(rm, c("xml")),
                         acc = match.arg(rm, c("text")),
                         fasta = match.arg(rm, c("xml", "text")),
                         seqid = match.arg(rm, c("text")),
                         gb = match.arg(rm, c("text", "xml")),
                         gbc = match.arg(rm, c("xml")),
                         gss = match.arg(rm, c("xml"))
                         
         ),
         protein = switch(rt %||% 'null',
                          null = match.arg(rm, c("text", "asn.1")),
                          native = match.arg(rm, c("xml")),
                          acc = match.arg(rm, c("text")),
                          fasta = match.arg(rm, c("xml", "text")),
                          seqid = match.arg(rm, c("text")),
                          gp = match.arg(rm, c("text", "xml")),
                          gpc = match.arg(rm, c("xml")),
                          ipg = match.arg(rm, c("xml")),
                          ft = match.arg(rm, c("text"))
                          
         ),
         popset = switch(rt %||% 'null',
                         null = match.arg(rm, c("text", "asn.1")),
                         native = match.arg(rm, c("xml")),
                         acc = match.arg(rm, c("text")),
                         fasta = match.arg(rm, c("xml", "text")),
                         seqid = match.arg(rm, c("text")),
                         gb = match.arg(rm, c("text", "xml")),
                         gbc = match.arg(rm, c("xml"))
                         
         ),
         pmc = switch(rt %||% 'null',
                      null = match.arg(rm, c("xml")),
                      medline = match.arg(rm, c("text"))
         ),
         pubmed = switch(rt %||% 'null',
                         null = match.arg(rm, c("xml", "asn.1")),
                         medline = match.arg(rm, c("text")),
                         uilist = match.arg(rm, c("text")),
                         abstract = match.arg(rm, c("text"))               
         ),
         sequences = switch(rt %||% 'null',
                            null = match.arg(rm, c("text")),
                            acc = match.arg(rm, c("text")),
                            fasta = match.arg(rm, c("text")),
                            seqid = match.arg(rm, c("text"))
         ),
         snp = switch(rt %||% 'null',
                      null = match.arg(rm, c("xml", "asn.1")),
                      flt = match.arg(rm, c("text")),
                      fasta = match.arg(rm, c("text")),
                      rsr = match.arg(rm, c("text")),
                      ssexemplar = match.arg(rm, c("text")),
                      chr = match.arg(rm, c("text")),
                      docset = match.arg(rm, c("text")),
                      uilist = match.arg(rm, c("xml", "text"))  
         ),
         sra = switch(rt,
                      full = match.arg(rm, c("xml"))
         ),
         taxonomy = switch(rt %||% 'null',
                           null = match.arg(rm, c("xml")),
                           uilist = match.arg(rm, c("xml", "text"))
         )
  )
}
