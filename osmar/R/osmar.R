#' @include osm-descriptors.R
{}



print_header <- function(what, x) {
  cat(sprintf("%s object\n", what))
  cat(paste(x, " ", names(x), sep = "", collapse = ", "))
}



print_content <- function(what, x) {
  for ( n in names(x) ) {
    cat(sprintf("..$%s data.frame:", n), "\n")
    cat(paste(strwrap(paste(x[[n]], sep = "", collapse = ", "),
                      indent = 4, exdent = 4), collapse = "\n"), "\n")
  }
}



abbrev <- function(x, nchar) {
  if ( is.null(nchar) )
    return(x)

  ifelse(nchar(x) <= nchar, x, sprintf("%s...", substr(x, 1, nchar)))
}



#' Summarize osmar objects
#'
#' Summaries of \code{osmar}, \code{nodes}, \code{ways}, and
#' \code{relations} objects. Use these methods to get an overview of
#' the content.
#'
#' @param object An object (\code{osmar}, \code{nodes}, \code{ways},
#'   or \code{relations} for which a summary is desired
#' @param ... Ignored
#'
#' @return
#'   \code{summary.osmar} returns a list with the summaries for nodes,
#'   ways, and relations.
#'
#'   \code{summary.nodes}, \code{summary.ways},
#'   \code{summary.relations} all return a list with
#'
#'   \describe{
#'
#'     \item{\code{key}}{A contingency table of the counts of each key
#'       label; in descending order}
#'
#'     \item{\code{val}}{A contingency table of the counts of each
#'       value label; in descending order}
#'
#'     \item{\code{keyval}}{A contingency table of the counts greater
#'       zero of each combination of key and value labels; in
#'       descending order}
#'
#'   }
#'
#' @seealso \code{\link{osmar}}
#'
#' @method summary osmar
#'
#' @S3method summary osmar
summary.osmar <- function(object, ...) {
  ret <- list(nodes = summary(object$nodes),
              ways = summary(object$ways),
              relations = summary(object$relations),
              n = sapply(object, function(y) nrow(y[[1]])))

  structure(ret, class = c("summary.osmar", class(ret)))
}



#' @param x The computed summary object to print
#' @param max.print Maximum number of shown tags
#' @param nchar.value Number of shown characters of the value column
#' @method print summary.osmar
#' @rdname summary.osmar
#' @S3method print summary.osmar
print.summary.osmar <- function(x, max.print = 3, nchar.value = 20, ...) {
  cat(print_header("osmar", x$n), "\n\n")
  print(x$nodes, max.print, nchar.value)
  cat("\n\n")
  print(x$ways, max.print, nchar.value)
  cat("\n\n")
  print(x$relations, max.print, nchar.value)
}



#' @S3method print osmar
print.osmar <- function(x, ...) {
  cat(print_header("osmar", dim(x)), "\n")
  invisible(x)
}




### Summarize nodes: #################################################


#' @method summary nodes
#' @rdname summary.osmar
#' @S3method summary nodes
summary.nodes <- function(object, ...) {
  ret <- list(n = NULL, bbox = NULL, content = NULL,
              key = NULL, val = NULL, keyval = NULL)

  ret$n <- c(nodes = nrow(object$attrs),
             tags = nrow(object$tags))

  if ( ret$n["nodes"] > 0 ) {
    ret$bbox <- cbind(lat = range(object$attrs$lat),
                      lon = range(object$attrs$lon))
    rownames(ret$bbox) <- c("min", "max")

    ret$content <- sapply(object, names)

    ret$key <- sort(table(object$tags$k), decreasing = TRUE)
    ret$key <- data.frame(Key = names(ret$key),
                          Freq = unname(ret$key),
                          stringsAsFactors = FALSE)
    rownames(ret$key) <- NULL

    ret$val <- sort(table(object$tags$v), decreasing = TRUE)
    ret$val <- data.frame(Value = names(ret$val),
                          Freq = unname(ret$val),
                          stringsAsFactors = FALSE)
    rownames(ret$val) <- NULL

    ret$keyval <- as.data.frame(table(Key = object$tags$k,
                                      Value = object$tags$v))
    if ( length(ret$keyval$Freq) > 0 ){
      ret$keyval <- ret$keyval[ret$keyval$Freq > 0, ]
      ret$keyval <- ret$keyval[order(-ret$keyval$Freq), ]
    }
    rownames(ret$keyval) <- NULL
  }

  structure(ret, class = c("summary.nodes", class(ret)))
}



#' @method print summary.nodes
#' @rdname summary.osmar
#' @S3method print summary.nodes
print.summary.nodes <- function(x, max.print = 10, nchar.value = 20, ...) {
  cat(print_header("osmar$nodes", x$n), "\n")
  if ( x$n["nodes"] > 0 ) {
    cat("\n")
    cat(print_content("osmar$nodes", x$content), "\n")
    cat("Bounding box:\n")
    print(x$bbox)
    cat("\nKey-Value contingency table:\n")
    keyval <- x$keyval[seq(min(max.print, nrow(x$keyval))), ]
    levels(keyval$Value) <- abbrev(levels(keyval$Value), nchar.value)
    print(keyval)
  }

  invisible(x)
}



#' @S3method print nodes
print.nodes <- function(x, ...) {
  n <- c(nodes = nrow(x$attrs),
         tags = nrow(x$tags))

  cat(print_header("osmar$nodes", n), "\n")
  if ( n["nodes"] > 0 ) {
    b <- cbind(lat = range(x$attrs$lat),
               lon = range(x$attrs$lon))
    rownames(b) <- c("min", "max")
    cat("\n")
    print(b)
  }

  invisible(x)
}



### Summarizing ways: ##################################################



#' @method summary ways
#' @rdname summary.osmar
#' @S3method summary ways
summary.ways <- function(object, ...) {
  ret <- list(n = NULL, key = NULL, content = NULL,
              val = NULL, keyval = NULL)

  ret$n <- c(ways = nrow(object$attrs),
             tags = nrow(object$tags),
             refs = nrow(object$refs))

  if ( ret$n["ways"] > 0 ) {
    ret$key <- sort(table(object$tags$k), decreasing = TRUE)
    ret$key <- data.frame(Key = names(ret$key),
                          Freq = unname(ret$key),
                          stringsAsFactors = FALSE)
    rownames(ret$key) <- NULL

    ret$content <- sapply(object, names)

    ret$val <- sort(table(object$tags$v), decreasing = TRUE)
    ret$val <- data.frame(Value = names(ret$val),
                          Freq = unname(ret$val),
                          stringsAsFactors = FALSE)
    rownames(ret$val) <- NULL

    ret$keyval <- as.data.frame(table(Key = object$tags$k,
                                      Value = object$tags$v))
    ret$keyval <- ret$keyval[ret$keyval$Freq > 0, ]
    ret$keyval <- ret$keyval[order(-ret$keyval$Freq), ]
    rownames(ret$keyval) <- NULL
  }

  structure(ret, class = c("summary.ways", class(ret)))
}



#' @method print summary.ways
#' @rdname summary.osmar
#' @S3method print summary.ways
print.summary.ways <- function(x, max.print = 10, nchar.value = 20, ...) {
  cat(print_header("osmar$ways", x$n), "\n")
  if ( x$n["ways"] > 0 ) {
    cat("\n")
    cat(print_content("osmar$ways", x$content), "\n")
    cat("Key-Value contingency table:\n")
    keyval <- x$keyval[seq(min(max.print, nrow(x$keyval))), ]
    levels(keyval$Value) <- abbrev(levels(keyval$Value), nchar.value)
    print(keyval)
  }
  invisible(x)
}



#' @S3method print ways
print.ways <- function(x, ...) {
  n <- c(ways = nrow(x$attrs),
         tags = nrow(x$tags),
         refs = nrow(x$refs))

  cat(print_header("osmar$ways", n), "\n")

  invisible(x)
}




### Summarizing relations: ###########################################


#' @method summary relations
#' @rdname summary.osmar
#' @S3method summary relations
summary.relations <- function(object, ...) {
  ret <- list(n = NULL, key = NULL, content = NULL,
              val = NULL, keyval = NULL)

  ret$n <- c(relations = nrow(object$attrs),
             tags = nrow(object$tags),
             node_refs = sum(object$refs$type == "node"),
             way_refs = sum(object$refs$type == "way"))

  if ( ret$n["relations"] > 0 ) {
    ret$key <- sort(table(object$tags$k), decreasing = TRUE)
    ret$key <- data.frame(Key = names(ret$key),
                          Freq = unname(ret$key),
                          stringsAsFactors = FALSE)
    rownames(ret$key) <- NULL

    ret$content <- sapply(object, names)

    ret$val <- sort(table(object$tags$v), decreasing = TRUE)
    ret$val <- data.frame(Value = names(ret$val),
                          Freq = unname(ret$val),
                          stringsAsFactors = FALSE)
    rownames(ret$val) <- NULL

    ret$keyval <- as.data.frame(table(Key = object$tags$k,
                                      Value = object$tags$v))
    ret$keyval <- ret$keyval[ret$keyval$Freq > 0, ]
    ret$keyval <- ret$keyval[order(-ret$keyval$Freq), ]
    rownames(ret$keyval) <- NULL
  }

  structure(ret, class = c("summary.relations", class(ret)))
}



#' @method print summary.relations
#' @rdname summary.osmar
#' @S3method print summary.relations
print.summary.relations <- function(x, max.print = 10, nchar.value = 20, ...) {
  cat(print_header("osmar$relations", x$n), "\n")
  if ( x$n["relations"] > 0 ) {
    cat("\n")
    cat(print_content("osmar$relations", x$content), "\n")
    cat("Key-Value contingency table:\n")
    keyval <- x$keyval[seq(min(max.print, nrow(x$keyval))), ]
    levels(keyval$Value) <- abbrev(levels(keyval$Value), nchar.value)
    print(keyval)
  }

  invisible(x)
}



#' @S3method print relations
print.relations <- function(x, ...) {
  n <- c(relations = nrow(x$attrs),
         tags = nrow(x$tags),
         node_refs = sum(x$refs$type == "node"),
         way_refs = sum(x$refs$type == "way"))

  cat(print_header("osmar$relations", n), "\n")
  invisible(x)
}



### Combining osmar objects: #########################################

#' Combine osmar objects
#'
#' Combine two or more \code{\link{osmar}} objects.
#'
#' @param ... \code{\link{osmar}} objects to be concatenated
#'
#' @return An \code{\link{osmar}} object based on the provided objects
#'
#' @examples
#'   \dontrun{
#'     muc <- get_osm(center_bbox(11.575278, 48.137222, 200, 200))
#'     o1 <- subset(muc, node_ids = find(muc, node(tags(v == "Marienplatz"))))
#'     o2 <- subset(muc, ids = find_down(muc, way(c(96619179, 105071000))))
#'
#'     o1
#'     o2
#'     c(o1, o2)
#'   }
#'
#' @method c osmar
#'
#' @S3method c osmar
c.osmar <- function(...) {
  ## TODO: object[[1]] attributes?
  objects <- list(...)

  stopifnot(are_osmar(objects))

  c_parts <- function(w1, w2) {
    do.call(rbind, lapply(objects, "[[", c(w1, w2)))
  }

  objects[[1]]$nodes$attrs <- unique(c_parts("nodes", "attrs"))
  objects[[1]]$nodes$tags <- unique(c_parts("nodes", "tags"))

  objects[[1]]$ways$attrs <- unique(c_parts("ways", "attrs"))
  objects[[1]]$ways$tags <- unique(c_parts("ways", "tags"))
  #objects[[1]]$ways$refs <- unique(c_parts("ways", "refs"))
  #  unique does make trouble with as_sp_polygons
  objects[[1]]$ways$refs <- c_parts("ways", "refs")

  objects[[1]]$relations$attrs <- unique(c_parts("relations", "attrs"))
  objects[[1]]$relations$tags <- unique(c_parts("relations", "tags"))
  objects[[1]]$relations$refs <- unique(c_parts("relations", "refs"))

  objects[[1]]
}



### Dimensions of osmar objects: ######################################

#' Dimension of osmar objects
#'
#' @param x An \code{\link{osmar}} object
#'
#' @return
#'   A named vector with the number of nodes, ways and relations.
#'
#' @examples
#'   \dontrun{
#'     muc <- get_osm(center_bbox(11.575278, 48.137222, 200, 200))
#'     dim(muc)
#'   }
#'
#' @method dim osmar
#'
#' @S3method dim osmar
dim.osmar <- function(x) {
  sapply(x, function(y) nrow(y[[1]]))
}

