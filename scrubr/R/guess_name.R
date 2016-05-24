guess_name <- function(x, name = NULL) {
  nms <- names(x)

  if (is.null(name)) {
    zz <- nms[grep(sprintf("^(%s)$", paste0(name_options, collapse = "|")), nms, ignore.case = TRUE)]

    if (length(zz) == 1) {
      message("Assuming '", zz, "' is the taxonomic name field")
      names(x)[names(x) %in% zz] <- name_var <-  "name"
    } else {
      stop("Couldn't infer taxonomic name column, please specify with the 'name' parameter",
           call. = FALSE)
    }
  } else {
    if (!any(names(x) %in% name)) stop("'", name, "' not found in your data", call. = FALSE)
    names(x)[names(x) %in% name] <- name_var <- "name"
  }
  structure(x, name_var = name_var)
}

name_options <- c("name", "species", "taxonomic_name", "taxon")
