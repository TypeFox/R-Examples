#' @export
#' @rdname get_gho_codes
filter_attrs <- function (x, ...) {
  .dots <- lazyeval::lazy_dots(...)

  tab_attrs <- attr(x, "attrs") %>%
    dplyr::filter_(.dots = .dots)

  codes <- unique(tab_attrs$code)

  build_gho(
    x[x %in% codes],
    labels = attr(x, "labels")[x %in% codes],
    attrs = tab_attrs
  )
}
