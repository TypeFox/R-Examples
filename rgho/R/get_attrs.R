#' Get Dimension Attributes
#'
#' @param xml_dim A parsed XML list of \code{Code} elements.
#'
#' @return A \code{data_frame} of attributes.
#'
get_attrs_ <- function(xml_dim) {

  codes <- xml_dim %>%
    xml2::xml_attr("Label")

  n_code <- unlist(lapply((xml_dim), xml2::xml_length)) - 1

  res <- dplyr::data_frame(
    code = rep(codes,
               n_code),
    key = xml_dim %>%
      xml2::xml_find_all("./Attr") %>%
      xml2::xml_attr("Category"),
    value = xml_dim %>%
      xml2::xml_find_all("./Attr/Value/Display") %>%
      xml2::xml_text()
  ) %>%
    tidyr::spread_(
      key_col = "key",
      value_col = "value"
    ) %>%
    dplyr::bind_rows(
      dplyr::data_frame(
        code = codes[n_code == 0]
      )
    )

  if (ncol(res) > 1) res else NULL
}

#' @rdname get_attrs_
get_attrs <- memoise::memoise(get_attrs_)
