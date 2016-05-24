#' Util

ESTAT_API_URL <- "http://api.e-stat.go.jp/"

#' @title e-Stat API
#'
#' @description Get Statistical Something From e-Stat API
#'
#' @param path API endpoint
#' @param appId application ID
#' @param ... other parameters
#'
#' @export
estat_api <- function(path, appId, ...)
{
  # convert params like cdCat01=c("001","002") to cdCat01="001,002"
  query <- flatten_query(list(appId = appId, ...))

  res <- httr::GET(
    ESTAT_API_URL,
    path = path,
    query = query
  )

  httr::warn_for_status(res)

  result_json <- httr::content(res)

  if(result_json[[1]]$RESULT$STATUS != 0) stop(result_json[[1]]$RESULT$ERROR_MSG)

  result_json
}

flatten_query <- function(x)
{
  purrr::map(x, ~ paste0(as.character(.), collapse = ","))
}

as_flattened_character <- function(x)
{
  purrr::map(purrr::flatten(x), as.character)
}

force_bind_rows <- function(x)
{
  if(is.list(x[[1]])) dplyr::bind_rows(x) else dplyr::as_data_frame(x)
}

get_class_info <- function(class_obj)
{
  class_info <- purrr::map(class_obj, ~ force_bind_rows(.$CLASS))
  names(class_info) <- purrr::map_chr(class_obj, ~ .$`@id`)
  class_info
}

merge_class_info <- function(value_df, class_info, name)
{
  info <- class_info[[name]] %>%
    dplyr::select_("`@code`", "`@name`")

  key <- sprintf("@%s", name)
  colnames(info) <- c(key, sprintf("%s_info", name))
  dplyr::left_join(value_df, info, by = key)
}
