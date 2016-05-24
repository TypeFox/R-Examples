#' Set or query otions related to the rtematresdata R package.
#'
#' This function is used to query and set the options used by the rtematres
#' package. For example you can set the URLs to your tematres server and
#' the API.
#'
#' @param \dots similar to \code{\link{options}}. see examples below.
#' @examples
#' #Tematres URLs
#' rtematres.options('tematres_url')
#' rtematres.options(tematres_url="http://www.example.com")
#' @export rtematres.options

rtematres.options = function(...) {
  lst = list(...)
  .rtematres.opts = .rtematres.env$.rtematres.opts
  if (length(lst)) {
    if (is.null(names(lst)) && !is.list(lst[[1]])) {
      lst = unlist(lst)
      if (length(lst) == 1) .rtematres.opts[[lst]] else .rtematres.opts[lst]
    }
    else {
      omf = .rtematres.opts
      if (is.list(lst[[1]]))
        lst = lst[[1]]
      if (length(lst) > 0) {
        .rtematres.opts[names(lst)] = lapply(lst, gsub, pattern = "\\s", replacement="")
        if (!is.null(lst$url)) {
          .rtematres.opts["url"] = sub(.rtematres.opts["url"], pattern = "(/)?$", replacement = "")
          .rtematres.opts["url"] = sub(.rtematres.opts["url"], pattern = "^(http://)?", replacement = "http://")
        }
        if (!is.null(lst$tematres_url)) {
          .rtematres.opts["tematres_service_url"] = sub(.rtematres.opts["tematres_url"], pattern = "index.php/?$", replacement = "services.php")
        }
        if (!is.null(lst$tematres_service_url)) {
          .rtematres.opts["tematres_url"] = sub(.rtematres.opts["tematres_service_url"], pattern = "services.php/?$", replacement = "index.php")
        }
        .rtematres.env$.rtematres.opts = .rtematres.opts
      }
      invisible(omf)
    }
  }
  else {
    .rtematres.opts
  }
}

.rtematres.env = new.env()
