#' Set or query otions related to the BEFdata R package.
#'
#' This function is used to query and set R BEFdata package specific options like
#' the URL to the BEFdata portal and the URL to the tematres thesaurus.
#'
#' @param \dots similar to \code{\link{options}}. see examples below.
#' @examples
#' # BEFdata URL
#' bef.options('url')
#' bef.options(url='http://www.example.com')
#' #Tematres URL
#' bef.options('tematres_url')
#' bef.options(tematres_url="http://www.example.com")
#'
#' @export

bef.options = function(...) {
  lst = list(...)
  .bef.opts = .bef.env$.bef.opts
  if (length(lst)) {
    if (is.null(names(lst)) && !is.list(lst[[1]])) {
      lst = unlist(lst)
      if (length(lst) == 1) .bef.opts[[lst]] else .bef.opts[lst]
    }
    else {
      omf = .bef.opts
      if (is.list(lst[[1]]))
        lst = lst[[1]]
      if (length(lst) > 0) {
        .bef.opts[names(lst)] = lapply(lst, gsub, pattern = "\\s", replacement="")
        if (!is.null(lst$url)) {
          .bef.opts["url"] = sub(.bef.opts["url"], pattern = "(/)?$", replacement = "")
          .bef.opts["url"] = sub(.bef.opts["url"], pattern = "^(http://)?", replacement = "http://")
        }
        if (!is.null(lst$tematres_url)) {
          .bef.opts["tematres_service_url"] = sub(.bef.opts["tematres_url"], pattern = "index.php/?$", replacement = "services.php")
        }
        if (!is.null(lst$tematres_service_url)) {
          .bef.opts["tematres_url"] = sub(.bef.opts["tematres_service_url"], pattern = "services.php/?$", replacement = "index.php")
        }
        .bef.env$.bef.opts = .bef.opts
      }
      invisible(omf)
    }
  }
  else {
    .bef.opts
  }
}

.bef.env = new.env()
