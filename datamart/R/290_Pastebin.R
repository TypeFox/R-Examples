#' Pastebin
#' 
#' This class exposes partially the Web API to the pastebin service.
#' 
#' @examples
#' getSlots("Pastebin")
#'
#' @seealso \code{\link{pastebin}}
#'
#' @name Pastebin-class
#' @rdname Pastebin-class
#' @exportClass Pastebin
setClass(
    Class="Pastebin", 
    representation=representation(
      api_dev_key = "character",
      api_user_key = "character",
      api_user_name= "character",
      curl.handle="ANY" # RCurl:::CURLHandle is not exported
    ),
    contains=c("Location", "Xdata")
)

#' Constructor for Pastebin Location object
#'
#' see Pastebin class for more information
#'
#' @param api_dev_key    API Dev Key, default getOption("pastebin.api_dev_key")
#' @param api_user_name  API User Name, default getOption("pastebin.api_user_name")
#' @param api_user_password API User password, default getOption("pastebin.api_user_password")
#' @param clss           Class name to initiate, default "Pastebin"
#'
#' @rdname Pastebin-class
#' @export
pastebin <- function(
  api_dev_key = getOption("pastebin.api_dev_key"),
  api_user_name = getOption("pastebin.api_user_name"),
  api_user_password = getOption("pastebin.api_user_password"),
  clss="Pastebin"
) {
  curlHandle = RCurl::getCurlHandle()
  api_user_key = RCurl::postForm("http://pastebin.com/api/api_login.php",
    api_dev_key = api_dev_key,
    api_user_name = api_user_name,
    api_user_password = api_user_password,
    curl = curlHandle
  )
  new(clss, api_user_key=api_user_key, api_dev_key=api_dev_key, api_user_name=api_user_name, curl.handle=curlHandle)
}


#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,Pastebin-method
setMethod(
  f="show",
  signature="Pastebin",
  definition=function(object) cat(sprintf("<%s @ Pastebin.com>\n", object@api_user_name))
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,Pastebin-method
setMethod(
  f="meta",
  signature="Pastebin",
  definition=function(self) {
    ans <- RCurl::postForm("http://pastebin.com/api/api_post.php",
      api_option="list",
      api_dev_key = self@api_dev_key,
      api_user_key = self@api_user_key,
      api_results_limit = 999,
      curl = self@curl.handle
    )
    ans <- XML::htmlParse(ans)
    ans <- XML::getNodeSet(ans, "//paste")
    ans <- sapply(ans, function(n) c(
      XML::xmlValue(n[["paste_key"]]),
      XML::xmlValue(n[["paste_date"]]),
      XML::xmlValue(n[["paste_title"]]),
      XML::xmlValue(n[["paste_size"]]),
      XML::xmlValue(n[["paste_expire_date"]]),
      XML::xmlValue(n[["paste_private"]]),
      XML::xmlValue(n[["paste_format_long"]]),
      XML::xmlValue(n[["paste_format_short"]]),
      XML::xmlValue(n[["paste_hits"]])
    ))
    ans <- as.data.frame(t(ans), stringsAsFactors=FALSE)
    colnames(ans) <- c("key", "date", "title", "size", "expire_date", "private", "format_long", "format_short", "hits")
    rownames(ans) <- ans[,"key"]
    ans <- ans[,-1]
    ans <- transform(ans, date=as.POSIXct(as.numeric(date), origin = "1970-01-01", tz = "GMT")) #maybe off one hour
    return(ans)
  }
)


#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,Pastebin,character-method
setMethod(
  f="query",
  signature=c(self="Pastebin", resource="character"),
  definition=function(self, resource, verbose=getOption("verbose"), ...) {
      dl <- RCurl::getURL(paste("http://pastebin.com/raw.php?i=", resource, sep=""))
      return(dl)
      if(resource %in% ls(envir=self@data_env))
        self@data_env[[resource]]
      else
        callNextMethod()
  }
)

