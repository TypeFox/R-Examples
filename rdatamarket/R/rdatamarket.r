library(zoo)
library(RCurl)
library(RJSONIO)

api_base <- "http://datamarket.com"
path_info <- "/api/v1/info.json"
path_series <- "/api/v1/series.csv"
path_list <- "/api/v1/list.csv"
short_url_services <- c(
  "http://data.is",
  "http://bit.ly",
  "http://is.gd",
  "http://t.co",
  "http://url.is"
)

#' Initialize DataMarket client with an API key
#'
#' An API key is required to access locked and non-public data. Access
#' to public datasets does not require an API key, for now at least.
#'
#' @param api.key a DataMarket API key, found in
#'   \href{http://datamarket.com/accounts/profile/}{your account profile},
#'   or NULL to forget the previously entered API key.
#' @export
#' @examples
#' dminit("1234567890abcdef1234567890abcdef")
#' dminit(NULL)
dminit <- function(api.key) {
  if (!missing(api.key)) {
    .rdatamarketEnv$api.key <- api.key
  }
}

#' Fetch information about a DataMarket dataset.
#'
#' This function makes an \code{info} API call to fetch metadata (including
#' dimension structure) of a DataMarket dataset.
#'
#' @param ds a dataset ID, DS string, URL query-string, or whole URL. The DS
#'           string to send is extracted from the URL as needed, and short URLs
#'           at data.is, bit.ly, is.gd, t.co and url.is are expanded.
#' @param .params extra GET parameters to pass along in the API request.
#' @return a structure of named lists representing the dataset metadata.
#'         If the resolved DS string is a multiple-dataset DS string (that is,
#'         contains / forward-slash characters), then the returned value is a
#'         named list of such structures, with `strsplit(ds, '/')` as names.
#' @export
#' @examples
#' dminfo("17tm")
#' dminfo("17tm!kqc=a")
#' dminfo("ds=17tm")
#' dminfo("ds=17tm!kqc=a")
#' dminfo("foo=bar&ds=17tm&baz=xyzzy")
#' dminfo("http://datamarket.com/api/v1/series.json?foo=bar&ds=17tm&baz=xyzzy")
#' dminfo("http://datamarket.com/data/set/17tm/#ds=17tm")
#' dminfo("http://datamarket.com/data/set/12r9/male-population-thousands")
dminfo <- function(ds, .params=list()) {
  return(dodminfo(ds, .params))
}

dodminfo <- function(ds, .params=list(), .curl=dmCurlHandle()) {
  origds <- ifelse(!is.na(ds) && grep("^https?:", ds), ds, NA)
  ctx <- interpret_ds(ds, .curl=.curl)
  if ("infos" %in% names(ctx)) {
    return(ctx$infos)
  }
  curlopts = list()
  if (!is.na(origds)) {
    curlopts$Referer <- origds
  }
  infojson <- getForm(
    paste(ctx$base, path_info, sep=""),
    curl=.curl,
    .opts=curlopts,
    .params=c(ctx$qs, callback="", .params=.params)
    )
  if (is.raw(infojson)) {
    infojson <- rawToChar(infojson)
  }
  infolist <- fromJSON(infojson, simplify=FALSE)
  names(infolist) <- lapply(infolist, FUN=function(i) i$ds)
  for (name in names(infolist)) {
    if (infolist[[name]]$status != "success") {
      class(infolist[[name]]) <- c('list', 'dmerror')
      next
    }
    names(infolist[[name]]$dimensions) <- lapply(infolist[[name]]$dimensions,
      FUN=function(dim) dim$id
    )
    for (dimname in names(infolist[[name]]$dimensions)) {
      class(infolist[[name]]$dimensions[[dimname]]) <- c("list", "dmdimension");
    }
    for (dimid in names(infolist[[name]]$dimensions)) {
      names(infolist[[name]]$dimensions[[dimid]]$values) <- lapply(
        infolist[[name]]$dimensions[[dimid]]$values,
        FUN=function(dimvalue) dimvalue[['id']]
      )
      class(infolist[[name]]$dimensions[[dimid]]$values) <- c("list",
                                                              "dmdimvalues");
      hier <- infolist[[name]]$dimension[[dimid]]$type == "hierarchical"
      if (hier) {
        depth <- list()
        for (dimvalueid in names(infolist[[name]]$dimensions[[dimid]]$values)) {
          class(infolist[[name]]$dimensions[[dimid]]$values[[dimvalueid]]) <-
            c("list", "dmhierarchicaldimvalue", "dmdimvalue")
          pid <- infolist[[name]]$dimensions[[dimid]]$values[[
            dimvalueid]][["parent_id"]]
          thisdepth <- ifelse(is.null(pid), 0, depth[[pid]] + 1);
          depth[[dimvalueid]] <- thisdepth
          attr(infolist[[name]]$dimensions[[dimid]]$values[[dimvalueid]],
               'depth') <- thisdepth
        }
      } else {
        for (dimvalueid in names(infolist[[name]]$dimensions[[dimid]]$values)) {
          class(infolist[[name]]$dimensions[[dimid]]$values[[dimvalueid]]) <-
            c("list", "dmdimvalue")
        }
      }
    }
    class(infolist[[name]]) <- c('list', 'dmdataset')
  }
  if (length(infolist) == 1) {
    infolist <- infolist[[1]];
  }
  return(infolist)
}

#' Fetch dimensions of a DataMarket dataset.
#'
#' Get a list of dataset dimension objects for the given dataset.
#'
#' @param ds a dataset ID, DS string, URL query-string, or whole URL. The DS
#'           string to send is extracted from the URL as needed, and short URLs
#'           at data.is, bit.ly, is.gd, t.co and url.is are expanded.
#'           If the DS string contains dimension filter specifications (the
#'           stuff after the ! character, so it's not just a dataset ID), these
#'           are preserved in the request to the API, but for normal DataMarket
#'           datasets they do not affect the response.
#' @param .params extra GET parameters to pass along in the API request.
#' @return named list of dataset dimension information. Each name is a dataset
#'         ID and each element is a named list of dimensions of that dataset.
#'         Each dimension is named for its dimension ID in that list, and is
#'         itself a named list of the four properties \code{id, title, type,
#'         values}. The first three of these properties are character strings,
#'         while \code{values} is a named list of dimension values. Each of
#'         these is a list of two properties \code{id, title}, and the \code{id}
#'         is also the name of the dimension value
#' @export
#' @examples
#' dmdims("17tm")
#' dmdims("17tm!kqc=a")
#' dmdims("ds=17tm")
#' dmdims("ds=17tm!kqc=a")
#' dmdims("foo=bar&ds=17tm&baz=xyzzy")
#' dmdims("http://datamarket.com/api/v1/series.json?foo=bar&ds=17tm&baz=xyzzy")
#' dmdims("http://datamarket.com/data/set/17tm/#ds=17tm")
dmdims <- function(ds, .params=list()) {
  infolist <- dminfo(ds, .params=.params)
  if ("dmdataset" %in% class(infolist)) {
    return(infolist$dimensions)
  } else {
    return(lapply(infolist, FUN=function(info) info$dimensions))
  }
}

#' Fetch timeseries from a DataMarket dataset.
#'
#' This performs a \code{series} API request at DataMarket.com, fetching the
#' requested data and wrapping it in a \code{zoo} object.
#'
#' @param ds a dataset ID, DS string, URL query-string, or whole URL. The DS
#'           string to send is extracted from the URL as needed, and short URLs
#'           at data.is, bit.ly, is.gd, t.co and url.is are expanded.
#' @param .params extra GET parameters to pass along in the API request.
#' @param ... named parameters whose names are dimension titles or IDs, and
#'            whose values are titles or IDs of values of those dimensions.
#'            E.g. if dataset \code{17tm} has a dimension named Country, then
#'            \code{dmseries("17tm", Country="Algeria")} filters on that
#'            dimension. If the dimension name includes spaces, it needs to be
#'            quoted: \code{dmlist("12rb", "Country or Area"="Afghanistan")}
#' @return a zoo object representing the fetched timeseries.
#' @export
#' @examples
#' dmseries("17tm")
#' dmseries("17tm!kqc=a")
#' dmseries("ds=17tm")
#' dmseries("ds=17tm!kqc=a")
#' dmseries("foo=bar&ds=17tm&baz=xyzzy")
#' dmseries("http://datamarket.com/api/v1/series.json?foo=bar&ds=17tm&baz=xyzzy")
#' dmseries("http://datamarket.com/data/set/17tm/#ds=17tm")
#' dmseries("17tm", Country="Algeria")
#' dmseries("17tm", Country=c("Algeria", "Angola"))
dmseries <- function(ds, .params=list(), ...) {
  origds <- ifelse(!is.na(ds) && grep("^https?:", ds), ds, NA)
  curl <- dmCurlHandle()
  ctx <- interpret_ds(ds, .curl=curl)
  if (!(identical(c(...), c()))) {
    infos <- dodminfo(ds, .curl=curl)
    ctx$qs$ds <- dimfilter(ctx$qs$ds, infos, ...)
  }
  csv <- get.datamarket.csv(ctx, path_series, curl, .params, origds)
  for (name in names(csv)) {
    if (all(is.na(csv[[name]]))) {
      csv[[name]] <- NULL
    }
  }
  timecolname <- NULL
  for (name in c("Date", "Year.and.week", "Year.and.month", "Year.and.quarter", "Year")) {
    if (name %in% names(csv)) {
      timecolname <- name
      break
    }
  }
  if (is.null(timecolname)) {
    stop(paste("No time column found in timeseries response. Columns are:",
               names(csv)))
  }
  zindex <- switch(timecolname,
    Year=as.integer(csv[[timecolname]]),
    Year.and.quarter=as.yearqtr(csv[[timecolname]]),
    Year.and.month=as.yearmon(csv[[timecolname]]),
    Date=as.Date(csv[[timecolname]]),
    csv[[timecolname]]
    )

  csv[[timecolname]] <- NULL
  zoo(csv, order.by=zindex)
}

#' Fetch data in long form as a data.frame from a DataMarket dataset.
#'
#' This performs a \code{list} API request at DataMarket.com, fetching the
#' requested data and wrapping it in a \code{data.frame} object.
#'
#' @param ds a dataset ID, DS string, URL query-string, or whole URL. The DS
#'           string to send is extracted from the URL as needed, and short URLs
#'           at data.is, bit.ly, is.gd, t.co and url.is are expanded.
#' @param .params extra GET parameters to pass along in the API request.
#' @param ... named parameters whose names are dimension titles or IDs, and
#'            whose values are titles or IDs of values of those dimensions.
#'            E.g. if dataset \code{17tm} has a dimension named Country, then
#'            \code{dmlist("17tm", Country="Algeria")} filters on that
#'            dimension. If the dimension name includes spaces, it needs to be
#'            quoted: \code{dmlist("12rb", "Country or Area"="Afghanistan")}
#' @return a data frame representing the fetched data. The data frame has one
#'         column for each of the dataset's dimensions, containing a factor
#'         whose levels are values of that dimension, and a `Value` column
#'         holding a numerical value for each combination of dimension values.
#' @export
#' @examples
#' dmlist("17tm")
#' dmlist("17tm!kqc=a")
#' dmlist("ds=17tm")
#' dmlist("ds=17tm!kqc=a")
#' dmlist("foo=bar&ds=17tm&baz=xyzzy")
#' dmlist("http://datamarket.com/api/v1/series.json?foo=bar&ds=17tm&baz=xyzzy")
#' dmlist("http://datamarket.com/data/set/17tm/#ds=17tm")
#' dmlist("17tm", Country="Algeria")
#' dmlist("17tm", Country=c("Algeria", "Angola"))
#' dmlist("12rb", "Country or Area"="Afghanistan")
dmlist <- function(ds, .params=list(), ...) {
  origds <- ifelse(!is.na(ds) && grep("^https?:", ds), ds, NA)
  curl <- dmCurlHandle()
  ctx <- interpret_ds(ds, .curl=curl)
  if (!(identical(c(...), c()))) {
    infos <- dodminfo(ds, .curl=curl)
    ctx$qs$ds <- dimfilter(ctx$qs$ds, infos, ...)
  }
  get.datamarket.csv(ctx, path_list, curl, .params, origds)
}

#' @S3method [ dmdimvalues
`[.dmdimvalues` <- function (x, i)  {
  y <- unclass(x)[i]
  class(y) <- class(x)
  return (y)
}

#' @S3method format dmhierarchicaldimvalue
format.dmhierarchicaldimvalue <- function(x, ...) {
  sprintf('%s"%s"',
          paste(rep("-> ", attr(x, "depth")), collapse=""),
          x[["title"]]
  )
}

#' @S3method format dmdimvalue
format.dmdimvalue <- function(x, ...) {
  sprintf('"%s"', x[["title"]])
}

#' @S3method format dmdimvalues
format.dmdimvalues <- function(x, ...) {
    paste(lapply(as.list(x),
        FUN=function(x) sprintf('%3s  %s', x[["id"]], format(x))
      ),
      collapse="\n");
}

#' @S3method format dmdimension
format.dmdimension <- function(x, ...) {
  N <- length(x$values);
  sprintf('"%s" (%d values):\n    %s%s',
    x$title,
    N,
    paste(lapply(as.list(x$values)[1:min(5,N)], FUN=format), collapse="\n    "),
    ifelse(N > 5, "\n    [...]", "")
    )
}

#' @S3method format dmdataset
format.dmdataset <- function(x, ...) {
  sprintf("Title: \"%s\"\nProvider: \"%s\"%s\nDimensions:\n  %s",
    x$title,
    x$meta$provider_title,
    ifelse(!is.null(x$meta$source_source) && x$meta$source_source != "",
           sprintf(" (citing \"%s\")", x$meta$source_source),
           ""),
    paste(lapply(x$dimensions, FUN=format), collapse="\n  ")
    )
}

#' @S3method print dmhierarchicaldimvalue
print.dmhierarchicaldimvalue <- function(x, ...) {
  cat(format(x, ...), "\n");
  invisible();
}

#' @S3method print dmdimvalue
print.dmdimvalue <- function(x, ...) {
  cat(format(x, ...), "\n");
  invisible();
}

#' @S3method print dmdimvalues
print.dmdimvalues <- function(x, ...) {
  cat(format(x, ...), "\n");
  invisible();
}

#' @S3method print dmdimension
print.dmdimension <- function(x, ...) {
  cat(format(x), "\n");
  invisible();
}

#' @S3method print dmdataset
print.dmdataset <- function(x, quote=FALSE, ...) {
  cat(format(x), "\n");
  invisible();
}

#' @S3method print dmerror
print.dmerror <- function(x, ...) {
  cat(sprintf("%s: %s\n", x$status, x$message))
  invisible();
}

#' Set extra RCurl options for full HTTP control.
#'
#' Use this to control the options with which RCurl connections are created, see
#' \href{http://www.omegahat.org/RCurl/installed/RCurl/html/curlOptions.html}{RCurl::curlOptions}.
#' A common use is setting a proxy to work with your company's firewall setup,
#' see example below.
#'
#' @param ... name-value pairs specifying curl options. See full list of options
#'            with `names(getCurlOptionsConstants())`, their types with
#'            `getCurlOptionTypes()`, and full documentation in
#'            \href{http://curl.haxx.se/docs/manpage.html}{the curl manpage}.
#' @param .opts a named list of options, or a previously created `CURLOptions`
#'              object. These are merged with the options specified in `...`.
#' @export
#' @examples
#' \dontrun{dmCurlOptions(proxy='http://outproxy.mycompany.com')}
dmCurlOptions <- function(..., .opts=list()) {
  .rdatamarketEnv$curlopts <- curlOptions(.opts=.opts, ...)
}

getUserAgent <- function() {
  if (is.null(.rdatamarketEnv$user.agent)) {
    rdmv <- tryCatch(sessionInfo('rdatamarket')$otherPkgs$rdatamarket$Version,
      warning=function(cond) '(devtools)')
    si <- sessionInfo("RCurl")
    .rdatamarketEnv$user.agent <- paste(
      "rdatamarket ",
      rdmv,
      "; ",
      si$R.version$version.string,
      "; RCurl ",
      si$otherPkgs$RCurl$Version,
      "; platform ",
      si$R.version$platform,
      sep=""
      )
  }
  return(.rdatamarketEnv$user.agent)
}

dmCurlHandle <- function() {
  getCurlHandle(
    .opts=.rdatamarketEnv$curlopts,
    httpheader=c(
      .rdatamarketEnv$curlopts$httpheader,
      'User-Agent'=getUserAgent(),
      'X-DataMarket-API-Key'=.rdatamarketEnv$api.key
    )
  )
}
