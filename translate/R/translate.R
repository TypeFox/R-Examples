##' Bindings for the Google Translate API v2
##' @name translate-package
##' @docType package
##' @title translate
##' @examples \dontrun{
##'   set.key('YOUR-API-KEY')
##'   translate('Hello, world!', 'en', 'de')
##' }
NULL

##' @import RCurl RJSONIO lisp functional
NULL

google.base <- 'https://www.googleapis.com/language/translate/v2%s?key=%s%s%s%s'

##' Pull the API key from \code{getOption('google.key')} or the
##' \code{GOOGLE_KEY} environment variable.
##' @return The API key
##' @export
get.key <- function() {
  option.key <- getOption('google.key')
  env.key <- Sys.getenv('GOOGLE_KEY')
  if (!is.null(option.key))
    option.key
  else if (!env.key == '')
    env.key
  else
    NULL
}

##' Set the API key in \code{options}.
##' @param key The API key
##' @export
set.key <- function(key) options(google.key=key)

##' Generate an API URL.
##' @param method One of \code{"translate"} (default),
##' \code{"detect"}, \code{"languages"}
##' @param query The text to invoke against
##' @param source The source language, e.g. \code{"en"}
##' @param target The target language, e.g. \code{"de"}
##' @param key The API key
##' @return An API URL
google.url <- function(method=NULL,
                       query=NULL,
                       source=NULL,
                       target=NULL,
                       key=get.key())
  sprintf(google.base,
          ifelse(is.null(method), '', sprintf('/%s', method)),
          key,
          ifelse(is.null(query), '', sprintf('&q=%s', curlEscape(query))),
          ifelse(is.null(source), '', sprintf('&source=%s', source)),
          ifelse(is.null(target), '', sprintf('&target=%s', target)))

##' Both submit the get request and parse the JSON.
##' @noRd
getJSON <- Compose(getURL, fromJSON)

##' URL for source detection
##' @noRd
detect.source.url <- Curry(google.url, method='detect')

##' Detect the source of a text.
##' @inheritParams google.url
##' @return A list of potential source languages
##' @export
detect.source <- function(query, key=get.key())
  Map(function(detection) detection$language,
      car(getJSON(detect.source.url(query=query,
                                    key=key))$data$detections))

##' URL for language mapping
##' @noRd
languages.url <- Curry(google.url, method='languages')

##' List the valid language mappings; optionally, from a given source
##' or to a given target.
##' @inheritParams google.url
##' @param key The API key
##' @export
languages <- function(source=NULL, target=NULL, key=get.key())
  unlist(Map(as.vector,
             getJSON(languages.url(source=NULL,
                                   target=NULL,
                                   key=key))$data$languages))

##' URL for translatian
##' @noRd
translate.url <- Curry(google.url, method=NULL)

##' Translate the given text from a source to a target language.
##' @inheritParams google.url
##' @export
translate <- function(query, source, target, key=get.key())
  Map(as.vector,
      getJSON(translate.url(query, source, target, key))$data$translations)
