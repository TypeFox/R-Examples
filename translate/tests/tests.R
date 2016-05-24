library(translate)

if (is.null(get.key())) {
  message("No API key is set; not running tests.")
} else {
  stopifnot(detect.source('hello') == 'en')
  stopifnot(length(languages(target='de')) == 53)
  stopifnot(translate('hello, world', 'en', 'de') == 'Hallo Welt!')
}
