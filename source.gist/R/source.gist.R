library(RCurl)
library(rjson)

source.gist <- function(id) {
  # Accept gist URLs as well as IDs.
  if (substr(id, 0, 24) == 'https://gist.github.com/') {
    id <- substring(id, 25)
  }

  gist.url <- sprintf("https://api.github.com/gists/%s", id)
  metadata <- fromJSON(getURL(gist.url, followlocation = TRUE))
  raw.url <- metadata$files[[1]]$raw_url
  data <- getURL(raw.url, followlocation = TRUE)
  eval(parse(text=data), envir= .GlobalEnv)
}
