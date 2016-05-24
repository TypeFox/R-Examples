## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(comment = "")
library(curl)

## ------------------------------------------------------------------------
req <- curl_fetch_memory("https://httpbin.org/get")
str(req)
parse_headers(req$headers)
cat(rawToChar(req$content))

## ------------------------------------------------------------------------
tmp <- tempfile()
curl_download("https://httpbin.org/get", tmp)
cat(readLines(tmp), sep = "\n")

## ------------------------------------------------------------------------
con <- curl("https://httpbin.org/get")
open(con)

# Get 3 lines
out <- readLines(con, n = 3)
cat(out, sep = "\n")

# Get 3 more lines
out <- readLines(con, n = 3)
cat(out, sep = "\n")

# Get remaining lines
out <- readLines(con)
close(con)
cat(out, sep = "\n")

## ------------------------------------------------------------------------
req <- curl_fetch_memory("https://httpbin.org/status/418")
print(req$status_code)

## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
invisible(gc())

## ------------------------------------------------------------------------
h <- new_handle()
handle_setopt(h, copypostfields = "moo=moomooo");
handle_setheaders(h,
  "Content-Type" = "text/moo",
  "Cache-Control" = "no-cache",
  "User-Agent" = "A cow"
)

## ------------------------------------------------------------------------
req <- curl_fetch_memory("http://httpbin.org/post", handle = h)
cat(rawToChar(req$content))

## ------------------------------------------------------------------------
con <- curl("http://httpbin.org/post", handle = h)
cat(readLines(con), sep = "\n")

## ------------------------------------------------------------------------
tmp <- tempfile()
curl_download("http://httpbin.org/post", destfile = tmp, handle = h)
cat(readLines(tmp), sep = "\n")

## ---- echo = FALSE, message = FALSE, warning=FALSE-----------------------
invisible(gc())

## ------------------------------------------------------------------------
# Start with a fresh handle
h <- new_handle()

# Ask server to set some cookies
req <- curl_fetch_memory("http://httpbin.org/cookies/set?foo=123&bar=ftw", handle = h)
req <- curl_fetch_memory("http://httpbin.org/cookies/set?baz=moooo", handle = h)
handle_cookies(h)

# Unset a cookie
req <- curl_fetch_memory("http://httpbin.org/cookies/delete?foo", handle = h)
handle_cookies(h)

## ------------------------------------------------------------------------
req1 <- curl_fetch_memory("https://httpbin.org/get", handle = new_handle())
req2 <- curl_fetch_memory("http://www.r-project.org", handle = new_handle())

## ------------------------------------------------------------------------
h <- new_handle()
system.time(curl_fetch_memory("https://api.github.com/users/ropensci", handle = h))
system.time(curl_fetch_memory("https://api.github.com/users/rstudio", handle = h))

## ------------------------------------------------------------------------
handle_reset(h)

## ------------------------------------------------------------------------
# Posting multipart
h <- new_handle()
handle_setform(h,
  foo = "blabla",
  bar = charToRaw("boeboe"),
  description = form_file(system.file("DESCRIPTION")),
  logo = form_file(file.path(Sys.getenv("R_DOC_DIR"), "html/logo.jpg"), "image/jpeg")
)
req <- curl_fetch_memory("http://httpbin.org/post", handle = h)

## ------------------------------------------------------------------------
library(magrittr)

new_handle() %>% 
  handle_setopt(copypostfields = "moo=moomooo") %>% 
  handle_setheaders("Content-Type" = "text/moo", "Cache-Control" = "no-cache", "User-Agent" = "A cow") %>%
  curl_fetch_memory(url = "http://httpbin.org/post") %$% content %>% rawToChar %>% cat

