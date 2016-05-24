as.WebResult <- function(x, ...) UseMethod("as.WebResult")

as.WebResult.WebResult <- function(x, ...) x

as.WebResult.default <- function(x, ...)
  WebResult("html", paste(as.character(x), collapse='\n'))

WebResult <- function(cmd="html", payload="", content.type="text/html; charset=utf-8", headers=character(0))
  structure(c(cmd, paste(as.character(payload), collapse="\n"), content.type, paste(c(as.character(.e$headers), as.character(headers)),collapse="\r\n")), class="WebResult")
