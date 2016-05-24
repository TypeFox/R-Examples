.run <- function(request, root, path) {
  as.character(try({
    .GlobalEnv$webapi <- 1.1
    .e$.out <- ''
    cmd <- 'html'
    ct <- 'text/html'
    setwd(file.path(root, "tmp"))

    # deprecated compatibility settings - should go away...
    requestURI <- request$uri
    remote.addr <- request$client.ip
    raw.cookies <- request$raw.cookies
    qs <- request$query.string

    # parse parameters in the order of precedence: query string, multipart, urlencoded
    pars <- list()
    if (request$c.type == 'application/x-www-form-urlencoded' && is.raw(request$body)) {
      ue <- rawToChar(request$body)
      for (x in strsplit(strsplit(ue,'&')[[1]], '=')) pars[[URLdecode(x[1])]] <- URLdecode(x[2])
    }
    if (grepl("^multipart", request$c.type)) pars <- parse.multipart()

    # add qs
    if (is.null(request$query.vector)) {
      for (x in strsplit(strsplit(qs,'&')[[1]], '=')) if (length(x) > 1L) pars[[URLdecode(x[1])]] <- URLdecode(x[2])
    } else {
      qs <- request$query.vector
      qn <- names(qs)
      for (i in seq.int(qs)) pars[[qn[i]]] <- qs[i]
    }
  
    # find the script
    request$path.info <- ''
    sfn <- sprintf("%s/web.R/%s.R", root, path)
    if (!file.exists(sfn)) { # if the file doesn't exist, we try to separate script name from PATH_INFO
      left <- path
      while (nzchar(left <- gsub("/[^/]*$", "", left)) && !file.exists(cand <- sprintf("%s/web.R/%s.R", root, left))) { if (!grepl("/", left)) left <- '' }
      if (!nzchar(left))
        return(c("html", paste("Script ", path, ".R not found", sep=''), "text/html", "Status: 404 Script not found"))
      request$path.info <- gsub(left, "", path, fixed=TRUE)
      sfn <- cand
    }

    oclear(TRUE, TRUE)
    .GlobalEnv$request <- request
    if(exists('init') && is.function(init)) init()

    source(sfn, local=TRUE)
    as.WebResult(do.call(run, pars))
  }, silent=TRUE))
}


## URLencode is *not* vectorized in R, believe it or not so we have to work around that ...
URLenc <- function(x) unlist(lapply(x, URLencode))

### this maps the Rhttpd/Rserve direct HTTP API into .run
.http.request <- function(url, query, body, headers) {
  root <- getOption("FastRWeb.root")
  if (is.null(root)) root <- "/var/FastRWeb"

  request <- list(uri=url, method='GET', c.type='', c.length=-1, body=NULL, client.ip='0.0.0.0', query.string='', raw.cookies='')
  if (length(query)) request$query.vector <- query ## back-door for parsed queries
  
  ## process headers to pull out request method (if supplied) and cookies
  if (is.raw(headers)) headers <- rawToChar(headers)
  if (is.character(headers)) {
    ## parse the headers into key/value pairs, collapsing multi-line values
    h.lines <- unlist(strsplit(gsub("[\r\n]+[ \t]+"," ", headers), "[\r\n]+"))
    h.keys <- tolower(gsub(":.*", "", h.lines))
    h.vals <- gsub("^[^:]*:[[:space:]]*", "", h.lines)
    names(h.vals) <- h.keys
    h.vals <- h.vals[grep("^[^:]+:", h.lines)]
    h.keys <- names(h.vals)

    if ("request-method" %in% h.keys) request$method <- c(h.vals["request-method"])
    if ("client-addr" %in% h.keys) request$client.ip <- c(h.vals["client-addr"])
    if ("cookie" %in% names(h.vals)) request$raw.cookies <- paste(h.vals[h.keys == "cookie"], collapse=" ")
  }

  ## this is a bit convoluted - the HTTP already parses the body - disable it where you can
  if (!is.raw(body)) {
    if (length(body)) {
      sb <- paste(unlist(lapply(names(body), function(x) paste(URLencode(x),"=",URLencode(as.character(body[[x]])),sep=''))),collapse='&')
      request$body <- charToRaw(sb)
      request$c.length <- length(request$body)
      request$c.type <- 'application/x-www-form-urlencoded'
    }
  } else {
    request$body <- body
    request$c.type <- attr(body, "content-type")
    request$c.length <- length(body)
  }

  r <- .run(request, root, url)
  if (length(r) < 2) return(list(r))
  cmd <- r[1]
  payload <- r[2]
  ct <- if (length(r) > 2) r[3] else "text/html"
  h <- if (length(r) > 3) r[4] else character(0)
  if (any(nchar(h) == 0L)) h <- h[nchar(h) > 0]
  if (cmd == "tmpfile" || cmd == "file") {
    fn <- paste(root, if (cmd == "tmpfile") "tmp" else "web", gsub("/", "_", payload, fixed=TRUE), sep='/')
    list(file=fn, ct, h)
  } else list(payload, ct, h)
}
