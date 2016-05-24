
interpret_ds <- function(ds, .curl=dmCurlHandle()) {
  base <- api_base
  if (class(ds) == 'list' && setequal(
      names(ds),
      c('status', 'dimensions', 'meta', 'title', 'ds', 'id')
      )) {
    infos <- list(ds)
    names(infos) <- c(ds$id)
    return(list(base=base, qs=list(ds=ds$id), infos=infos))
  } else if (class(ds) == 'list' && length(ds) >= 1 && setequal(
      names(ds[[1]]),
      c('status', 'dimensions', 'meta', 'title', 'ds', 'id')
      )) {
    infos <- ds
    names(infos) <- sapply(ds, FUN=function(s) s$id)
    return(list(
      base=base,
      qs=list(ds=paste(sapply(ds, FUN=function(s) s$id), collapse='/')),
      infos=infos
      ))
  }
  if (grepl('^https?:', ds)) {
    spliturl <- urlsplit(ds)
    base <- spliturl$base
    path <- spliturl$path
    anchor <- parse_qs(spliturl$anchor)
    qs <- parse_qs(spliturl$qs)
    if ('ds' %in% names(anchor)) {
      qs <- anchor
    } else if ('ds' %in% names(qs)) {
    } else if (base %in% short_url_services) {
      h=basicTextGatherer()
      getURL(ds, headerFunction=h$update, curl=.curl)
      conn <- textConnection(paste(h$value(NULL)[-1], collapse=""))
      headers=as.list(read.dcf(conn)[1,])
      close(conn)
      if ('Location' %in% names(headers)) {
        return(interpret_ds(headers$Location, .curl=.curl))
      } else {
        stop("No redirect found for URL ", ds)
      }
    } else if (grepl('(?:/[a-zA-Z]{2})?/data/set/(?:[0-9a-z]+)(?:/|$)', path)) {
      qs <- list(ds=sub('.*/data/set/([0-9a-z]+).*', '\\1', path))
    } else {
      stop("Can't make sense of URL ", ds)
    }
  } else if (grepl('&|^ds=', ds)) {
    qs <- parse_qs(ds)
  } else {
    qs <- list(ds=ds)
  }
  list(base=base, qs=qs)
}

dimfilter <- function(ds, infos, ...) {
  args <- list(...)
  newds <- c()
  if ("dmdataset" %in% class(infos)) {
    infosds <- infos$ds;
    infos <- list(infos);
    names(infos) <- list(infosds);
  }
  for (dsid in strsplit(ds, '\\/')[[1]]) {
    dimspec <- c()
    if (dsid %in% names(infos)) {
      for (name in names(args)) {
        matchdim <- NULL
        for (dim in infos[[dsid]]$dimensions) {
          if (dim$title == name) {
            matchdim <- dim
            break
          }
        }
        if (is.null(matchdim) && name %in% names(infos[[dsid]]$dimensions)) {
          matchdim <- infos[[dsid]]$dimensions[[name]]
        }
        if (!(is.null(matchdim))) {
          valueids <- c()
          for (value in args[[name]]) {
            valueid <- NULL
            for (prospect in dim$values) {
              if (prospect[['title']] == value) {
                valueid <- prospect[['id']]
                break
              }
            }
            if (is.null(valueid) && value %in% names(dim$values)) {
              valueid <- value
            }
            if (!(is.null(valueid))) {
              valueids <- c(valueids, valueid)
            }
          }
          if (identical(valueids, c())) {
            stop(paste("No match found for '", name, "'='", args[name], "'",
                       sep=""))
          }
          dimspec <- c(dimspec, paste(
            matchdim$id, '=', paste(valueids, collapse='.'),
            sep=''
            ))
        }
      }
    }
    newds <- c(newds,
      ifelse(
        identical(dimspec, c()),
        dsid,
        paste(dsid, '!', paste(dimspec, collapse=':'), sep='')
      ))
  }
  return(paste(newds, collapse='/'))
}

urlsplit <- function(url) {
  uri_and_query_and_anchor <- strsplit(url, '\\#')[[1]]
  uri_and_query <- strsplit(uri_and_query_and_anchor[[1]], '\\?')[[1]]
  uri <- uri_and_query[[1]]
  slashinds <- gregexpr('/', uri)[[1]]
  if (length(slashinds) >= 3 &&
      (slashinds[1:2]==c(6,7) || slashinds[1:2] == c(7,8)) &&
      substring(uri, 1, 4) == 'http'
    ) {
    slashind <- slashinds[3]
  } else {
    slashind <- slashinds[1]
  }
  list(
    base=substring(uri, 1, slashind - 1),
    path=substring(uri, slashind),
    qs=uri_and_query[2],
    anchor=uri_and_query_and_anchor[2]
  )
}

parse_qs <- function(qs) {
  if (class(qs) == 'list') {
    return(qs)
  }
  if (!is.na(qs) && grepl('^!', qs)) {
    qs <- substr(qs, 2, 100000)
  }
  l <- lapply(as.list(strsplit(qs, '&', fixed=TRUE)[[1]]), FUN=function(pair) {
    keyval <- as.list(strsplit(pair, '=', fixed=TRUE)[[1]])
    if (is.na(keyval[[1]]) || keyval[[1]] %in% c("display", "s", "e", "f")) {
      return (NA)
    }
    val <- ifelse(length(keyval) > 1, paste(keyval[-1], collapse='='), '')
    names(val) <- keyval[[1]]
    return(val)
  })
  l <- Filter(function(x) !is.na(x), l)
  names(l) <- sapply(l, names)
  return(l)
}

get.datamarket.csv <- function(ctx, path, curl, .params, origds=NA) {
  curlopts = list()
  if (!is.na(origds)) {
    curlopts$Referer <- origds
  }
  content <- getForm(
    paste(ctx$base, path, sep=""),
    curl=curl,
    .opts=curlopts,
    .params=c(ctx$qs, split_time=0, callback="", .params)
    )
  if (is.raw(content)) {
    content <- rawToChar(content)
  }
  conn <- textConnection(content)
  csv <- read.csv(conn, header=TRUE)
  close(conn)
  if (names(csv)[[1]] == 'request.error') {
    stop(paste("Request error: ", csv[[1]]))
  }
  if (names(csv)[[1]] == 'server.error') {
    stop(paste("Server error:", csv[[1]]))
  }
  return(csv)
}

