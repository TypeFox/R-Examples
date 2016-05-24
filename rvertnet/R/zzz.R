# Wrapper for search functions vertsearch, searchbyterm, spatialsearch and bigsearch

vertwrapper <- function(fxn = "", args = NULL, lim = NULL, rfile = NULL, email = NULL,
                        compact = TRUE, verbose = TRUE, ...){
  mssg(verbose, "Processing request...")
  if (fxn == "bigsearch") {
    tt <- GET(vdurl(), query = list(q = make_bigq(args, email, rfile)), ...)
    stop_for_status(tt)
    mssg(verbose, "\nThank you! Download instructions will be sent by email.")
  } else {
    ress <- vert_GET(fxn, args, lim, verbose, ...)
    
    # Remove columns populated fully by NAs
    if (compact) { 
      r <- ress$data
      ress$data <- r[,!sapply(r, function(x) all(is.na(x)))]
    }
    
    # Return results
    if (NROW(ress$data) == 0) {
      NULL
      mssg(verbose, "No records match this search request")
      if (fxn == "spatialsearch") { 
        mssg(verbose, "Check signs on decimal longitude and latitude") 
      }
    } else {
      ress
    }
  }
}

mssg <- function(v, ...) if (v) message(...)

get_terms <- function(){
  url <- "https://raw.githubusercontent.com/tdwg/dwc/master/downloads/SimpleDwCTermsList.txt"
  termlist <- read.table(text = content(GET(url), as = "text", encoding = "UTF-8"), stringsAsFactors = FALSE)
  # Strip embedded header from termlist and deal with upper vs. lower case in termlist vs. out$recs
  if (grep("term", tolower(termlist[1,1]))) termlist <- as.data.frame(termlist[-1,1], stringsAsFactors = FALSE)
  fullr <- as.data.frame(matrix(NA, 1, length(termlist[,1]))) # Create a full data frame to populate
  colnames(fullr) <- tolower(termlist[,1]) # Sync case to facilitate merge
  list(termlist = termlist, fullr = fullr)
}

vert_GET <- function(fxn="searchbyterm", args, limit = 1000, verbose = TRUE, ...){
  cursor <- NULL
  allres <- 0
  result <- list()
  i <- 0
  while (allres < limit) {
    i <- i + 1
    tt <- GET(vurl(), query = list(q = make_q(fxn, args, cursor, getlim(limit, allres))), ...)
    stop_for_status(tt)
    txt <- content(tt, "text", encoding = "UTF-8")
    out <- jsonlite::fromJSON(txt)
    avail <- out$matching_records
    cursor <- out$cursor
    result[[i]] <- out$recs
    allres <- sum(vapply(result, NROW, 1))
    if (char2num(avail) <= allres) allres <- limit
  }
  df <- if (sum(sapply(result, NROW)) == 0) data.frame(NULL, stringsAsFactors = FALSE) else rbind_all(result)
  names(df) <- tolower(names(df))
  res <- get_terms()
  df <- merge(res$fullr, df, all = TRUE)[, tolower(res$termlist[,1]) ]
  df <- df[ -NROW(df), ]
  mssg(verbose, paste("\nLast Query URL: \"", tt$url, "\"", sep = ""))
  mssg(verbose, paste("\nMatching records:", NROW(df), "returned,", avail, "available", sep = " "))
  list(meta = make_meta(out), data = tbl_df(df))
}

make_q <- function(fxn, x, cursor = NULL, limit=1000){
  if (fxn == "vertsearch") x <- paste0(unname(unlist(x)), collapse = " ")
  if (fxn == "spatialsearch") x <- sprintf("distance(location,geopoint(%s,%s))<%s", x$lat, x$long, x$radius)
  if (!is.null(limit)) {
    if (!is.null(cursor)) {
      ff <- sprintf('{"q":"%s","l":%s,"c":"%s"}', noc(gsub('\"|\\{|\\}', "", jsonlite::toJSON(x, auto_unbox = TRUE)), fxn), limit, cursor)
    } else {
      ff <- sprintf('{"q":"%s","l":%s}', noc(gsub('\"|\\{|\\}', "", jsonlite::toJSON(x, auto_unbox = TRUE)), fxn), limit)
    }
  } else {
    ff <- sprintf('{"q":"%s"}', noc(gsub('\"|\\{|\\}', "", jsonlite::toJSON(x, auto_unbox = TRUE)), fxn))
  }
  tmp <- gsub(":>", ">", gsub(":<", "<", gsub(":=", "=", ff)))
  gsub("year\\.[0-9]", "year", tmp)
}

vurl <- function() "http://api.vertnet-portal.appspot.com/api/search"
vdurl <- function() "http://api.vertnet-portal.appspot.com/api/download"

make_meta <- function(x){
  tmp <- x[ !names(x) %in% "recs" ]
  plyr::rename(tmp, c(cursor = "last_cursor"))
}

getlim <- function(x, y){
  stopifnot(is.numeric(x), is.numeric(y))
  toget <- x - y
  if (toget < 0) 0 else toget
}

char2num <- function(x) as.numeric(strextract(x, "[0-9]+"))
strextract <- function(str, pattern) regmatches(str, regexpr(pattern, str))
strtrim <- function(str) gsub("^\\s+|\\s+$", "", str)

noc <- function(x, fxn){
  if (fxn == "spatialsearch") {
    x
  } else {
    gsub(",", " ", x)
  }
}

make_bigq <- function(x, email, rfile){
  ff <- sprintf('{"q":"%s","n":"%s","e":"%s"}', noc(gsub('\"|\\{|\\}', "", jsonlite::toJSON(x, auto_unbox = TRUE)), ""), rfile, email)
  tmp <- gsub(":>", ">", gsub(":<", "<", gsub(":=", "=", ff)))
  gsub("year\\.[0-9]", "year", tmp)
}

combyr <- function(x) {
  if (!is.null(x) && length(x) > 1) {
    setNames(as.list(x), rep("year", length(x)))
  } else {
    list(year = x)
  } 
}

checkfourpkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    stop("Please install ", x, call. = FALSE)
  } else {
    invisible(TRUE)
  }
}
  