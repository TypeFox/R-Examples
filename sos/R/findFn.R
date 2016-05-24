findFn <- function(string,
                   maxPages = 20,
                   sortby = NULL,
                   verbose = 1, ...) {
  ####################################################################
  ##
  ## RSiteSearch(string, "fun")
  ##
  ## returning a data.frame with one row for each hit
  ## giving package and "function" / entry name,
  ## sorted with the most frequently selected package first
  ##
  ## with the following attributes
  ## matches = total number of hits found by the search
  ## summary = sort(table(ans$Package));
  ##      summary contains results for only the first maxPages,
  ##      so sum(summary) may be less than hits.
  ##
  ##
  ##
  ## 1.  Define internal function
  ## internal functions
  parseLinks <- function(links) {
    lnk <- sub("<dt>.*<strong><a href=\\\"(.*\\.html)\\\">.*$", "\\1",
               links, useBytes = TRUE)
    desc0 <- sub("<dt>.*<strong><a href=\\\".*\\\">R.*:(.*)</a>.*$", "\\1",
                 links, useBytes = TRUE)
    desc <- gsub("(<strong class=\"keyword\">)|(</strong>)|^[ ]+|[ ]+$", "",
                 desc0, useBytes = TRUE)
    list(link = lnk, description = desc)
  }
  parseHTML <- function(href) {
    link <- try(url(href))
    on.exit(close(link))
    if (inherits(link, "try-error")) {
      warning("An error occurred opening ", href,
              "\nfindFn needs Internet access;  is it available?")
      ch0 <- character(0)
      ans <- data.frame(Package = ch0,
                        Function = ch0,
                        Date = ch0,
                        Score = numeric(0),
                        Description = ch0,
                        Link = ch0, stringsAsFactors=FALSE)
      attr(ans, "matches") <- 0
      return(ans)
    }
    html <- try(readLines(link))
    if (inherits(html, "try-error")) {
      warning("An error occurred in readLine(link), link = ", link,
              "\nfindFn needs Internet access;  is it available?")
      ch0 <- character(0)
      ans <- data.frame(Package = ch0,
                        Function = ch0,
                        Date = ch0,
                        Score = numeric(0),
                        Description = ch0,
                        Link = ch0, stringsAsFactors=FALSE)
      attr(ans, "matches") <- 0
      return(ans)
    }
#   Find the hit count
    hitPattern <- "^.*<!-- HIT -->(.*)<!-- HIT -->.*$"
    hitRows <- html[grep(hitPattern, html, useBytes = TRUE)]
#
    Hits <- as.numeric(sub(hitPattern, "\\1", hitRows, useBytes = TRUE))
#   Find dates
    findDates <- grep("<strong>Date</strong>", html, useBytes = TRUE)
    dates <- html[findDates]
    links <- html[findDates - 2]
    pattern <-
      "^.*http://finzi.psych.upenn.edu/R/library/(.*)/html/(.*)\\.html.*$"
    pac <- sub(pattern, "\\1", links, useBytes = TRUE)
    fun <- sub(pattern, "\\2", links, useBytes = TRUE)
    scoreLoc <- regexpr("score:", links, useBytes = TRUE)
    scorePattern <- "^.*\\(score: ([0-9]+)\\).*$"
    scoreCh <- sub(scorePattern, "\\1", links, useBytes = TRUE)
    score <- as.numeric(scoreCh)
    Date <- sub("^.*<em>(.*)</em>.*$", "\\1", dates, useBytes = TRUE)
    pLinks <- parseLinks(links)
#    if (length(pac) < 1 && length(Date) > 0) {
    if (length(pac) < 1) {
      countDocs <- grep("Too many documents hit. Ignored",
                        html, useBytes = TRUE)
#      tooMany <- length(tooMany) > 0
      tooMany <- (length(countDocs)>0)
      if (tooMany) {
        Hits <- Inf
        warning("Too many documents hit.  Ignored")
      } else {
        if(length(Date)>0){
            op <- 'SOFTWARE PROBLEM:  dates found without'
            oops <- paste(op, 'content;  ignored.')
            warning(oops)
        }
      }
      Date <- Date[numeric(0)]
    }
    Ans <- data.frame(Package = pac,
                      Function = fun,
                      Date = strptime(Date, "%a, %d %b %Y %T"),
                      Score = score,
                      Description = pLinks$description,
                      Link = pLinks$link, stringsAsFactors=FALSE)
    oops <- (substring(pac, 1, 1)=='<')
    ans <- Ans[!oops,]
    attr(ans, "matches") <- Hits
    ans
  }
  ##
  ## end internal functions
  ##
  quiet <- (verbose < 2)
  ##
  ## 2.  Set up query
  ##
  if (substr(string, 1, 1) != "{") {
    string <- gsub(" ", "+", string)
  } else {
    ## scan(url(...)) fails with spaces
    string <- gsub(" ", "%20", string)
  }
  fmt <- paste("http://search.r-project.org/cgi-bin/namazu.cgi?",
          "query=%s&max=20&result=normal&sort=score&idxname=functions",
               sep = "")
  href <- sprintf(fmt, string)
  ##
  ## 3.  Query
  ##
  ##  3.1.  Set up
  ans <- parseHTML(href)
  hits. <- attr(ans, "matches")
  if (length(hits.) < 1) {
    warning("HIT not found in HTML;  processing one page only.")
    hits. <- nrow(ans)
    attr(ans, "matches") <- hits.
  } else {
    if (length(hits.) > 1) {
      warning("HIT found more than once in first HTML page; ",
              "first 2 = ", hits.[1], ", ", hits.[2],
              ";  processing one page only")
      hits. <- nrow(ans)
      attr(ans, "matches") <- hits.
    }
  }
  if (is.na(hits.)) {
    warning("HIT found, not numeric, in the first HTML page; ",
            "processing one page only")
    hits. <- nrow(ans)
    attr(ans, "matches") <- hits.
  }
  if (verbose) {
    es <- if (hits. == 1) "" else "es"
    cat("found ", hits., " match", es, sep = "")
  }
  ##  3.2.  Retrieve
  n <- min(ceiling(hits./20), maxPages)
  if((hits.<Inf) && (nrow(ans) < hits.)) {
    if (verbose)
      cat(";  retrieving", n, c("page", "pages")[1 + (n > 1)])
    if (verbose) {
      if ((20 * n) < hits.) {
        cat(",", 20 * n, "matches.\n")
      } else {
        cat("\n")
      }
    }
    for(i in seq(2, length = n - 1)) {
      if (!quiet) {
        cat("retrieving page ", i, " of ", n, "\n", sep = "")
      } else if (verbose > 0) {
        cat(i, "")
        if((i%%10)==0) cat('\n')
        flush.console()
      }
      href.i <- sprintf("%s&whence=%d", href, 20 * (i - 1))
      ans <- rbind(ans, parseHTML(href.i))
    }
    if (verbose>0) cat("\n")
  } else {
    cat("\n")
  }
  ##
  ## 4.  Compute Summary
  ##
  ans$Score <- as.numeric(as.character(ans$Score))
  pkgSum <- PackageSummary(ans)
  if((hits.<Inf) && (hits.>0)){
      nlk <- sum(pkgSum$Count)
      cat('Downloaded', nlk, 'links in', nrow(pkgSum), 'packages.\n')
  }
  ##
  ## 5.  Sort order
  ##
  s0 <- c("Count", "MaxScore", "TotalScore", "Package",
          "Score", "Function", "Date", "Description", "Link")
  s0. <- tolower(s0)
  if (is.null(sortby)) {
    sortby <-  s0
  } else {
    s1 <- match.arg(tolower(sortby), s0., TRUE)
    s1. <- c(s1, s0.[!(s0. %in% s1)])
    names(s0) <- s0.
    sortby <- s0[s1.]
  }
  ##
  ## 6.  Merge(packageSum, ans)
  ##
  packageSum <- pkgSum
  rownames(pkgSum) <- as.character(pkgSum$Package)
  pkgSum$Package <- NULL
  pkgSum$Date <- NULL
  pkgS2 <- pkgSum[as.character(ans$Package), , drop = FALSE]
  rownames(pkgS2) <- NULL
  Ans <- cbind(as.data.frame(pkgS2), ans)
  ##
  ## 7.  Sort Ans by "sort."
  ##
  Ans.num <- Ans[, c("Count", "MaxScore", "TotalScore", "Score")]
  ans.num <- cbind(as.matrix(Ans.num), Date=as.numeric(Ans$Date) )
  Ans.ch <- Ans[, c("Package","Function", "Description", "Link")]
  ans.ch <- as.data.frame(as.matrix(Ans.ch))
  ansKey <- cbind(as.data.frame(-ans.num), ans.ch)
  ##
  oSch <- do.call("order", ansKey[sortby])
  AnSort <- Ans[oSch, ]
  ##
  ## 8.  attributes
  ##
  rownames(AnSort) <- NULL
  ##
  attr(AnSort, "matches") <- hits.
  attr(AnSort, "PackageSummary") <- packageSum
  attr(AnSort, "string") <- string
  attr(AnSort, "call") <- match.call()
  class(AnSort) <- c("findFn", "data.frame")
  AnSort
}
