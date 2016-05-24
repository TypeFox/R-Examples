library(XML)
htmlErrorHandler <- function(msg, code, domain, line, col, level, filename) {
  if (level > 2)
    stop("Failed to Parse HTML [", line, ":", col, "]: ", msg)
}

setwd("~/research/src/qt-x11-opensource-src-4.5.0/doc/html/")

parseTypes <- function(section) {
  files <- dir(pattern="^q[^3]")
  ## ind <- 1
  ## if (params)
  ##   ind <- "position() > 1"
  strs <- unlist(lapply(files, function(f) {
    dom <- htmlTreeParse(f, useInternalNodes=TRUE, error=htmlErrorHandler)
    path <- paste("/html//h3[text() = '", section,
                  "']/following-sibling::ul[1]/li/text()", sep = "")
    sapply(getNodeSet(dom, path), xmlValue)
  }))

  ## wee bit of cleanup
  strs <- gsub("[ )(,]", "", strs)
  strs <- strs[!grepl("[=\302]", strs)]
  strs <- strs[strs != "const" & nchar(strs)]
  strs <- sub("virtual", "", sub("const", "const ", strs))
  
  sort(table(strs), decreasing=TRUE)
}

countArgs <- function(section) {
  files <- dir(pattern="^q[^3]")
  strs <- unlist(lapply(files, function(f) {
    dom <- htmlTreeParse(f, useInternalNodes=TRUE, error=htmlErrorHandler)
    path <- paste("/html//h3[text() = '", section,
                  "']/following-sibling::ul[1]/li/text()",
                  sep = "")
    sapply(getNodeSet(dom, path), xmlValue)
  }))
  strs <- gsub("[ )(,]", "", strs)
  strs <- strs[!grepl("[=\302]", strs)]
  strs <- strs[strs != "const" & nchar(strs)]
  counts <- diff(grep("void", strs)) - 1
  counts
}

formatTable <- function(x) {
  cat(paste(" ", names(x), x, collapse="\n"))
}
  
