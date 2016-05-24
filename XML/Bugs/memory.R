library(XML)
f = system.file("exampleData", "mtcars.xml", package="XML")
curxmltext <- paste(readLines(f),collapse="\n")

memoryLeakDoc <- function() {
  require(XML)
  for(i in 1:10000) {
    gg <- xmlParse(curxmltext,
                   asText=TRUE)
    rm(gg)
  }
}

memoryLeak <- function(file = f, useText = FALSE, n = 100) {
  require(XML)
  if(useText)
    f = readLines(file)
  
  for(i in 1:n) {
    gg <- xmlParse(f, asText = useText)
    ret <- xpathApply(doc = gg,
                      path = "/blah",
                      fun = xmlValue)

    rm(ret)
#    free(gg)
    rm(gg)
  }
  gc(); gc()
}


memoryLeak3 <- function(file = f, useText = FALSE, n = 100) {
  require(XML)
  if(useText)
    f = readLines(file)

  ans = vector("list", n)  
  for(i in 1:n) {
    gg <- xmlParse(f, asText = useText)
    ret <- xpathApply(doc = gg, path = "//record[1]")
    ans[[i]] = ret
    rm(ret)
    rm(gg)
    gc()
  }
  gc(); gc()
  ans
}




require(XML)
memoryLeak2 <- function() {

  for(i in 1:10000) {
    gg <- xmlParse(curxmltext, asText=TRUE)
    xmlRoot(gg)
    rm(gg)
  }
  gc(); gc()
}


if(FALSE) {
 library(XML)
 f = system.file("exampleData", "mtcars.xml", package="XML")
 doc = xmlParse(f)
 invisible(getNodeSet(doc, "/blah"))
 rm(doc)
 gc()
}

if(FALSE) {
  library(XML)
  f = system.file("exampleData", "mtcars.xml", package="XML")
  doc = xmlParse(f)
  invisible(xmlRoot(doc))
  gc(); gc()
  rm(doc)
  gc()
}
