# Function to get species IUCN species code
# Used inside lets.iucn functions
# Bruno Vilela

.getcode <- function(input) {
  
  binomialerror <- length(unlist(strsplit(input, " "))) == 2
  input <- gsub(as.matrix(input), pattern = " ", replacement = "-")
  
  h <- try(htmlParse(paste("http://api.iucnredlist.org/go/",
                           input, sep = "")),
           silent = TRUE)
  if (class(h)[1] != "try-error" & binomialerror) {
    b <- try(xpathSApply(h, '//div', xmlValue), silent = TRUE)[1]
    c <- as.numeric(gsub("\\D", "", b))
    
    
    # Subsecies control
    http <- "http://www.iucnredlist.org/details/summary/"
    h1 <- htmlParse(paste(http, c, "/0", sep = ""))
    links <- xpathSApply(h1, "//a/@href")
    links <- strsplit(links, "\n")
    parents <- xpathSApply(h1, "//a")
    
    # function to transform xml in character list
    tocharacter <- function(x) {
      do.call(paste, as.list(capture.output(x)))
    }
    parents2 <- sapply(parents, tocharacter)
    #menos <- length(parents2) - length(links)
    
    posParent <- grep("_parent", parents2)
    if (length(posParent) == 1) {
      #cpar <- gsub("\\D", "", (links[posParent - menos]))
      cpar <- gsub("\\D", "", (parents2[posParent]))
      c <- as.numeric(substr(cpar, 1, nchar(cpar) - 1))
    }
    ################################################
    return(c)
  } else {
    return(NULL)
  }
}