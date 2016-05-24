#
# To get Keywords associated with a sequence.
#

getKeyword <- function(object, ...) UseMethod("getKeyword")

getKeyword.default <- function(object, ...)
  stop(paste("no getKeyword method for objects of class:", class(object)))

getKeyword.list <- function(object, ...)
  lapply(seq_len(length(object)), function(i) getKeyword(object[[i]], ...))

getKeyword.SeqAcnucWeb <- function(object, ..., socket = autosocket()){
  getKeywordsocket <- function(socket, name){
  #modif simon
  writeLines(paste("isenum&name=", name, sep = ""), socket, sep = "\n")
  res <- readLines(socket, n = 1)
  number <- parser.socket(res)[1] 

  writeLines(paste("readsub&num=", number, sep = ""), socket, sep = "\n")
  res2 <- readLines(socket, n = 1) 
  rr <- parser.socket(res2)

  writeLines(paste("readshrt&num=", rr[7], sep = ""), socket, sep = "\n")
  res3 <- readLines(socket, n = 1)
  #modif simon   

  # Get the nb of kw (not used here)
  # nbkws <- parser.socket(res3)[2]

  #recupere la liste de paires val, next 
  tmpl <- unlist(strsplit(res3, "&"))
  #transforme en liste
  tmpl <- unlist(strsplit(tmpl[3],","))
  kwl <- unlist(tmpl)[c(TRUE, FALSE)]

  lapply(kwl, function(x){
    writeLines(paste("readkey&num=", x, sep = ""), socket, sep = "\n")  
    res4 <- readLines(socket, n = 1)
    res <-parser.socket(res4)[2]
    substring(res[1], 2, nchar(res[1]) - 1)
  })

} 
  unlist(getKeywordsocket(socket, name = object))
}

getKeyword.qaw <- function(object, ...) getKeyword(object$req, ...)

getKeyword.logical <- function (object, ...)
  object # so that NA is returned for virtual lists
