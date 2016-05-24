simple <- function(src,...){#title a simple Parser Function
  noquotes <- gsub("([\"'`]).*\\1","",src)
  print(noquotes)
  comments <- grep("#",noquotes,value=TRUE)
  pat <- "[^#]*#([^ ]*) (.*)"
  tags <- gsub(pat,"\\1",comments)
  docs <- as.list(gsub(pat,"\\2",comments))
  names(docs) <- tags
  docs[tags!=""]#value all the tags with a single pound sign
}

.parsers <- list(simple=forfun(simple))

testfun <- function(x,y,z){#item{x} the first arg
  a <- (x+y)*z # just a regular comment
  a #value the sum of the first two times the third
  #description a useless formula
}

.result <- list(testfun=list(`item{x}`="the first arg",
                  value="the sum of the first two times the third",
                  description="a useless formula"),
                simple=list(title="a simple Parser Function",
                  value="all the tags with a single pound sign"))

.dontcheck <- TRUE # because of use of forfun here
