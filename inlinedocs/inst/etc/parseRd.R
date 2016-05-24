library(tools)
parsed <- parse_Rd("../man/package.skeleton.dx.Rd")
cleantags <- function(x)gsub("\\\\","",tools:::RdTags(x))
names(parsed) <- cleantags(parsed)
## all base-level TEXT tags are just whitespace...
stopifnot(all(parsed[names(parsed)=="TEXT"]=="\n"))
tl <- parsed[names(parsed)!="TEXT"]
arg <- tl[[which(names(tl)=="arguments")]]
names(arg) <- cleantags(arg)
arg <- arg[names(arg)!="TEXT"]
docs <-
  lapply(arg,function(L)paste(tools:::as.character.Rd(L[[2]]),collapse=""))
names(docs) <-
  sprintf("item{%s}",sapply(arg,function(L)tools:::as.character.Rd(L[[1]])))
docs <- c(docs,lapply(tl[names(tl)!="arguments"],paste,collapse=""))
library(inlinedocs)
ext <- extract.docs.file("../R/package.skeleton.dx.R")$package.skeleton.dx
## ignore definition and format
for(N in names(ext)){
  cat(N,":")
  if(N %in% names(docs)){
    if(ext[[N]]!=docs[[N]])cat("MISMATCH!!!")
  }else cat("does not exist in Rd!!!")
  cat("\n")
}
