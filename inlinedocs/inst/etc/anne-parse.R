library(inlinedocs)
str_match_all_perl <- function(string,pattern){
  parsed <- gregexpr(pattern,string,perl=TRUE)
  lapply(seq_along(parsed),function(i){
    r <- parsed[[i]]
    starts <- attr(r,"capture.start")
    if(r[1]==-1)return(matrix(nrow=0,ncol=1+ncol(starts)))
    names <- attr(r,"capture.names")
    lengths <- attr(r,"capture.length")
    full <- substring(string[i],r,r+attr(r,"match.length")-1)
    subs <- substring(string[i],starts,starts+lengths-1)
    m <- matrix(c(full,subs),ncol=length(names)+1)
    colnames(m) <- c("",names)
    m
  })
}

anne.setMethod <- function(code,...){
  re <- paste('#\\s*(?<description>.*?)\\s*',
              'set(?<replace>Replace)?Method[(]\\s*',
              'f\\s*=\\s*"(?<f>.*?)"\\s*,',
              '\\s*signature\\s*=\\s*"(?<signature>.*?)"',
              sep="")
  txt <- paste(code,collapse="\n")
  extracted <- str_match_all_perl(txt,re)[[1]][,-1]
  rownames(extracted) <-
    ifelse(extracted[,"replace"]=="Replace","method","replace")
  print(extracted)
  list("A_--methods"=as.list(extracted["replace",]),
       "A-methods"=as.list(extracted["method",]))
}

custom.parsers <- list(anne.setMethod=anne.setMethod)

## we assume anne.R is in the current directory, with some setMethod
## definitions.
extract.docs.file("anne.R",custom.parsers)

# then run inlinedocs on your package like this
package.skeleton.dx("annepkg",c(custom.parsers))
