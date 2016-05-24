.createDepthStructure <- function(l){
  # Transform a nested list to a vector of objects with two elements: $depth and $str
  # $depth refers to the nested level
  # $str refers to the string representation for that element in the nested list
  s <- list()
  for(i in 1:length(l))
  {
    if(is.list(l[[i]])){
      recX <- .createDepthStructure(l[[i]])
      recX <- lapply(recX,function(x){
        x$depth <- x$depth+1
        x
      })
      s <- do.call(c,list(s,recX))
    }
    else{
      x <- list(depth=0, str = l[[i]])
      s <- c(s,list(x))
    }
  }
  s
}
.nestedList2String <- function(l, numbered=TRUE){
  # A toString method por nested lists.
  x <- .createDepthStructure(l)
  s <- ""
  for(i in 1:length(x)){
    sep <- ifelse(numbered, paste0(i," ) "), '*) ')
    s <- paste(s,do.call(function(...){ paste(...,sep="") },as.list(c(rep('\t',x[[i]]$depth),sep,x[[i]]$str,'\n'))),sep="")
  }
  s
}

.nestedList2HTML <- function(l, numbered=TRUE){
  # A toString method por nested lists.
  x <- .createDepthStructure(l)
  depth <- 0
  sepA <- ifelse(numbered, "<ol>\n", "<ul>\n")
  sepB <- ifelse(numbered, "\n</ol>", "\n</ul>")
  s <- sepA
  for (i in 1:length(x)){
    
    if (x[[i]]$depth > depth)
      s <- paste0(s, sepA)
    if (x[[i]]$depth < depth)
      s <- paste0(s, sepB)
    
    depth <- x[[i]]$depth
    s <- paste0(s, sprintf("<li>%s</li>\n", x[[i]]$str))
  }
  s <- paste0(s, sepB)
  s
}

.nestedList2Latex <- function(l){
  # A toString method por nested lists.
  x <- .createDepthStructure(l)
  depth <- 0
  s <- "\\begin{enumerate}\n"
  for (i in 1:length(x)){
    
    if (x[[i]]$depth > depth)
      s <- paste0(s, "\\begin{enumerate}\n")
    if (x[[i]]$depth < depth)
      s <- paste0(s, "\n\\end{enumerate}")
    
    depth <- x[[i]]$depth
    s <- paste0(s, sprintf("\\item %s\n", x[[i]]$str))
  }
  s <- paste0(s, "\n\\end{enumerate}")
  s
}