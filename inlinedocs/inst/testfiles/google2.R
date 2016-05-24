CalculateSampleCovariance <- function(x, y, verbose = TRUE) {
  # Computes the sample covariance between two vectors.
  #
  # Args:
  #   x: One of two vectors whose sample covariance is to be calculated.
  #   y: The other vector. x and y must have the same length, greater than one,
  #      with no missing values.
  #   verbose: If TRUE, prints sample covariance; if not, not. Default is TRUE.
  #
  # Returns:
  #   The sample covariance between x and y.
  n <- length(x)
  # Error handling
  if (n <= 1 || n != length(y)) {
    stop("Arguments x and y have invalid lengths: ",
         length(x), " and ", length(y), ".")
  }
  if (TRUE %in% is.na(x) || TRUE %in% is.na(y)) {
    stop(" Arguments x and y must not have missing values.")
  }
  covariance <- var(x, y)
  if (verbose)
    cat("Covariance = ", round(covariance, 4), ".\n", sep = "")
  return(covariance)
}

str_match_perl <- function(string,pattern){
  parsed <- regexpr(pattern,string,perl=TRUE)
  names <- attr(parsed,"capture.names")
  captured.text <- substr(string,parsed,parsed+attr(parsed,"match.length")-1)
  captured.text[captured.text==""] <- NA
  captured.groups <- do.call(rbind,lapply(seq_along(string),function(i){
    if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(names)))
    st <- attr(parsed,"capture.start")[i,]
    substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
  }))
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",names)
  result
}

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


google <- function(src,...){
  # Extract docs from google header comments.
  #
  # Args:
  # src: lines of code of the function source.
  #
  # Returns:
  # An inner Documentation List.
  pre <- "\\W*#\\W+"
  lines.pat <- paste("(?:",pre,".*?\n)+",sep="")
  google.doc.pattern <-
    paste("(?<description>",lines.pat,")",
          pre,"Args:\\W*\n",
          "(?<args>",lines.pat,")",
          pre,"Returns:\\W*\n",
          "(?<value>",lines.pat,")",
          sep="")
  lines <- paste(src,collapse="\n")
  parsed <- str_match_perl(lines,google.doc.pattern)
  if(is.na(parsed[1]))return(list())
  docs <- list(description=gsub(pre,"",parsed[,"description"]),
               value=gsub(pre,"",parsed[,"value"]))
  arg.pat <-
    paste(pre,"(?<name>[a-z]+): ",
          "(?<content>.*?\n",
          "(?:",pre,"[^:]*?\n)*)",
          sep="")
  m <- str_match_all_perl(parsed[,"args"],arg.pat)[[1]]
  argsList <- as.list(gsub(paste("\n",pre,sep=""),"\n",m[,"content"]))
  names(argsList) <- sprintf("item{%s}",m[,"name"])
  lapply(c(docs,argsList),function(x)gsub("\n$","",x))
}

.parsers <- list(google=forfun(google))

##src <- getSource(CalculateSampleCovariance)
## src <- getSource(google)
## src <- getSource(str_match_all_perl)
## google(src)

.dontcheck <- TRUE

.result <-
  list(google=list(description="Extract docs from google header comments.",
         `item{src}`="lines of code of the function source.",
         value="An inner Documentation List."),
       str_match_all_perl=list(),
       str_match_perl=list(),
       CalculateSampleCovariance=list(
         description="Computes the sample covariance between two vectors.",
         `item{x}`="One of two vectors whose sample covariance is to be calculated.",
         `item{y}`="The other vector. x and y must have the same length, greater than one,\nwith no missing values.",
         `item{verbose}`="If TRUE, prints sample covariance; if not, not. Default is TRUE.",
         value="The sample covariance between x and y."))
