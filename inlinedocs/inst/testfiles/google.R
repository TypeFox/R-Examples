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

google <- function(src,...){
  # Extract docs from google header comments.
  #
  # Args:
  # src: lines of code of the function source.
  #
  # Returns:
  # An inner Documentation List.
  i <- 2
  docs <- list()
  mode <- "description"
  comment.pattern <- "\\W*#\\W*"
  while(grepl(comment.pattern,line <- src[i])){
    content <- gsub(comment.pattern,"",line)
    if(grepl("Returns:",line)){
      mode <- "value"
    }else if(mode=="value"){
      docs$value <- c(docs$value,content)
    }
    if(grepl("Args:",line)){
      mode <- "args"
    }else if(mode=="args"){
      prefix <- gsub(":.*","",content)
      if(grepl(":",content))lname <- sprintf("item{%s}",prefix)
      suffix <- gsub(".*:\\s*","",content)
      docs[[lname]] <- c(docs[[lname]],suffix)
    }
    if(mode=="description"){
      docs$description <- c(docs$description,content)
    }
    i <- i+1
  }
  lapply(docs,function(x)gsub("\\s*$","",paste(x,collapse="\n")))
}
.parsers <- list(google=forfun(google))

#src <- getSource(CalculateSampleCovariance)
#google(src)
.dontcheck <- TRUE

.result <-
  list(google=list(description="Extract docs from google header comments.",
         `item{src}`="lines of code of the function source.",
         value="An inner Documentation List."),
       CalculateSampleCovariance=list(
         description="Computes the sample covariance between two vectors.",
         `item{x}`="One of two vectors whose sample covariance is to be calculated.",
         `item{y}`="The other vector. x and y must have the same length, greater than one,\nwith no missing values.",
         `item{verbose}`="If TRUE, prints sample covariance; if not, not. Default is TRUE.",
         value="The sample covariance between x and y."))
