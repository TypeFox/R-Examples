`autoFill` <-
function(x, squash=FALSE)
{
  ## description: fill out blanks of a vector with preceeding label
  ## eg. c("a","","","b","") becomes c("a","a","a","b","b")

  ## Arguments:
  ## x: data frame of two character vectors
  ## squash: remove spaces from labels (default: TRUE)
  ## NB: NO ERROR CHECKING SO BE CAREFUL

  ## require("gdata")  # obsolete in current R as in Depends/Imports
  
  ## set NA's to "" otherwise a problem
  if (is.factor(x))
    levels(x) <- c(levels(x),"")
  x[is.na(x)] <- ""
  ## and of course - what about leading/trailing spaces or just spaces...
  if (is.factor(x)) {
    levels(x) <- gdata::trim(levels(x))
    } else {
      x <- gdata::trim(x)
    }
  
  if (x[1] == "") cat("Warning: first value of vector is blank!\n")

  for (i in 2:length(x)) {
    if (x[i] == "")
      x[i] <- x[i-1]
  }

  if (squash) {
    if (is.factor(x)) {
      levels(x) <- gsub(" ","",levels(x))
      x <- factor(x)
    }
    if (is.character(x))
      x <- gsub(" ","",x)
  }
  x
}

