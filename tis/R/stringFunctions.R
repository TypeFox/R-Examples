blanks <- function(n) sprintf(paste("%", round(n), "s", sep = ""), "")

stripBlanks <- function(strings){
  ## strips leading and trailing blanks
  ans <- gsub("^ *(.*)", "\\1", strings)
  notEmpty <- nchar(ans) > 0
  ans[notEmpty] <- gsub(" *$", "", ans[notEmpty])
  ans
}

