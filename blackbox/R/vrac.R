"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
"%innc%" <- function(x, y) (tolower(x) %in% tolower(y)) # a nocase version of %in%
"%==nc%" <- function(x, y) (tolower(x) == tolower(y)) # a nocase version of ==
