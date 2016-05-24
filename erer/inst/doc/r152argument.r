# A1. Too FEW arguments defined
# "x" is defined in the global environment. 
# Running "pa" in another session may not work if "x" is not defined.
x <- 10
pa <- function(y) {y + x} 
pa(y = 30)

# A2. Too MANY arguments defined
# "w3" is not used in the body at all; it works but can be confusing.
pb <- function(w1, w2, w3) {w1 + w2}
pb(w1 = 10, w2 = 40, w3 = 100)
pb(w1 = 10, w2 = 40)

# B. Typical formats for arguments
library(erer); data(daBed) 
pc <- function(y, c1 = NULL, c2 = FALSE, c3 = "increase", 
  c4 = c("year", "month"), c5 = c("Mon", "Wed", "Fri"), 
  c6 = 10 + nrow(y), ...) { 
  # Examine arguments without default value
  if (!inherits(y, "matrix")) {stop("y should be 'matrix').\n")}

  # Examine arguments with default value
  if (!is.null(c1)) {w.c1 <- 50} else {w.c1 <- 10}
  c5 <- match.arg(arg = c5, choices = c("Mon", "Wed", "Fri"))
  
  # Examine arguments with ellipsis
  c7 <- list(...)
  w.eps <- round(mean(y[, 3]), ...)
  
  result <- listn(y, c1, w.c1, c2, c3, c4, c5, c6, c7, w.eps)
  return(result)
}

bed <- daBed[1:3, 1:5]; bed 
ra <- pc(y = bed); ra
rb <- pc(y = bed, c3 = "any", c5 = "Wed", c6 = 43); rb
rc <- pc(y = bed, c3 = "all", c5 = "Fri", c6 = 43, digits = 3); rc
rc$y; rc$c1; rc$c4; rc$w.eps
rc$c7
unlist(rc[-c(1, 2, 6, 9, 10)])