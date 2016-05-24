# A. Minimum information for a function
test <- function() {3}; test()
(function() {})()

# B. Function properties
dog <- function(x) {
  y <- x + 10
  z <- x * 3
  w <- paste("x times 3 is equal to", z, sep = " ")
  result <- list(x = x, y = y, z = z, w = w)
  return(result)
}
class(dog); args(dog); formals(dog); body(dog)
dog(8)              # default printing for all
res <- dog(8); res  # assignment and selected printing
res$x; res$w  

# C. Anonymous function
ga <- lapply(X = 1:3, FUN = seq); ga
my.seq <- function (y) {seq(from = 1, to = y, by = 1)}
gb <- lapply(X = 1:3, FUN = my.seq)
gc <- lapply(X = 1:3, FUN = function(y) {seq(from = 1, to = y, by = 1)})
gc
identical(ga, gb); identical(gb, gc)