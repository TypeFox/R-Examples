# Indexing
x <- c(30, 100, 20); names(x) <- c("a", "b", "c")
x2 <- x[1]; x3 <- x[c("a", "c")]
y <- data.frame(aa = 50:52, bb = c(TRUE, FALSE, FALSE))
y2 <- y[c(1, 3), "bb"]; y3 <- y[2:3, ]; x3; y3

# Flow control: if(), for()
my.pound <- 100; diagno <- target <- NULL
if (my.pound > 80) {
  diagno <- "overweight"
  target <- my.pound * 0.9
}
diagno; target

my.pro <- 1; my.sum <- 0
for (i in 1:5) {
  my.pro <- my.pro * i
  my.sum <- my.sum + i
}
my.pro; my.sum 

# User-defined new functions
my.pound <- 100
my.fun <- function(x) {
  x * 0.454
}
my.kilo <- my.fun(x = my.pound); ls()
class(my.pound); my.pound
class(my.fun); my.fun
class(my.kilo); my.kilo