# A. Multiple return statements
test1 <- function(x) {
  y <- x + 1
  if (x < 10) return(list(x = x, y = y))
  z <- x + 2
  if (x >= 10) return(list(x = x, z = z))
}
unlist(test1(x = 7)); unlist(test1(x = 15))

# B1. Use list to hold reults: define it at the end
test2 <- function(aa) {
  # Body: transformation
  bb <- aa + 83
  if (aa > 10) {cc <- aa * 3} 
  
  # Body: output 
  if (aa < 10) {
    result <- list(aa = aa, bb = bb)  # list at the end
  } else {
    result <- list(aa = aa, bb = bb, cc = cc)
  }    
  return(result)
}
unlist(test2(aa = 8)); unlist(test2(aa = 11)) 

# B1. Use list to hold reults: define it at the beginning
test3 <- function(aa) {
  # Body: transformation and output selection  
  result <- list(); result$aa <- aa
  bb <- aa + 83; result$bb <- bb
  if (aa > 10) {result$cc <- aa * 3}
  
  # Body: output 
  return(result)
}
unlist(test3(aa = 8)); unlist(test3(aa = 11))

# C1. Output: a value is returned and can be printed without saving
# The following three treatments are the same.
test4 <- function(x) {y <- x + 1; return(y)}; test4(10)
test5 <- function(x) {y <- x + 1; y}; test5(10)
test6 <- function(x) {x + 1}; test6(10)

# C2. Output: a value is returned, but cannot be printed directly
test7 <- function(x) {y <- x + 1} 
test8 <- function(x) {y <- x + 1; invisible(y)}
test9 <- function(x) {y <- x + 1; invisible(list(y, x))}
test7(10); test8(10); test9(10)  # No, cannot print the output 
(r7 <- test7(10)); (r8 <- test8(10)); (r9 <- test9(10))  # Yes, can print

# comparison: with or without invisible()
identical(r7, r8)  # TRUE 

# D. Functions without returning a value
plot.default; write.table
out <- plot(1:10); out  # out is NULL