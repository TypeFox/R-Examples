# a. Correct use of if() without the else part
x <- 10
if (x > 5) { 
  m <- x + 200
}
x; m

# b. Correct use of if() with the else part
if (x > 5) {
  y <- x + 1
} else {
  y <- x + 2
}
x; y 

# c. Incorrect use of if/else; a mismatch of curly braces
if (x > 5) {
  z <- x + 1}
else {
  z <- x + 2 
} 
x; z 