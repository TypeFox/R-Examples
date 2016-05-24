# Loops A, B, C, and D generate very similar results.
# A. for loop
x1 <- 2; x2 <- NULL
for (i in 5:7) {
  x2 <- c(x2, x1 * i)
}
x2; i

# B. while loop with "cond" being evaluated
y1 <- 2; y2 <- NULL; j <- 5 
while (j < 8) {
  y2 <- c(y2, y1 * j)
  j <- j + 1
}
y2; j

# C. while loop with break
z1 <- 2; z2 <- NULL; k <- 5  
while (TRUE) { 
  z2 <- c(z2, z1 * k)
  k <- k + 1  
  if (k >= 8) break
}
z2; k

# D. repeat loop with break
p1 <- 2; p2 <- NULL; n <- 5  
repeat { 
  p2 <- c(p2, p1 * n)
  n <- n + 1  
  if (n >= 8) break
}
p2; n

# E. repeat loop with an unknown number of loops
input <- 0.1; tolerance <- 0.0001
iteration <- 0; max.iteration <- 100
repeat {
  output <- input * input
  input <- input - output
  iteration <- iteration + 1
  if (abs(output) < tolerance || iteration > max.iteration) break
}
input; output; iteration