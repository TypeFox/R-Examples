# Index on vector
v <- c(10, 20, 30, 40, 50); names(v) <- LETTERS[1:length(v)]
v1 <- v[c(2:3, 5)]; v1 <- v[ -c(1, 4)]; v1 <- v[c(0, 2, 3, 0, 5)]
v2 <- v[c(1:5, 1:3)]
v1 <- v[c("B", "C", "E")] 
v1 <- v[-match(x = c("A", "D"), table = names(v))] # character to num index
v1 <- v[c(FALSE, TRUE, TRUE, FALSE, TRUE)]
v1 <- v[c(FALSE, TRUE, TRUE)]
v3 <- v[v > 30]; v3 <- v[which(v > 30)]  # logical to numeric index
v; v1; v2; v3

# Index on data frame
d <- data.frame(y = 11:13, lett = I(letters[1:3])) 
d1 <- d[1:2, 2, drop = FALSE] 
d1 <- d[1:2, "lett", drop = FALSE]
d1 <- d[d$y < 13, "lett", drop = FALSE]                # mix logic/char
d1 <- d[which(d$y < 13), "lett", drop = FALSE]
d2 <- d[c(1:3, 1), c(2, 1, 2, 1, 2), drop = FALSE]     # larger new object
d3 <- d[rev(1:3), c(2,1)]                              # reorder
d1 <- subset(x = d, subset = y < 13, select = "lett")
d1 <- subset(x = d, subset = y < 13, select = lett)
d1 <- subset(x = d, subset = y < 13, select = 2)
d1 <- subset(x = d, subset = c(TRUE, TRUE, FALSE, FALSE, FALSE), 
  select = lett)
d; d1; d2

# Index on list
b <- list(dog = v, chick = d, pig = letters[1:4])
b1 <- b[[1]]; b1 <- b[["dog"]]; b1 <- b$dog; b1  # a new vector
b2 <- b[1:2]; b2 <- b[c("dog", "chick")]; b2     # a new list
b3 <- b[["chick"]][3, ]; b3                      # two layers of indexes

# Index on array
(a <- array(data = 1:30, dim = c(3, 5, 2)))
(a1 <- a[2, 3, 1]); (a2 <- a[, 3, 2])