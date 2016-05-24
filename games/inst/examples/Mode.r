x <- c(1, 3, 3, 4)
Mode(x)  # 3
x.char <- letters[x]
Mode(x.char)  # "c"
x.factor <- as.factor(x.char)
Mode(x.factor)  # "c", with levels "a", "c", "d"
x.logical <- x > 3
Mode(x.logical)  # FALSE

## Behavior with ties
y <- c(3, 3, 1, 1, 2)
Mode(y)  # 3
z <- c(2, 1, 1, 3, 3)
Mode(z)  # 1
