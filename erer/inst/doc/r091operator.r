# A binary operator can be reformatted as a function
aa <- 3 + c(5, 6)
bb <- "+"(3, c(5, 6))

# Four types of operators: arithmetic, assignment, relational, logical
dog <- c(1, 3, 5, 2 ^ 4, 70, 100 %/% 8)
pig <- c(1, 2, 6) + 1
cow <- 70
r1 <- dog == pig; r2 <- dog < cow
r3 <- r1 &  r2;   r4 <- r1 |  r2
r5 <- r1 && r2;   r6 <- r1 || r2
r7 <- all(r3);    r8 <- any(r4)
pig + c(10, 20)

# Defining a new binary operator
"%score%" <- function(x, y) {40 + x + y}
exam.a <- c(23, 26, 29); exam.b <- c(18, 14, 10)
exam.total <- exam.a %score% exam.b
exam.total

# Understanding precedence rules of operators 
 3 ^ 2 + 4 > 60  &  10 != 2.5   # FALSE
(3 ^ 2 + 4 > 60) & (10 != 2.5)  # FALSE; clearer with the parentheses

# Rounding errors
71 / 3 * 3 == 71              # TRUE
all.equal(71 / 3 * 3, 71)     # TRUE
all.equal(10 / 3, 3.3333333)  # TRUE
identical(10 / 3, 3.3333333)  # FALSE
all.equal(pi, 355 / 113)