## package and data
library("betareg")
data("ReadingSkills", package = "betareg")

## augment with random noise
set.seed(1071)
n <- nrow(ReadingSkills)
ReadingSkills$x1 <- rnorm(n)
ReadingSkills$x2 <- runif(n)
ReadingSkills$x3 <- factor(sample(0:1, n, replace = TRUE))

## fit beta regression tree
rs_tree <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
  data = ReadingSkills, minsplit = 10)

## methods
print(rs_tree)
summary(rs_tree)
coef(rs_tree)
sctest(rs_tree)
