
library("basefun")

n <- 100
x <- runif(n)
y <- runif(n)
g <- gl(4, n / 4)
d <- data.frame(y = y, x = x, g = g)
de <- d[-(1:nrow(d)),]

b1 <- as.basis(~ x + g, data = d)
b2 <- Bernstein_basis(numeric_var("y", support = c(0.0, 1.0)), 
                      order = 4, ui = "incre")

b12 <- b(b1 = b1, b2 = b2)
c12 <- c(b1 = b1, b2 = b2)

X1 <- model.matrix(b(b2 = b2, b1 = b1), data = d)
X2 <- model.matrix(b(b2 = b2, b1 = b1, sumconstr = TRUE), data = d)
all.equal(X1[,], X2[,])
attr(X1, "constraint")$ui
attr(X2, "constraint")$ui

dim(model.matrix(b12, d))
# nparm(b12, d)
dim(model.matrix(c12, d))
# nparm(c12, d)

tmp <- c(b12 = b12, c12 = c12)
class(tmp)
dim(model.matrix(tmp, d))
# nparm(tmp, d)

xd <- data.frame(x = x)
b <- as.basis(~ scale(x), data = xd)
stopifnot(all.equal(b(xd)[1:10,], b(xd[1:10,,drop = FALSE])[,]))

