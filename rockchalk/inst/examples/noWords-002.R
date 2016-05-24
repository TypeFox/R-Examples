cat("\n\n Basics: vector, matrix, attribute, attr, dim, structure \n\n, t, %*%, outer, crossprod")


x <- c(13, 11, 18, 2, 134, 154, 2, 3, 4, 12)
is.vector(x)
is.matrix(x)
attributes(x)
?attributes

matrix(x, nrow=2)
?matrix

X1m1 <- matrix(x, nrow=2)
class(X1m1)
X1m2 <- matrix(x, ncol=5)
X1m1 == X1m2

identical(X1m1, X1m2)
?identical

attributes(X1m1)
attributes(X1m2)

X1m3 <- x
attr(X1m3, "dim") <- c(2,5)
?attr
X1m3
is.vector(X1m3)
is.matrix(X1m3)
is.data.frame(X1m3)
attributes(X1m3)
identical(X1m1, X1m3)


X1m4 <- x
dim(X1m4) <- c(2,5)
?dim
X1m4
attributes(X1m4)

identical(X1m1, X1m4)


X1m5 <- structure(x, dim=c(2,5))
?structure
X1m5

identical(X1m1, X1m5)


X1m6 <- rbind(x[c(1,3,5,7,9)], x[c(2,4,6,8,10)])
class(x[1:5])
class(X1m6)
identical(X1m1, X1m6)

X1m7 <- cbind(x[1:2], x[3:4], x[5:6], x[7:8],x[9:10])
identical(X1m1, X1m7)


x1m7 <- X1m7
identical(x, x1m7)
dim(x1m7) <- c(1,10)
x1m7
identical(x, x1m7) #humphf!

x1m7new1 <- drop(x1m7)
?drop
x1m7new1
identical(x, x1m7new1) #whew


x1m7new2 <-  as.vector(x1m7)
identical(x, x1m7new2)



?t
tX1m1 <- t(X1m1)
tX1m1

X2m1 <- matrix(x, nrow=5, byrow=TRUE)
X2m1
dim(X2m1)
identical(X2m1, t(X1m1))



X2m2 <- X1m1
X2m2
dim(X2m2)
dim(X2m2) <- c(5,2)
dim(X2m2)
X2m2
identical(X1m1, X2m2)
X2m1 == X2m2 #ack!



y <- c(12, 12, 41, 22)
x %*% y #heck

y <- rep(c(1,2,3,4,5), 2)
x %*% y
t(x) %*% y
y %*% x


outer(x,y)
outer(x,y, FUN="*")
outer(x,y, FUN="+")
outer(x,y, FUN="/")
outer(x,y, FUN=function(a, b) {b - 2*a})


X3 <- matrix(rnorm(10), nrow=2)

crossprod(X1m1, X3)
t(X1m1) %*% X3
?crossprod



## Puzzle

x <- 1:6

X <- data.frame(x1 = x, x2 = x, x3 = x)
Y <- data.frame(z1 = rnorm(6), z2 = rnorm(6), x3 = rnorm(6))

Z <- cbind(X, Y)
colnames(Z)

Z2 <- data.frame(X,Y)
colnames(Z2)


