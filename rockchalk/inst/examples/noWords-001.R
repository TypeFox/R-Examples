cat("\n\n Basics: R math, workspace, indexes, class, seq, rep, identical, all.equal \n\n")

1 + 4
exp(7)
sqrt(81)
log(4)
4^2

x <- c(1,2,3,4,5,6)
x

ls()
objects() #same

x2 <- seq(1, 6, by = 1)
x2
x3 <- seq(1, 6)
x3
x4 <- 1:6
x4
identical(x, x2, x3, x4)

## caution: why identical if class not same?
class(x)
class(x2)
class(x3)
class(x4)

xint <- as.integer(x)
class(xint)

x2 <- seq(1L, 6L, by = 1L)
class(x2)

x2 <- seq.int(1,6)
class(x2)


ls()
rm(x, x2, x3, x4)
ls()

x <- seq(0.5, 4.1, by = 0.3)
x

x^2
exp(x)
sqrt(x)
log(x)
1/x


x[4]
x[4:6]
x[c(1,4,5)]
x[ x > 4 ]
x[ -2 ]
x[ -c(2,3) ]


attributes(x)
is.null(x)
is.vector(x)
is.numeric(x)
is.integer(x)
is.character(x)
is.data.frame(x)
is.logical(x)
is.matrix(x)
is.list(x)

rm(x)

y <- c(80, 90, 100, 110, 120, 130)
x <- seq_along(y)
x
x + y
x - y
y * x
x * y
y^x
3*x + 4*y

w <- 2*x + sqrt(y)
w
mxyw <- cbind(x,y,w)
class(mxyw)
mxyw[ , 2:3]
mxyw[ , "w"]
mxyw[ , c("y", "w")]
mxyw[ 2:3, "w"]


z <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
is.logical(z)
which(z == TRUE)
which(z)
which(z != TRUE)
which(!z)

print("Exclamation mark (!) means NOT")
cat("Exclamation mark (!) means NOT \n")


zTrueIndex <- which(z == TRUE)
x[zTrueIndex]
y[zTrueIndex]
x[zTrueIndex] + y[zTrueIndex]

xy <- data.frame(x, y, w, z)
xy
colnames(xy)
dim(xy)
fix(xy) #Hit "Quit"

xy[1, ]
xy[ ,2]
xy[1:3, ]
xy[c(1,4,6), ]
xy[ , -3]
xy[-c(1,2), ]
xy[ , "w"]
xy[ , c("x","y")]
xy[ , c("z")]
xy[["z"]]
xy[xy$z, ]
xy[!xy$z, ]


newx1 <- xy[ , "x"]
newx1
is.data.frame(newx1)
class(newx1)
is.vector(newx1)

newx2 <- xy[ , "x", drop=F] # drop: R magic/curse?
newx2
is.data.frame(newx2)
class(newx2)
is.vector(newx2)


all.equal(newx1, newx2)
all.equal(newx1, as.vector(newx2))
identical(newx1, as.vector(newx2))



subset(xy, subset = x < 3)
xy[ x<3, ]
xy[ x<3 & z == TRUE, ]
xy[ x<3 & z, ]
xy[ x<3 & !z, ]





x1 <- c(1,1,1,1,1,1,1,1,1,1)
x2 <- rep(1, times = 10)
x2
x2 <- rep(1, 10)
x1 == x2
identical(x1, x2)



rep(c(1, 5, 8), length.out = 9)
rep(c(1, 5, 8), each = 3)
rep(c(1, 5, 8), each = 3, length.out = 9)
rep(c(1, 5, 8), each = 3, length.out = 10)
rep(c(1, 5, 8), times = 3)
rep(c(1, 5, 8), each = 3, times = 3)





rm(x)
x <- vector(12, mode="integer")
x
x[4] <- 13L
class(x)
is.integer(x)
x[5] <- -5L
is.integer(x)
x[9] <- -11
is.integer(x) # type promotion b/c -11 not integer
class(x)
x <- as.integer(x)
is.integer(x)
class(x)
rm(x)

x <- vector(12, mode="integer")
x
x[4] <- 13
class(x)
is.integer(x)
x[5] <- -5L
is.integer(x)
x[9] <- -11.1
is.integer(x)
x <- as.integer(x) #rounds
is.integer(x)
class(x)


xs1 <- seq(from = 5, to = 20, by=1)
xs2 <- seq(from = 5, to = 20, length.out=16)
xs1 == xs2
identical(xs1, xs2)

xs3 <- 5:20
xs1 == xs3
identical(xs1, xs3) ##hmmm. huh?
is.vector(xs3)
is.vector(xs1)
length(xs1)
length(xs3)
all.equal(xs1, xs3) #ok!

which( xs1 != xs3 )
which( xs1 == xs3 )

for( i in seq_along(xs3)){
  print(xs1[i]-xs3[i])
}
## hmmm. again
is.integer(xs1)
is.numeric(xs1)
is.numeric(xs3)
is.integer(xs3) ## aha!
class(xs1)
class(xs2)
class(xs3)

identical(xs3, as.integer(xs1)) ##puzzle solved


