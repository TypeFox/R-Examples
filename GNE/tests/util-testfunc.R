if(!require("GNE"))stop("this test requires package GNE.")

f <- function(x) x
f2 <- function(x, y) x+y
f3 <- function(x, y, z) f2(x,y)+z
f4 <- function(x, y, z, a) f3(x,y,z)+a
f5 <- function(x, y, z, a, b) f4(x,y,z,a)+b

g <- function(x, y, z) GNE:::evalwitharglist(f3, x, list(y, z))
g2 <- function(x, y, z) GNE:::evalwitharglist(f3, x, c(list(y), z))

h <- function(x, z, y) x+z+y$a+y$b


args(GNE:::testfunc)
GNE:::testfunc(f, echo=2, errmess="1")
GNE:::testfunc(f, 1, echo=2, errmess="2")
try( GNE:::testfunc(f, 1, 1:10, echo=2, errmess="3") )
GNE:::testfunc(f2, 1, 1:10, echo=2, errmess="4")
GNE:::testfunc(f3, 1, 1:10, echo=2, errmess="5")
GNE:::testfunc(f3, 1, NULL, echo=2, errmess="6")
GNE:::testfunc(f4, 1, echo=2, errmess="7")
GNE:::testfunc(f5, 1, echo=2, errmess="8")

GNE:::testfunc(g, 1, 1:10, echo=2, errmess="9")
GNE:::testfunc(g2, 1, list(1:10), echo=2, errmess="10")

GNE:::testfunc(h, 1, list(a=1, b=pi), echo=2, errmess="11")
GNE:::testfunc(h, arg=list(a=1, b=pi), echo=2, errmess="12")

