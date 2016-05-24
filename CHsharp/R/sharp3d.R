"sharp3d" <-
function (x, y, z, hspace=1, htime=1, v=1) 
{
n <- length(x)
xsharp <- x
ysharp <- y
zsharp <- z
hsharp <- hspace
z <- .Fortran("sharp3d", as.integer(n), as.double(hsharp), as.double(htime),
as.double(x), as.double(y), as.double(z), as.double(xsharp),
as.double(ysharp), as.double(zsharp), as.integer(v), PACKAGE="CHsharp")

names(z) <- c("n", "hsharp", "htime", "x", "y", "z", "xsharp", "ysharp", "zsharp", "v")
data.frame(x.sharp=z$xsharp, y.sharp=z$ysharp, z.sharp=z$zsharp)
}

"sharp3dB" <-
function (x, y, z, hspace=1, htime=1, v=1) 
{
n <- length(x)
xsharp <- x
ysharp <- y
zsharp <- z
hsharp <- hspace
z <- .Fortran("sharp3dB", as.integer(n), as.double(hsharp), 
as.double(htime),
as.double(x), as.double(y), as.double(z), as.double(xsharp),
as.double(ysharp), as.double(zsharp), as.integer(v), PACKAGE="CHsharp")

names(z) <- c("n", "hsharp", "htime", "x", "y", "z", "xsharp", "ysharp", "zsharp", "v")
data.frame(x.sharp=z$xsharp, y.sharp=z$ysharp, z.sharp=z$zsharp)
}

"sharp2d" <-
function (x, y, hspace=1, htime=1, v=1) 
{
n <- length(x)
xsharp <- x
ysharp <- y
hsharp <- hspace
z <- .Fortran("sharp2d", as.integer(n), as.double(hsharp), as.double(htime),
as.double(x), as.double(y), as.double(xsharp),
as.double(ysharp), as.integer(v), PACKAGE="CHsharp")

names(z) <- c("n", "hsharp", "htime", "x", "y", "xsharp", "ysharp", "v")
data.frame(x.sharp=z$xsharp, y.sharp=z$ysharp)
}

"sharp1d" <- function(x, h, v=1) {
   h <- h/sqrt(2)^v
   for (j in 1:v) {
       NW.numerator <- apply(outer(x,x,function(x,y) dnorm(x-y,sd=h)*y), 1, sum)
       NW.denominator <- apply(outer(x,x,function(x,y) dnorm(x-y, sd=h)), 1, sum)
       x <- NW.numerator/NW.denominator
   }
   x.sharp <- x
   x.sharp
}

