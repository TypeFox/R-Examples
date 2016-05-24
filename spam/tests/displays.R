# This is file ../spam/tests/displays.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








options( echo=FALSE)
library( spam, warn.conflict=FALSE)

# This script illustrates the plotting capacities


# the following function is form fields (and should be made default)  
tim.colors <- function (n = 64)
{
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
            "#AF0000", "#9F0000", "#8F0000", "#800000")
    if (n == 64)
        return(orig)
    rgb.tim <- t(col2rgb(orig))
    temp <- matrix(NA, ncol = 3, nrow = n)
    x <- seq(0, 1, , 64)
    for (k in 1:3) {
      # the original function uses splint here.
        hold <- spline(x, rgb.tim[, k], n)$y
        hold[hold < 0] <- 0
        hold[hold > 255] <- 255
        temp[, k] <- round(hold)
    }
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}


m <- 10
n <- 7

set.seed(124)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0
ss <- as.spam(tt)



par(mfcol=c(1,2))
spam.options(imagesize=10)
display(ss,cex=1)

spam.options(imagesize=n*m+1)
display(ss)



par(mfcol=c(1,2))
plot(tt)
plot(ss)


plot(      tt[,1])
plot(as.spam( tt[,1]))

plot(      t( tt[,2]))
plot(as.spam( t( tt[,2])))



nl <- length(ss)  #ok
ss@entries <- 1:nl
z <- ss
br <- c(seq(0.1,max(z),l=nl),max(z))
par(mfcol=c(1,2))
spam.options( imagesize=1000)
image(z, breaks=br,col=tim.colors(nl))
spam.options( imagesize=10)
image(z, breaks=br,col=tim.colors(nl),cex=1)


nl <- length(ss)
ss@entries <- 1:nl
par(mfcol=c(1,2))
spam.options( imagesize=1000)
image(ss, col=tim.colors(nl))

spam.options( imagesize=10)
image(ss, col=tim.colors(nl),cex=1)


# very large sample
nz <- 128
ln <- nz^2
ss <- spam(0,ln,ln)
for (i in 1:nz) ss[sample(ln,1),sample(ln,1)] <- i
par( mfcol=c(1,1))
image(ss, col=tim.colors(nl),cex=100)





