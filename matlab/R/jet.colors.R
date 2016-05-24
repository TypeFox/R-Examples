###
### $Id: jet.colors.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Create a vector of colors beginning with dark blue, ranging through
### shades of blue, cyan, green, yellow and red, and ending with dark red.
###

DEFAULT.COLORMAP.SIZE <- 64


##-----------------------------------------------------------------------------
jet <- function(m = DEFAULT.COLORMAP.SIZE) {
    n <- matlab::ceil(m / 4)
    u <- c(seq(1, n) / n,
           if (n > 1) {
               matlab::ones(1, n-1)
           } else {
               NULL
           },
           seq(n, 1, by = -1) / n)

    g <- matlab::ceil(n / 2) - (matlab::mod(m, 4) == 1) + (1:length(u))
    r <- g + n
    b <- g - n

    g <- g[!(g > m)]
    r <- r[!(r > m)]
    b <- b[!(b < 1)]

    J <- matlab::zeros(m, 3)
    J[r, 1] <- u[seq(along = r)]
    J[g, 2] <- u[seq(along = g)]
    J[b, 3] <- u[seq(length(u)-length(b)+1, length(u))]

    return(J)
}


##-----------------------------------------------------------------------------
jet.colors <- function(n) {
    if ((n <- as.integer(n[1])) > 0) {
        ans <- jet(n)

        return(rgb(red = ans[,1], green = ans[,2], blue = ans[,3]))
    }
    else {
        return(character(0))
    }
}

