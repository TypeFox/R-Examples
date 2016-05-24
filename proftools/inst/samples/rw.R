rw <- function(n = 100) {
    xpos <- ypos <- numeric(n)
    truefalse <- c(TRUE, FALSE)
    plusminus1 <- c(1, -1)
    for(i in 2:n) {
        ## Decide whether we are moving horizontally
        ## or vertically.
        if (sample(truefalse, 1)) {
            xpos[i] <- xpos[i-1] + sample(plusminus1, 1)
            ypos[i] <- ypos[i-1]
        }
        else {
            xpos[i] <- xpos[i-1]
            ypos[i] <- ypos[i-1] + sample(plusminus1, 1)
        }
    }
    list(x = xpos, y = ypos)
}

