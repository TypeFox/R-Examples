col.alpha <- function (col, alpha=0.5) 
{
        COL <- col2rgb(col)/255
        rgb(red=COL[1], green=COL[2], blue=COL[3], alpha=alpha)
}

