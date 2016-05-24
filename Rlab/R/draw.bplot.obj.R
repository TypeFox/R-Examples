"draw.bplot.obj" <-

function (obj, width, xpos, outlier = TRUE, horizontal = FALSE) 

{

    N <- obj$N

    bb <- obj$bb

    mid <- xpos

    low <- mid - width * 0.5

    high <- mid + width * 0.5

    if (N > 5) {

        y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], bb[2], 

            bb[4])

        x <- c(high, low, NA, mid, mid, NA, high, low, low)

        y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], bb[5], 

            bb[5], bb[5])

        x <- c(x, high, high, high, low, NA, mid, mid, high, 

            low)

        if (horizontal) {

            lines(y, x)

        }

        else {

            lines(x, y)

        }

    }

    outs <- obj$out

    olen <- length(outs)

    if ((olen > 0) & outlier) {

        if (horizontal) {

            points(outs, rep(mid, olen))

        }

        else {

            points(rep(mid, olen), outs)

        }

    }

}

