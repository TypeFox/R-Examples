"draw.bplot" <-

function (temp, width, xpos, outlier = TRUE, style = "tukey") 

{

    if (temp$N < 1) 

        return()

    if (style == "quantile") {

        temp <- temp[!is.na(temp)]

        quant <- c(0.05, 0.25, 0.5, 0.75, 0.95)

        bb <- quantile(temp, quant)

        mid <- xpos

        low <- mid - width * 0.5

        high <- mid + width * 0.5

        if (length(temp) > 5) {

            y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], 

                bb[2], bb[4])

            x <- c(high, low, NA, mid, mid, NA, high, low, low)

            y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], 

                bb[5], bb[5], bb[5])

            x <- c(x, high, high, high, low, NA, mid, mid, high, 

                low)

            lines(x, y)

        }

        if (length(temp) > 5) {

            outs <- temp[(temp < bb[1]) | (temp > bb[5])]

        }

        else outs <- temp

        olen <- length(outs)

        if ((olen > 0) & outlier) 

            points(rep(mid, olen), outs)

    }

    if (style == "tukey") {

        temp <- temp[!is.na(temp)]

        quant <- c(0.05, 0.25, 0.5, 0.75, 0.95)

        bb <- quantile(temp, quant)

        iqr <- bb[4] - bb[2]

        mid <- xpos

        low <- mid - width * 0.5

        high <- mid + width * 0.5

        bb[1] <- min(temp[temp >= bb[2] - 1.5 * iqr])

        bb[5] <- max(temp[temp <= bb[4] + 1.5 * iqr])

        if (length(temp) > 5) {

            y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], 

                bb[2], bb[4])

            x <- c(high, low, NA, mid, mid, NA, high, low, low)

            y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], 

                bb[5], bb[5], bb[5])

            x <- c(x, high, high, high, low, NA, mid, mid, high, 

                low)

            lines(x, y)

        }

        if (length(temp) > 5) {

            outs <- temp[(temp < bb[2] - 3 * iqr) | (temp > bb[4] + 

                3 * iqr)]

        }

        else outs <- temp

        olen <- length(outs)

        if ((olen > 0) & outlier) 

            points(rep(mid, olen), outs)

    }

}

