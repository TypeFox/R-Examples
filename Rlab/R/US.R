"US" <-

function (xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = FALSE, 

    shift = FALSE, asp = 1, ...) 

{

    if (!exists("US.dat")) 

        data(US.dat)

    if (shift) {

        US.dat$x[US.dat$x < 0] <- US.dat$x[US.dat$x < 0] + 360

    }

    if (!add) {

        plot(US.dat$x, US.dat$y, ylim = ylim, xlim = xlim, xlab = "", 

            ylab = "", type = "n", xaxt = "n", yaxt = "n", asp = asp, ...)

    }

    lines(US.dat$x, US.dat$y, err = -1, ...)

    invisible()

}

