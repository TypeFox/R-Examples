test.formatDistData <- function() {
    dat <- data.frame(distance=1:100, site=gl(5, 20),
                      visit=factor(rep(1:4, each=5)))
    cutpt <- seq(0, 100, by=25)
    y <- formatDistData(dat, "distance", "site", cutpt)
    checkEqualsNumeric(y, matrix(c(20,   0,   0,   0,
                                    5,  15,   0,   0,
                                    0,  10,  10,   0,
                                    0,   0,  15,   5,
                                    0,   0,   0,  20), 5, 4, byrow=TRUE))
    dat.bad <- dat
    dat.bad$distance <- as.character(dat$distance)
    checkException(formatDistData(dat.bad, "distance", "site", cutpt))

    dat.bad <- dat
    dat.bad$site <- as.character(dat$site)
    y2 <- formatDistData(dat.bad, "distance", "site", cutpt)
    checkEqualsNumeric(y2, matrix(c(20,   0,   0,   0,
                                    5,  15,   0,   0,
                                    0,  10,  10,   0,
                                    0,   0,  15,   5,
                                    0,   0,   0,  20), 5, 4, byrow=TRUE))

    y3 <- formatDistData(dat, "distance", "site", cutpt, "visit")
    checkEqualsNumeric(y3, matrix(c(
5, 0, 0, 0,   5, 0, 0, 0,   5, 0, 0, 0,   5, 0, 0, 0,
5, 0, 0, 0,   0, 5, 0, 0,   0, 5, 0, 0,   0, 5, 0, 0,
0, 5, 0, 0,   0, 5, 0, 0,   0, 0, 5, 0,   0, 0, 5, 0,
0, 0, 5, 0,   0, 0, 5, 0,   0, 0, 5, 0,   0, 0, 0, 5,
0, 0, 0, 5,   0, 0, 0, 5,   0, 0, 0, 5,   0, 0, 0, 5), 5, 16, byrow=TRUE))

}
