
"sm.sphere" <- function (lat, long, kappa = 20, hidden = FALSE, sphim = FALSE,
          addpoints = FALSE, ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 32)
    ngrid <- opt$ngrid
    panel <- opt$panel
    phi <- opt$phi
    theta <- opt$theta
    kap <- kappa
    invis <- plot2d(lat, long, theta, phi)
    sphdraw(theta, phi)
    if (!opt$panel) {
        if (hidden)
            hidplot(invis, theta, phi)
        if (sphim)
            sphimage(lat, long, kap, theta, phi, ngrid)
        if (sphim & addpoints)
            addplot(lat, long, theta, phi)
    }
    else {
        items <- c("Set theta and phi", "  - increase theta",
            "  - decrease theta", "  - increase phi", "  - decrease phi",
            "Add hidden points", "Add density estimate", "  - increase s.p.",
            "  - decrease s.p.", "  - add data points", "Exit")
        ind <- menu(items, graphics = TRUE, title = "Sphere")
        while (items[ind] != "Exit") {
            if (items[ind] == "Set theta and phi") {
                a <- change(theta, phi)
                theta <- a$theta
                phi <- a$phi
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - increase theta") {
                theta <- inctheta(theta, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - decrease theta") {
                theta <- inctheta(theta, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - increase phi") {
                phi <- incphi(phi, 30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "  - decrease phi") {
                phi <- incphi(phi, -30)
                invis <- plot2d(lat, long, theta, phi)
                sphdraw(theta, phi)
            }
            else if (items[ind] == "Add hidden points") {
                hidplot(invis, theta, phi)
            }
            else if (items[ind] == "Add density estimate") {
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - increase s.p.") {
                kap <- kap * 2
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - decrease s.p.") {
                kap <- kap/2
                sphimage(lat, long, kap, theta, phi, ngrid)
            }
            else if (items[ind] == "  - add data points") {
                par(pch = "*")
                addplot(lat, long, theta, phi)
            }
            else if (items[ind] == "Add 2nd data set") {
                par(pch = "x")
                addplot(lat2, long2, theta, phi)
            }
            ind <- menu(items, graphics = TRUE, title = "Sphere")
        }
    }
    par(pty = "m")
    invisible(list(theta = theta, phi = phi, kappa = kap))
}

"sphdraw" <- function (theta, phi) {
    a1 <- 0
    a2 <- 30
    a3 <- 60
    a4 <- 90
    a5 <- 120
    a6 <- 150
    b1 <- (-90)
    b2 <- (-60)
    b3 <- (-30)
    b4 <- 0
    b5 <- 30
    b6 <- 60
    b7 <- 90
    latlines(b1, theta, phi)
    latlines(b2, theta, phi)
    latlines(b3, theta, phi)
    latlines.e(b4, theta, phi)
    latlines(b5, theta, phi)
    latlines(b6, theta, phi)
    latlines(b7, theta, phi)
    longlines.e(a1, theta, phi)
    longlines(a2, theta, phi)
    longlines(a3, theta, phi)
    longlines(a4, theta, phi)
    longlines(a5, theta, phi)
    longlines(a6, theta, phi)
    circle(1)
    }

"sphimage" <- function (latitude, longitude, kap, theta, phi, ngrid = 32) {
    values <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
    xgrid <- rep(values, rep(ngrid, ngrid))
    ygrid <- rep(values, ngrid)
    dvec <- rep(0, ngrid^2)
    xlong <- longitude * pi/180
    xlat <- latitude * pi/180
    n <- length(longitude)
    radtheta <- theta * pi/180
    radphi <- phi * pi/180
    xgrid[xgrid^2 + ygrid^2 >= 1] <- NA
    ygrid[xgrid^2 + ygrid^2 >= 1] <- NA
    za <- -xgrid * sin(radtheta) - ygrid * cos(radtheta) * sin(radphi)
    zb <- cos(radphi) * cos(radtheta) * sqrt(1 - xgrid^2 - ygrid^2)
    z <- za + zb
    if ((theta == 90) | (theta == 270))
        x <- -ygrid * sin(radtheta) * sin(radphi) + cos(radphi) *
            sqrt(1 - ygrid^2 - z^2)
    else x <- (xgrid + z * sin(radtheta))/cos(radtheta)
    if (phi == 90)
        y <- sqrt(1 - x^2 - z^2)
    else if (phi == -90)
        y <- -sqrt(1 - x^2 - z^2)
    else y <- (ygrid + (x * sin(radtheta) + z * cos(radtheta)) *
        sin(radphi))/cos(radphi)
    xyzok <- (((x/sqrt(x^2 + z^2)) * (sqrt(1 - y^2)) * sin(radtheta) *
        cos(radphi)) + (y * sin(radphi)) - ((-z/sqrt(x^2 + z^2)) *
        (sqrt(1 - y^2)) * cos(radphi) * cos(radtheta)))
    other <- !is.na(xyzok) & xyzok < 0
    z[other] <- (za - zb)[other]
    x[other] <- ((xgrid + (z * sin(radtheta)))/cos(radtheta))[other]
    y[other] <- ((ygrid + ((x * sin(radtheta)) + (z * cos(radtheta))) *
                   sin(radphi))/cos(radphi))[other]
    xj <- cos(xlong) * cos(xlat)
    yj <- sin(xlat)
    zj <- -sin(xlong) * cos(xlat)
    dvec <- exp(kap * cbind(x, y, z) %*% rbind(xj, yj, zj)) %*% rep(1/n, n)
    dvec[is.na(xgrid)] <- 0
    dvec <- dvec/max(dvec)
    fmat <- matrix(dvec, ngrid, ngrid, byrow = TRUE)
    x <- seq(-1 + 1/ngrid, 1 - 1/ngrid, length = ngrid)
    y <- x
    image(x, y, fmat, add = TRUE)
    angle <- seq(0, pi/2, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 0, 1, 1), c(yy, 1, 1, 0), col = "white", border = "white")
    angle <- seq(pi/2, pi, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, -1, -1, 0), c(yy, 0, 1, 1), col = "white", border = "white")
    angle <- seq(pi, 3 * pi/2, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 0, -1, -1), c(yy, -1, -1, 0), col = "white", border = "white")
    angle <- seq(3 * pi/2, 2 * pi, length = 50)
    xx <- cos(angle)
    yy <- sin(angle)
    polygon(c(xx, 1, 1, 0), c(yy, 0, -1, -1), col = "white", border = "white")
    sphdraw(theta, phi)
}

"addplot" <- function (d, f, theta, phi) {
    a <- (f * pi)/180
    b <- (d * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    llong <- a[xyzcheck >= 0]
    llat <- b[xyzcheck >= 0]
    if (length(llat) == 0) {
        break
    }
    X <- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) *
        cos(llat) * sin(radtheta))
    Y <- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) *
        ((sin(llong) * cos(radtheta)) - (cos(llong) * sin(radtheta))))
    par(pty = "s")
    points(X, Y)
}

"change" <-
function (th, ph)
{
    cat("Theta =", th, "\n")
    cat("Phi =", ph, "\n")
    scan(n = 1)
    cat("Change theta to ? \n")
    theta <- scan(n = 1)
    if (theta >= 360) theta <- theta - 360
    cat("\n", "Change phi to ? \n")
    phi <- scan(n = 1)
    if (phi > 90) phi <- 90
    if (phi < -90) phi <- -90
    cat("Theta =", theta, "\n")
    cat("Phi =", phi, "\n")
    list(theta = theta, phi = phi)
}

"circle" <- 
function (r)
{
    angle <- seq(0, 7, by = 0.1)
    x <- r * cos(angle)
    y <- r * sin(angle)
    par(lty = 1)
    lines(x, y)
}

"hidplot" <-
function (invis, theta, phi)
{
    invislong <- invis$invislong
    invislat <- invis$invislat
    par(pch = "O")
    a <- (invislong * pi)/180
    b <- (invislat * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    if (length(invislat) == 0) {
        points(0, 0, type = "n")
        break
    }
    X <- (cos(invislong) * cos(invislat) * cos(radtheta)) +
        (sin(invislong) * cos(invislat) * sin(radtheta))
    Y <- (sin(invislat) * cos(radphi)) + ((cos(invislat) * sin(radphi)) *
        ((sin(invislong) * cos(radtheta)) - (cos(invislong) *
            sin(radtheta))))
    points(X, Y)
}

"incphi" <-
function (ph, inc)
{
    phi <- ph + inc
    if (phi > 90) phi <- 90
    if (phi < -90) phi <- -90
    cat("Phi =", phi, "\n")
    phi
}

"inctheta" <-
function (th, inc)
{
    theta <- th + inc
    if (theta >= 360) theta <- theta - 360
    theta
}

"latlines" <-
function (beta, theta, phi)
{
    if (beta < (phi - 90) | beta > (phi + 90)) return()
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radbeta <- (beta * pi)/180
    radphi <- (phi * pi)/180
    alpha <- seq(0, (2 * pi), by = 0.05)
    xyzcheck <- ((cos(alpha) * cos(radbeta) * sin(radtheta) *
        cos(radphi)) + (sin(radbeta) * sin(radphi)) - (sin(alpha) *
        cos(radbeta) * cos(radphi) * cos(radtheta)))
    alphaplot <- alpha[xyzcheck >= 0]
    X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) + (sin(alphaplot) *
        cos(radbeta) * sin(radtheta))
    Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) * cos(radtheta)) -
        (cos(alphaplot) * sin(radtheta))) * cos(radbeta) * sin(radphi))
    points(X, Y)
}

"latlines.e" <-
function (beta, theta, phi)
{
    if (beta < (phi - 90) | beta > (phi + 90)) return()
    par(lty = 2)
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radbeta <- (beta * pi)/180
    radphi <- (phi * pi)/180
    alpha <- seq(0, (2 * pi), by = 0.005)
    xyzcheck <- ((cos(alpha) * cos(radbeta) * sin(radtheta) *
        cos(radphi)) + (sin(radbeta) * sin(radphi)) - (sin(alpha) *
        cos(radbeta) * cos(radphi) * cos(radtheta)))
    alphaplot <- alpha[xyzcheck >= 0]
    X <- (cos(alphaplot) * cos(radbeta) * cos(radtheta)) + (sin(alphaplot) *
        cos(radbeta) * sin(radtheta))
    Y <- (sin(radbeta) * cos(radphi)) + (((sin(alphaplot) * cos(radtheta)) -
        (cos(alphaplot) * sin(radtheta))) * cos(radbeta) * sin(radphi))
    points(X, Y)
}

"longlines" <-
function (alpha, theta, phi)
{
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radalpha <- (alpha * pi)/180
    radphi <- (phi * pi)/180
    beta <- seq(0, (2 * pi), by = 0.05)
    xyzcheck <- ((cos(radalpha) * cos(beta) * sin(radtheta) *
        cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) *
        cos(beta) * cos(radphi) * cos(radtheta)))
    betaplot <- beta[xyzcheck >= 0]
    X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) + (sin(radalpha) *
        cos(betaplot) * sin(radtheta))
    Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) * cos(radtheta)) -
        (cos(radalpha) * sin(radtheta))) * cos(betaplot) * sin(radphi))
    points(X, Y)
}

"longlines.e" <-
function (alpha, theta, phi)
{
    par(lty = 2)
    par(pch = ".")
    radtheta <- (theta * pi)/180
    radalpha <- (alpha * pi)/180
    radphi <- (phi * pi)/180
    beta <- seq(0, (2 * pi), by = 0.005)
    xyzcheck <- ((cos(radalpha) * cos(beta) * sin(radtheta) *
        cos(radphi)) + (sin(beta) * sin(radphi)) - (sin(radalpha) *
        cos(beta) * cos(radphi) * cos(radtheta)))
    betaplot <- beta[xyzcheck >= 0]
    X <- (cos(radalpha) * cos(betaplot) * cos(radtheta)) + (sin(radalpha) *
        cos(betaplot) * sin(radtheta))
    Y <- (sin(betaplot) * cos(radphi)) + (((sin(radalpha) * cos(radtheta)) -
        (cos(radalpha) * sin(radtheta))) * cos(betaplot) * sin(radphi))
    points(X, Y)
}

"plot2" <-
function (latitude2, longitude2, theta, phi)
{
    par(pch = "x")
    a <- (longitude2 * pi)/180
    b <- (latitude2 * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    long2 <- a[xyzcheck >= 0]
    lat2 <- b[xyzcheck >= 0]
    if (length(lat2) == 0) {
        points(0, 0, type = "n")
        text(0.6, -1.2, labels = "Data set:")
        break
    }
    X <- (cos(long2) * cos(lat2) * cos(radtheta)) + (sin(long2) *
        cos(lat2) * sin(radtheta))
    Y <- (sin(lat2) * cos(radphi)) + ((cos(lat2) * sin(radphi)) *
        ((sin(long2) * cos(radtheta)) - (cos(long2) * sin(radtheta))))
    points(X, Y)
}

"plot2d" <-
function (d, f, theta, phi)
{
    par(pch = "*")
    a <- (f * pi)/180
    b <- (d * pi)/180
    radtheta <- (theta * pi)/180
    radphi <- (phi * pi)/180
    xyzcheck <- ((cos(a) * cos(b) * sin(radtheta) * cos(radphi)) +
        (sin(b) * sin(radphi)) - (sin(a) * cos(b) * cos(radphi) *
        cos(radtheta)))
    llong <- a[xyzcheck >= 0]
    llat <- b[xyzcheck >= 0]
    invislong <- a[xyzcheck < 0]
    invislat <- b[xyzcheck < 0]
    if (length(llat) == 0) {
        par(pty = "s")
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
            xlim = c(-1, 1), ylim = c(-1, 1))
        list(invislong = invislong, invislat = invislat)
        break
    }
    X <- (cos(llong) * cos(llat) * cos(radtheta)) + (sin(llong) *
        cos(llat) * sin(radtheta))
    Y <- (sin(llat) * cos(radphi)) + ((cos(llat) * sin(radphi)) *
        ((sin(llong) * cos(radtheta)) - (cos(llong) * sin(radtheta))))
    par(pty = "s")
    plot(X, Y, axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1),
         ylim = c(-1, 1))
    list(invislong = invislong, invislat = invislat)
}
