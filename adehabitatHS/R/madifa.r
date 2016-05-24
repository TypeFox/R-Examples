#######################################################################
#######################################################################
#######                                                          ######
#######                                                          ######
#######                       MADIFA                             ######
#######                                                          ######
#######                                                          ######
#######################################################################
#######################################################################


madifa <- function(dudi, pr, scannf = TRUE, nf = 2)
{
    ## Verifications
    if (!inherits(dudi, "dudi"))
        stop("object of class dudi expected")
    call <- match.call()
    if (any(is.na(dudi$tab)))
        stop("na entries in table")
    if (!is.vector(pr))
        stop("pr should be a vector")

    ## Bases for the analysis
    prb <- pr
    pr <- pr/sum(pr)
    row.w <- dudi$lw
    col.w <- dudi$cw
    Z <- as.matrix(dudi$tab)
    n <- nrow(Z)
    f1 <- function(v) sum(v * pr)
    center <- apply(Z, 2, f1)
    Z <- sweep(Z, 2, center)
    f2 <- function(v) sum((v^2) * pr)

    ## take into account different weights for the columns
    Ze <- sweep(Z, 2, sqrt(col.w), "*")

    ## Inertia matrices S and G
    DpZ <- apply(Ze, 2, function(x) x*pr)
    Se <- crossprod(Ze, DpZ)
    Ge <- crossprod(Ze, apply(Ze,2,function(x) x*row.w))


    ## S^(-1/2)
    eS <- eigen(Se)
    S12 <- eS$vectors %*% diag(eS$values^(-0.5)) %*% t(eS$vectors)

    ## Eigen structure
    W <- S12 %*% Ge %*% S12
    s <- eigen(W)$values

    ## number of eigenvalues
    if (scannf) {
        barplot(s)
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0 | nf > ncol(Ze))
      nf <- 1

    ## Coordinates of the columns
    tt <- as.data.frame((S12 %*% eigen(W)$vectors))
    ww <- apply(tt, 2, function(x) x/sqrt(col.w))
    norw <- sqrt(diag(t(as.matrix(tt))%*%as.matrix(tt)))
    co <- sweep(ww, 2, norw, "/")

    ## scores of the rows in the distorted space
    li <- Z %*% apply(co, 2, function(x) x*col.w)


    co <- as.data.frame(co)
    li <- as.data.frame(li)

    ## Coordinates of the rows in the distorted space
    varus <- apply(li,2,function(x) sum((x^2)*pr))
    l1 <- sweep(li, 2, sqrt(unlist(varus)), "/")

    ## Mahalanobis distances (total)
    mahasu <- apply(l1, 1, function(x) sqrt(sum(x^2)))

    ## For output
    co <- data.frame(co[,1:nf])
    li <- data.frame(li[,1:nf])
    l1 <- data.frame(l1[,1:nf])
    names(co) <- paste("Axis", (1:nf), sep = "")
    row.names(co) <- dimnames(dudi$tab)[[2]]
    names(li) <- paste("Comp.", (1:nf), sep = "")
    names(l1) <- paste("Comp.", (1:nf), sep = "")
    row.names(li) <- dimnames(dudi$tab)[[1]]
    row.names(l1) <- dimnames(dudi$tab)[[1]]

    ## Correlation with these axes
    corav <- cor(dudi$tab, li)

    ## Output
    madifa <- list(call = call, tab = data.frame(Z), pr = prb, cw = col.w,
                   nf = nf, eig = s, lw = row.w, li = li, l1 = l1,
                   co = co, mahasu = mahasu, cor = corav)
    class(madifa) <- "madifa"
    return(madifa)
  }


print.madifa <- function (x, ...)
{
    if (!inherits(x, "madifa"))
        stop("Object of class 'madifa' expected")
    cat("MADIFA")
    cat("\n$call: ")
    print(x$call)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...")
    cat("\n$nf:", x$nf, "axes saved")
    cat("\n")
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("vector", "length",
                                            "mode", "content")))
    sumry[1, ] <- c("$pr", length(x$pr), mode(x$pr), "vector of presence")
    sumry[2, ] <- c("$mahasu", length(x$mahasu), mode(x$mahasu), "squared Mahalanobis distances")
    sumry[3, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[4, ] <- c("$cw", length(x$lw), mode(x$lw), "column weights")
    sumry[5, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow",
                                            "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$li), ncol(x$li), "row normed scores (variance weighted by $pr = 1)")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$cor", nrow(x$cor), ncol(x$cor), "cor(habitat var., scores) for available points")
    class(sumry) <- "table"
    print(sumry)
    if (length(names(x)) > 15) {    cat("\nother elements: ")
                                    cat(names(x)[16:(length(x))], "\n")
                                }
}



scatter.madifa <- function (x, xax = 1, yax = 2, pts = FALSE, percent = 95,
                            clabel = 1, side = c("top", "bottom", "none"),
                            Adensity, Udensity, Aangle, Uangle, Aborder,
                            Uborder, Acol, Ucol, Alty,
                            Ulty, Abg, Ubg, Ainch, Uinch, ...)
{
    ## Verifications
    side <- match.arg(side)
    if (!inherits(x, "madifa"))
        stop("Object of class 'madifa' expected")

    ## Graphical settings
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow=c(1,2))

    ## the bases for the graphs
    x1 <- x$l1[, xax]
    x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
    xlim <- range(x1)
    y1 <- x$l1[, yax]
    y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
    ylim <- range(y1)

    ## background graph
    scatterutil.base(dfxy = x$l1[, c(xax, yax)], xax = 1, yax = 2,
                     xlim = xlim, ylim = ylim, grid = TRUE,
                     addaxes = FALSE,
                     cgrid = 1, include.origin = TRUE, origin = c(0, 0),
                     sub = "",
                     csub = 1.25, possub = "bottomleft",
                     pixmap = NULL, contour = NULL,
                     area = NULL, add.plot = FALSE)

    ## adds the points
    if (pts) {

        ## graphical settings
        if (missing(Acol))
            Acol <- gray(0.8)
        if (missing(Ucol))
            Ucol <- "black"
        if (missing(Abg))
            Abg <- gray(0.8)
        if (missing(Ubg))
            Ubg <- "black"
        if (missing(Ainch))
            Ainch <- 0.03
        if (missing(Uinch))
            Uinch <- Ainch * max(x$pr)

        ## the points
        symbols(x$l1[, c(xax, yax)], circles = rep(1, length(x$pr)),
                fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
        symbols(x$l1[x$pr > 0, c(xax, yax)], circles = x$pr[x$pr > 0],
                fg = Ucol, bg = Ubg,
                inches = Uinch, add = TRUE)
        abline(v = 0)
        abline(h = 0)
    } else {

        ## graphical settings
        if (missing(Adensity))
            Adensity <- NULL
        if (missing(Udensity))
            Udensity <- NULL
        if (missing(Aangle))
            Aangle <- 45
        if (missing(Uangle))
            Uangle <- 45
        if (missing(Aborder))
            Aborder <- NULL
        if (missing(Uborder))
            Uborder <- NULL
        if (missing(Acol))
            Acol <- gray(0.95)
        if (missing(Ucol))
            Ucol <- gray(0.6)
        if (missing(Alty))
            Alty <- NULL
        if (missing(Ulty))
            Ulty <- NULL

        ## adds mcps
        pcff <- function(xy)
        {
            mo <- apply(xy,2,mean)
            dis <- apply(xy, 1, function(x) sum((x-mo)^2))
            xy <- xy[dis < quantile(dis, percent/100),]
            return(xy[chull(xy[,1], xy[,2]),])
        }
        mcpA <- pcff(x$l1[, c(xax, yax)])
        mcpU <- pcff(x$l1[rep(1:length(x$pr), x$pr), c(xax, yax)])
        polygon(mcpA, density = Adensity, angle = Aangle,
                border = Aborder, col = Acol, lty = Alty)
        polygon(mcpU, density = Udensity, angle = Uangle,
                border = Uborder, col = Ucol, lty = Ulty)
        abline(v = 0)
        abline(h = 0)
    }

    ## Bases for the scores for the columns
    dfarr <- x$co[, c(xax, yax)]
    born <- par("usr")
    k1 <- min(dfarr[, 1])/born[1]
    k2 <- max(dfarr[, 1])/born[2]
    k3 <- min(dfarr[, 2])/born[3]
    k4 <- max(dfarr[, 2])/born[4]
    k <- c(k1, k2, k3, k4)
    dfarr <- 0.75 * dfarr/max(k)

    ## Legend
    xax <- paste("Axis", xax)
    yax <- paste("Axis", yax)
    if (side != "none") {
        tra <- paste(" xax =", xax, "\n yax =", yax)
        wt <- strwidth(tra, cex = 1)
        ht <- strheight(tra, cex = 1) * 1.5
        xl <- par("usr")[1]
        yu <- par("usr")[4]
        yd <- par("usr")[3]
        if (side == "top") {
            rect(xl, yu - ht, xl + wt, yu, col = "white", border = 0)
            text(xl + wt/2, yu - ht/2, tra, cex = 1)
        }
        if (side == "bottom") {
            rect(xl, yd + ht, xl + wt, yd, col = "white", border = 0)
            text(xl + wt/2, yd + ht/2, tra, cex = 1)
        }
    }
    box()

    ## column scores
    s.arrow(x$co, clabel = clabel)
    box()
}




hist.madifa <- function (x, scores = TRUE, type = c("h", "l"),
                         adjust = 1, Acol,
                         Ucol, Aborder, Uborder, Alwd = 1, Ulwd = 1, ...)
{
    ## verifications
    type <- match.arg(type)
    if (!inherits(x, "madifa"))
        stop("Object of class 'madifa' expected")

    ## which histogram should be drawn?
    if (scores)
        tab <- x$li
    else tab <- x$tab
    pr <- x$pr

    ## Graphical settings
    if (missing(Acol)) {
        Acol <- NULL
        Acolf <- "white"
        Acold <- "black"
    } else {
        Acold <- Acol
        Acolf <- Acol
    }
    if (missing(Aborder))
        Aborder <- "black"
    if (missing(Ucol)) {
        Ucol <- gray(0.8)
        Ucold <- gray(0.8)
    }
    else Ucold <- Ucol
    if (missing(Uborder))
        Uborder <- gray(0.8)
    clas <- rep("", ncol(tab))


    ## Quantitative or factor?
    for (j in 1:ncol(tab)) {
        w1 <- "q"
        if (is.factor(tab[, j]))
            w1 <- "f"
        clas[j] <- w1
    }
    if (any(clas == "f") & type == "l")
        warning("type = 'l' is not possible for factors, type = 'h' used instead.\n")

    ## Graphical settings, again
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(0.5, 0.5, 2, 0.5))
    par(mfrow = rev(n2mfrow(ncol(tab))))

    ## The function used for plotting each variable
    f1 <- function(j) {

        ## Use and availability
        tmpU <- rep(tab[, j], pr)
        tmpA <- tab[, j]
        name <- names(tab)[j]

        if (clas[j] == "f") {

            ## in case of factors
            par(mar = c(3, 0.5, 2, 0.5))
            mat <- t(cbind(table(tmpA), table(tmpU)))
            mat <- lapply(1:2, function(i) mat[i, ]/sum(mat[i,
                                                            ]))
            mat <- rbind(mat[[1]], mat[[2]])
            max <- max(mat)
            max <- max + max/20
            ylim <- c(0, max)
            barplot(mat, col = c(Acolf, Ucol), border = c(Aborder,
                                               Uborder), ylim = ylim, main = name, ylab = NULL,
                    axes = FALSE, beside = TRUE, ...)
            par(mar = c(0.5, 0.5, 2, 0.5))
        }
        else {

            ## In case of continuous variables

            if (type == "h") {

                ## If an histogram is desired
                xrange <- range(tmpA)
                H <- hist(tmpU, plot = FALSE, br = seq(min(xrange),
                                              max(xrange), length = 15))
                G <- hist(tmpA, plot = FALSE, br = seq(min(xrange),
                                              max(xrange), length = 15))
                yrange <- c(0, max(H$density, G$density))
                plot(H, freq = FALSE, col = Ucol, border = Uborder,
                     xlim = xrange, ylim = yrange, main = name,
                     xlab = NULL, ylab = "Density", axes = FALSE,
                     ...)
                plot(G, freq = FALSE, col = Acol, border = Aborder,
                     add = TRUE)
            }
            if (type == "l") {

                ## if a smoothing is wanted
                densA <- density(tmpA, adjust = adjust)
                densU <- density(tmpU, adjust = adjust, from = min(densA$x),
                                 to = max(densA$x))
                max <- max(densU$y, densA$y)
                max <- max + max/20
                ylim <- c(0, max)
                plot(densU, col = Ucol, ylim = ylim, type = "l",
                     lwd = Ulwd, main = name, xlab = NULL, ylab = "Density",
                     axes = FALSE, ...)
                lines(rep(mean(tmpU), 2), c(0, densU$y[512 -
                                                       sum(densU$x >
                                                           mean(tmpU))]),
                      col = Ucol, lty = 2,
                      lwd = Ulwd)
                lines(densA, col = Acold, lwd = Alwd)
                lines(rep(mean(tmpA), 2), c(0, densA$y[512 -
                                                       sum(densA$x >
                                                           mean(tmpA))]),
                      col = Acold, lty = 2,
                      lwd = Alwd)
            }
        }
        box()
    }

    ## applies f1 for each variable
    lapply(1:ncol(tab), f1)
    return(invisible(NULL))
}




predict.madifa <- function (object, map, nf, ...)
{
    ## Verifications
    if (!inherits(object, "madifa"))
      stop("object should be of class 'madifa'")
    if ((missing(nf)) || (nf > object$nf))
        nf <- object$nf

    if (!inherits(map, "SpatialPixelsDataFrame"))
        stop("should be an object of class SpatialPixelsDataFrame")
    gridded(map) <- TRUE
    gr <- gridparameters(map)
    if (nrow(gr) > 2)
        stop("map should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")

    ## sum of the squared scores of the variables
    ll <- apply(data.frame(object$l1[, 1:nf]),1,
                function(x) sum(object$cw * (x^2)))
    ll <- data.frame(Approx.MD=ll)
    coordinates(ll) <- coordinates(map)
    gridded(ll) <- TRUE
    return(ll)
  }



s.madifa <- function(x, xax=1, yax=2, cgrid = 1, clab=1, ...)
{
    ## Verifications
    if (!inherits(x, "madifa"))
        stop("Object of class 'madifa' expected")
    co <- x$co
    cw <- x$cw
    out <- seq(1,-1,length=200)

    ## graphical settings
    opar <- par(mar=c(0,0,0,0))
    on.exit(par(opar))

    ## trouvcoo finds the coordinates of the vertices of the ellipse
    trouvcoo <- function(co, cw, z, xax=1,yax=2)
    {
        x <- co[,xax]*sqrt(cw)
        y <- co[,yax]*sqrt(cw)
        mat <- rbind(z, sqrt(1-(z^2)))
        matb <- rbind(z[length(z):1], -sqrt(1-(z[length(z):1]^2)))
        mat <- cbind(mat,matb)
        u2 <- c(sum(x*y), sqrt(1 - ((sum(x*y))^2)))
        out <- cbind(c(z,z[length(z):1]), apply(mat,2,function(x) sum(x*u2)))
        return(out)
    }

    ## These coordinates are in yy
    yy <- trouvcoo(co, cw, out, xax, yax)

    ## The "outside" polygon, masking what is outside the ellipse
    pol <- data.frame(c(1,1,-1,-1,1,1),c(0,1,1,-1,-1,0))
    po <- rbind(as.matrix(pol),yy)

    ## Draws the plot
    s.arrow(co, xax=xax, yax=yax, xlim=c(-1,1), ylim=c(-1,1), clabel=clab)
    polygon(yy)
    polygon(po, col="white")
    s.arrow(co, xax=xax, yax=yax, xlim=c(-1,1), ylim=c(-1,1), add.plot=TRUE,
            clabel=clab, ...)

    ## scale box
    xaxp <- par("xaxp")
    ax <- (xaxp[2] - xaxp[1])/xaxp[3]
    yaxp <- par("yaxp")
    ay <- (yaxp[2] - yaxp[1])/yaxp[3]
    a <- min(ax, ay)
    cha <- paste(" d = ", a, " ", sep = "")
    cex0 <- par("cex") * cgrid
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    x1 <- par("usr")[2]
    y1 <- par("usr")[4]
    rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
    text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)

}




plot.madifa <- function(x, map, xax=1, yax=2, cont=FALSE,...)
{

    ## Verifications
    if (!inherits(x, "madifa"))
        stop("Object of class 'madifa' expected")

    if (!inherits(map, "SpatialPixelsDataFrame"))
        stop("should be an object of class SpatialPixelsDataFrame")
    gridded(map) <- TRUE
    gr <- gridparameters(map)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")

    ## Graphical settings
    opar <- par(mfrow=c(3,3))
    on.exit(par(opar))

    ## The eigenvaloue diagram
    scatterutil.eigen(x$eig, wsel = c(xax, yax))

    ## Column scores
    s.madifa(x, xax, yax, cgrid=2, clab=1.25)


    ## function to draw the niche
    foo <- function()
    {
        ## Some bases for the plot
        opar2 <- par(mar = c(0.1, 0.1, 0.1, 0.1))
        x1 <- x$l1[, xax]
        x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
        xlim <- range(x1)
        y1 <- x$l1[, yax]
        y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
        ylim <- range(y1)

        ## Background plot
        scatterutil.base(dfxy = x$l1[, c(xax, yax)], xax = 1, yax = 2,
                         xlim = xlim, ylim = ylim,
                         grid = TRUE, addaxes = FALSE,
                         cgrid = 2, include.origin = TRUE,
                         origin = c(0, 0), sub = "",
                         csub = 1.25, possub = "bottomleft",
                         pixmap = NULL, contour = NULL,
                         area = NULL, add.plot = FALSE)

        ## Graphical settings
        Acol <- gray(0.8)
        Ucol <- "black"
        Abg <- gray(0.8)
        Ubg <- "black"
        Ainch <- 0.03
        Uinch <- Ainch * max(x$pr)

        ## adds available and used points
        symbols(x$l1[, c(xax, yax)], circles = rep(1, length(x$pr)),
                fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
        symbols(x$l1[x$pr > 0, c(xax, yax)], circles = x$pr[x$pr > 0],
                fg = Ucol, bg = Ubg,
                inches = Uinch, add = TRUE)

        ## axes
        abline(v = 0)
        abline(h = 0)
        box()
        par(opar2)
    }
    foo()

    ## the maps
    ka <- data.frame(Maha=x$mahasu, mod=apply(x$l1[,c(xax,yax)],1,
                                    function(x) sqrt(sum(x^2))),
                     xa=x$l1[,xax], ya=x$l1[,yax])
    coordinates(ka) <- coordinates(map)
    gridded(ka) <- TRUE

    u <- par(mar=c(0.1,0.1,2,0.1))
    image(ka, 1)
    title(main="Mahalanobis distances")
    if (cont)
        contour(ka,1, add=TRUE)
    box()
    image(ka,2, axes=FALSE)
    title(main="from the analysis")
    if (cont)
        contour(ka,2,add=TRUE)
    box()
    par(u)

    ## Correlation with the environmental variables
    s.arrow(x$cor, xax=xax,yax=yax,
            sub="Cor(habitat var., scores) available",
            clabel=1.25, csub=2, cgrid=2, xlim=c(-1,1), ylim=c(-1,1))
    u <- par(mar=c(0.1,0.1,2,0.1))

    ## Again the maps
    image(ka,3,axes=FALSE)
    title(main="Axis 1")
    if (cont)
        contour(ka,3, add=TRUE)
    box()
    image(ka,4, axes=FALSE)
    title(main="Axis 2")
    if (cont)
        contour(ka,4, add=TRUE)
    box()
    par(u)
}
