##' Illustrate multicollinearity in regression, part 1.
##'
##' @description
##' This is a set of functions that faciliates the examination
##' of multicollinearity. Suppose the "true" relationship is
##'  y[i] = 0.2 * x1[i] + 0.2 * x2[i] + e
##' where e is Normal(0, stde^2).
##'
##' mcGraph1 draws the 3D regression space, but all of the points
##' are illustrated "in" the flat plane of x1-x2 variables.
##'
##' @details
##' These functions are specialized to a particular purpose. If you
##' just want to draw 3D regressions, look at plotPlane.
##' @name mcGraph1
##' @param x1 a predictor vector
##' @param x2 a predictor vector
##' @param y  the dependent variable
##' @param x1lab label for the x1 axis, (the one called "xlab" inside persp)
##' @param x2lab label for the x2 axis, (the one called "ylab" inside persp)
##' @param ylab  label for the y (vertical) axis (the one called "zlab" inside persp)
##' @param ... additional parameters passed to persp
##' @export
##' @rdname mcGraph
##' @return mcGraph1 and mcGraph2 return only the perspective matrix
##' from persp. It is returned so that users can add additional embellishments on the 3D plot (can be used with trans3d)
##' @keywords regression hplot
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @examples
##' set.seed(12345)
##' ## Create data with x1 and x2 correlated at 0.10
##' dat <- genCorrelatedData(rho=.1, stde=7)
##'
##' mcGraph1(dat$x1, dat$x2, dat$y, theta=20, phi=8, ticktype="detailed", nticks=10)
##'
mcGraph1 <- function (x1, x2, y, x1lab, x2lab, ylab, ...){
    x1range <- magRange(x1, 1.25)
    x2range <- magRange(x2, 1.25)
    yrange <- magRange(y, 1.5)
    
    if (missing(x1lab)) x1lab <- gsub(".*\\$", "", deparse(substitute(x1)))
    if (missing(x2lab)) x2lab <- gsub(".*\\$", "", deparse(substitute(x2)))
    if (missing(ylab)) ylab  <- gsub(".*\\$", "", deparse(substitute(y)))
    
    res <- perspEmpty(x1 = plotSeq(x1range, 5),
                      x2 = plotSeq(x2range,5),
                      y = yrange, x1lab = x1lab,
                      x2lab = x2lab, ylab = ylab, ...)
    
    yMinimum <- rep(yrange[1] , length(x1))
    mypoints1 <- trans3d(x1, x2, yMinimum, pmat = res)
    points( mypoints1, pch = 16, col = "blue")
    invisible(res)
}
NULL

 
##' mcGraph2 draws a 3-D representation of a scatterplot with shadows
##' in the x1-x2 plane.  The observations are represented by blue
##' points floating above the x1-x2 plane. If scaley=1, the end result
##' is a scatterplot "cloud" of the y points above the x1-x2 plane,
##' and gray shadows of the points are cast down from the cloud onto
##' the x1-x2 plane itself. This uses persp to make the actual
##' drawing.
##'
##' @param rescaley a single scalar value or a vector of the same
##' length as y.
##' @export
##' @rdname mcGraph
##' @examples
##' set.seed(12345)
##' ## Create data with x1 and x2 correlated at 0.10
##' dat <- genCorrelatedData(rho=.1, stde=7)
##' ## This will "grow" the "cloud" of points up from the
##' ## x1-x2 axis
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.0, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.1, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.2, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.3, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.4, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.5, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.6, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.7, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.8, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 0.9, theta = 0)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 1, theta = 0)
##'
##' ##rotate this
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 1, theta = 20)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 1, theta = 40)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 1, theta = 60)
##' mcGraph2(dat$x1, dat$x2, dat$y, rescaley = 1, theta = 80)
##'
##' ## once they reach the top, make them glitter a while
##' for(i in 1:20){
##'   mcGraph2(dat$x1, dat$x2, dat$y, rescaley = runif(length(dat$x1), .9,1.1), theta = 0)
##' }
##'
mcGraph2 <- function(x1, x2, y, rescaley = 1, drawArrows = TRUE, x1lab, x2lab, ylab, ...){
    x1range <- magRange(x1, 1.25)
    x2range <- magRange(x2, 1.25)
    yrange <- magRange(y, 1.5)

    if (missing(x1lab)) x1lab <- gsub(".*\\$", "", deparse(substitute(x1)))
    if (missing(x2lab)) x2lab <- gsub(".*\\$", "", deparse(substitute(x2)))
    if (missing(ylab)) ylab  <- gsub(".*\\$", "", deparse(substitute(y)))

    res <- perspEmpty(x1 = plotSeq(x1range,5), x2 = plotSeq(x2range,5), y = yrange, x1lab = x1lab, x2lab = x2lab, ylab = ylab, ...)

    mypoints1 <- trans3d ( x1, x2 ,yrange[1], pmat = res )
    newy <- rescaley * (y - yrange[1]) + yrange[1]
    mypoints2 <- trans3d ( x1 , x2 , newy , pmat = res )

    points( mypoints1, pch = 16, col=gray(0.8))
    points( mypoints2, pch = 1, col= "blue")
    mypoints2s <- trans3d(x1, x2, (0.8)*newy, pmat = res)
    if (drawArrows) arrows(mypoints1$x , mypoints1$y , mypoints2s$x , mypoints2s$y , col="red" , lty = 2, lwd = 0.3, length = 0.1)
    invisible(res)
}
NULL

##' With mcGraph3, the observations are scattered in 3-dimensions, the
##' fitted values are represented by a mesh, and their shadows in the
##' x1-x2 plane are also represented.
##'
##' @param interaction a TRUE or FALSE request for inclusion of the x1-x2 interaction in the regression calculation
##' @param drawArrows TRUE or FALSE, do you want arrows from the plane to observed y?
##' @return mcGraph3 returns a list of 2 objects. 1) the fitted
##' regression model2) the perspective matrix used with persp to draw
##' the image.
##' @export
##' @rdname mcGraph
##' @examples
##' set.seed(12345)
##' ## Create data with x1 and x2 correlated at 0.10
##' dat <- genCorrelatedData(rho=.1, stde=7)
##'
##' mcGraph3(dat$x1, dat$x2, dat$y, theta = 0)
##'
##' dat2 <- genCorrelatedData(rho = 0, stde = 7)
##'
##' mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = 0, phi = 10)
##' mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = 30, phi = 10)
##' mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = -30, phi = 10)
##' mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = -30, phi = -10)
##' mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = -30, phi = -15)
##'
##' ## Run regressions with not-strongly correlated data
##' modset1 <- list()
##' for(i in 1:20){
##'   dat2 <- genCorrelatedData(rho = .1, stde = 7)
##'   summary(lm( y ~ x1 + x2 , data = dat2))
##'   modset1[[i]] <- mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = -30)
##' }
##'
##'
##' ## Run regressions with strongly correlated data
##' modset2 <- list()
##' for(i in 1:20){
##'   dat2 <- genCorrelatedData(rho = .981, stde = 7)
##'   summary(lm( y ~ x1 + x2 , data = dat2))
##'   modset2[[i]] <- mcGraph3(dat2$x1, dat2$x2, dat2$y, theta = -30)
##' }
##'
##' dat3 <- genCorrelatedData(rho = .981, stde = 100, beta=c(0.1, 0.2, 0.3, -0.1))
##' mcGraph3(dat3$x1, dat3$x2, dat3$y, theta=-10, interaction = TRUE)
mcGraph3 <- function(x1, x2, y, interaction = FALSE, drawArrows = TRUE, x1lab, x2lab, ylab, ...){
    x1range <- magRange(x1, 1.25)
    x2range <- magRange(x2, 1.25)
    yrange <- magRange(y, 1.5)

    if (missing(x1lab)) x1lab <- gsub(".*\\$", "", deparse(substitute(x1)))
    if (missing(x2lab))  x2lab <- gsub(".*\\$", "", deparse(substitute(x2)))
    if (missing(ylab)) ylab  <- gsub(".*\\$", "", deparse(substitute(y)))

    res <- perspEmpty(x1 = plotSeq(x1range, 5), x2 = plotSeq(x2range, 5),  y = yrange, x1lab = x1lab, x2lab = x2lab, ylab = ylab, ...)

    mypoints1 <- trans3d( x1, x2, yrange[1], pmat = res )
    points(mypoints1, pch = 16, col = gray(0.8))

    mypoints2 <- trans3d( x1, x2, y, pmat = res )
    points(mypoints2, pch = 1, col = "blue")

    if (interaction) m1 <- lm(y ~ x1 * x2)
    else m1 <- lm(y ~ x1 + x2)

    x1seq <- plotSeq (x1range, length.out = 20)
    x2seq <- plotSeq (x2range, length.out = 20)

    zplane <- outer (x1seq, x2seq, function(a, b) { predict(m1,
                      newdata = data.frame(x1 = a, x2 = b))})

    for( i in 1:length(x1seq) ){
        lines(trans3d(x1seq[i], x2seq, zplane[i,], pmat = res), lwd = 0.3)
    }
    for( j in 1:length(x2seq) ){
        lines(trans3d(x1seq, x2seq[j], zplane[,j], pmat = res), lwd = 0.3)
    }

    mypoints4 <- trans3d (x1 , x2 , fitted(m1) , pmat = res)
    ##  points(mypoints4)
    newy <- ifelse(fitted(m1) < y, fitted(m1)+ 0.8*(y-fitted(m1)),
                   fitted(m1) + 0.8 * (y-fitted(m1)))
    mypoints2s <- trans3d(x1, x2, newy, pmat = res)

    if (drawArrows) arrows(mypoints4$x , mypoints4$y , mypoints2s$x , mypoints2s$y , col = "red" , lty = 4, lwd = 0.3, length = 0.1)
    invisible(list(lm = m1, res = res))
}
NULL

