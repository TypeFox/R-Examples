plot.cv.npmr <-
function(x, feature.names = TRUE, ...) {
    plot(x$fit, x$lambda.min, feature.names)
}
