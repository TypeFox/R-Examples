predict.cv.npmr <-
function(object, newx, ...) {
    predict(object$fit, newx)[,,1]
}
