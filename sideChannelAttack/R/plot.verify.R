plot.verify.loo <-
function (x, ...) 
{
    object = x
    Y = (object$TP + object$TN)/(object$TP + object$TN + object$FP + 
        object$FN)
    X = object$dim
    plot(x = X, y = Y, type = "l", main = "% of good answers", 
        xlab = "number of features", ylab = "success rate")
    points(x = X, y = Y)
}
plot.verify.ho <-
function (x, ...) 
{
    object = x
    Y = (object$TP + object$TN)/(object$TP + object$TN + object$FP + 
								 object$FN)
    X = object$dim
    plot(x = X, y = Y, type = "l", main = "% of good answers", 
		 xlab = "number of features", ylab = "success rate")
    points(x = X, y = Y)
}
plot.verify.cv <-
function (x, ...) 
{
    object = x
    Y = (object$TP + object$TN)/(object$TP + object$TN + object$FP + 
        object$FN)
    X = object$dim
    plot(x = X, y = Y, type = "l", main = "% of good answers", 
        xlab = "number of features", ylab = "success rate")
    points(x = X, y = Y)
}