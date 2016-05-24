`SS.calib` <- structure(function (Blev, k, gamm, GL) 
{
    .expr1 <- GL^gamm
    .expr4 <- ifelse(GL == 0, 0, log(GL))
    .expr5 <- .expr1 * .expr4
    .value <- Blev + k * .expr1
    .grad <- array(0, c(length(.value), 3), list(NULL, c("Blev", 
        "k", "gamm")))
    .hessian <- array(0, c(length(.value), 3, 3), list(NULL, 
        c("Blev", "k", "gamm"), c("Blev", "k", "gamm")))
    .grad[, "Blev"] <- 1
    .grad[, "k"] <- .expr1
    .hessian[, "k", "k"] <- 0
    .hessian[, "k", "gamm"] <- .hessian[, "gamm", "k"] <- .expr5
    .grad[, "gamm"] <- k * .expr5
    .hessian[, "gamm", "gamm"] <- k * (.expr5 * .expr4)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
, initial = function(mCall, data, LHS) {
	xy <- sortedXyData(mCall[["GL"]], LHS, data)
	Blev <- min(xy[["y"]])
	k <- max(xy[["y"]])
	gamm <- 2.5
	value <- c(Blev, k, gamm)
	names(value) <- mCall[c("Blev", "k", "gamm")]
	value
	}
, pnames = c("Blev", "k", "gamm"), class = "selfStart")
