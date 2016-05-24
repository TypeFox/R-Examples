`SS.RGBcalib` <-structure(function (Blev, Br, Bg, Bb, gamm, Rgun, Ggun, Bgun) 
{
    .expr1 <- Rgun^gamm
    .expr4 <- Ggun^gamm
    .expr7 <- Bgun^gamm
    .expr10 <- ifelse(Rgun == 0, 0, log(Rgun))
    .expr11 <- .expr1 * .expr10
    .expr12 <- ifelse(Ggun == 0, 0, log(Ggun))
    .expr13 <- .expr4 * .expr12
    .expr14 <- ifelse(Bgun == 0, 0, log(Bgun))
    .expr15 <- .expr7 * .expr14
    .value <- Blev + Br * .expr1 + Bg * .expr4 + Bb * .expr7
    .grad <- array(0, c(length(.value), 5), list(NULL, c("Blev", 
        "Br", "Bg", "Bb", "gamm")))
    .hessian <- array(0, c(length(.value), 5, 5), list(NULL, 
        c("Blev", "Br", "Bg", "Bb", "gamm"), c("Blev", "Br", 
            "Bg", "Bb", "gamm")))
    .grad[, "Blev"] <- 1
    .grad[, "Br"] <- .expr1
    .hessian[, "Br", "Br"] <- 0
    .hessian[, "Br", "Bg"] <- .hessian[, "Bg", "Br"] <- 0
    .hessian[, "Br", "Bb"] <- .hessian[, "Bb", "Br"] <- 0
    .hessian[, "Br", "gamm"] <- .hessian[, "gamm", "Br"] <- .expr11
    .grad[, "Bg"] <- .expr4
    .hessian[, "Bg", "Bg"] <- 0
    .hessian[, "Bg", "Bb"] <- .hessian[, "Bb", "Bg"] <- 0
    .hessian[, "Bg", "gamm"] <- .hessian[, "gamm", "Bg"] <- .expr13
    .grad[, "Bb"] <- .expr7
    .hessian[, "Bb", "Bb"] <- 0
    .hessian[, "Bb", "gamm"] <- .hessian[, "gamm", "Bb"] <- .expr15
    .grad[, "gamm"] <- Br * .expr11 + Bg * .expr13 + Bb * .expr15
    .hessian[, "gamm", "gamm"] <- Br * (.expr11 * .expr10) + 
        Bg * (.expr13 * .expr12) + Bb * (.expr15 * .expr14)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
, initial = function(mCall, data, LHS) {
	Lum <- eval(asOneSidedFormula("Lum")[[2]], data)
	Rgun <- eval(asOneSidedFormula("Rgun")[[2]], data)
	Ggun <- eval(asOneSidedFormula("Ggun")[[2]], data)
	Bgun <- eval(asOneSidedFormula("Bgun")[[2]], data)
	Blev <-  min(Lum)
	Br <- max(Lum[Rgun > 0])
	Bg <- max(Lum[Ggun > 0])
	Bb <- max(Lum[Bgun > 0])
	gamm <- 2.5
	value <- c(Blev, Br, Bg, Bb, gamm)
	names(value) <- mCall[c("Blev", "Br", "Bg", "Bb", "gamm")]
	value
	}
, pnames = c("Blev", "Br", "Bg", "Bb", "gamm"), class = "selfStart")
