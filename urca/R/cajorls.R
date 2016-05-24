cajorls <- function(z, r = 1, reg.number = NULL){
    if (!(class(z) == "ca.jo") && !(class(z) == "cajo.test")) {
        stop("\nPlease, provide object of class 'ca.jo' or 'cajo.test' as 'z'.\n")
    }
    P <- ncol(z@Z0)
    r <- as.integer(r)
    if((r < 1) || (r > P - 1)){
      stop(paste("Please, provide a cointegration rank 'r' in the interval of 1 to", P - 1, ", in accordance with 'z@V'.\n", sep = " "))
    }
    beta <- matrix(z@V[, 1:r], ncol = r)
    C1 <- diag(r)
    C2 <- matrix(0, nrow = nrow(beta) - r, ncol = r)
    C <- rbind(C1, C2)
    betanorm <- beta %*% solve(t(C) %*% beta)
    ECT <- z@ZK %*% betanorm
    colnames(ECT) <- paste("ect", 1:r, sep = "")
    colnames(betanorm) <- colnames(ECT)
    rownames(betanorm) <- colnames(z@ZK)
    data.mat <- data.frame(z@Z0, ECT, z@Z1)
    text <- colnames(data.mat)[-c(1:P)]
    text1 <- paste(text, "", sep = "+", collapse = "")
    text2 <- paste("~", substr(text1, 1, nchar(text1) - 1))
    if (!is.null(reg.number)) {
        reg.number <- as.integer(reg.number)
        if (reg.number > ncol(z@Z0) || reg.number < 1) {
            stop("\nPlease, provide a valid number of the regression within \n the VECM, numbering from 1 to ", P, ".\n")
        }
        form1 <- formula(paste("z@Z0[, reg.number]", text2, "-1"))
        rlm <- lm(substitute(form1), data = data.mat)
    }
    else if (is.null(reg.number)) {
        form1 <- formula(paste("z@Z0", text2, "-1"))
        rlm <- lm(substitute(form1), data = data.mat)
    }
    return(list(rlm = rlm, beta = betanorm))
}
