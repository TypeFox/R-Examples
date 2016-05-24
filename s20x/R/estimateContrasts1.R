estimateContrasts1 <- function (contrast.matrix, fit, alpha = 0.05,L) 
{
    if (!inherits(fit, "lm")) 
        stop("Second input is not an \"lm\" object")
    if (length(dimnames(fit$model)[[2]]) != 2) 
        stop("Require only one factor in \"lm\" model!")
    else if (!is.factor(fit$model[[2]])) 
        stop("Explanatory variable in \"lm\" should be a factor!")
    if (alpha < 0 || alpha > 1) 
        alpha <- 0.05
     k <- fit$rank
    contrast.matrix <- if (!is.matrix(contrast.matrix)) 
        matrix(contrast.matrix, ncol = length(contrast.matrix))
    else contrast.matrix
    if (ncol(contrast.matrix) != k) 
        stop("Wrong number of columns in contrast matrix")
    y <- fit$model[, 1]
    f1 <- fit$model[, 2]
    means <- as.vector(tapply(y, f1, mean))
    lengths <- as.vector(tapply(y, f1, length))
    rms <- summary(fit)$sigma
    contrasts <- contrast.matrix %*% means
    if (is.null(L))
    	L <- length(contrasts)
    qqt <- qt(1 - alpha/(2 * length(contrasts)), fit$df)
    stderr <- sqrt(contrast.matrix^2 %*% (1/lengths)) * rms
    tstat <- round(contrasts/stderr,4)
    tquant <- qtukey(1 - alpha,ncol(contrast.matrix), fit$df)
    tl <- round(contrasts - stderr * tquant/sqrt(2), 4)
    tu <- round(contrasts + stderr * tquant/sqrt(2), 4)
    unadj.p <- round(2*(1 - pt(abs(tstat),fit$df)),4) 
    tukey.p <- 1-round((ptukey(abs(tstat)*sqrt(2),ncol(contrast.matrix),fit$df)),4)
    bonf.p <- pmin(round(L* 2*(1 - pt(abs(tstat),fit$df)), 4),1)
   # drop unadjusted probs 2005
   #outmat <- cbind(contrasts, tl, tu, tukey.p, unadj.p)
    outmat<-cbind(contrasts, tl, tu, tukey.p)
    contrast.names <- if (is.null(dimnames(contrast.matrix))) 
        paste("Contrast", 1:length(contrasts))
    else dimnames(contrast.matrix)[[1]]
    # drop unadjusted probs 2005
    #dimnames(outmat) <- list(contrast.names,c("Estimate","Tukey.L","Tukey.U","Tukey.p","Unadj.p"))
    dimnames(outmat) <- list(contrast.names,c("Estimate","Tukey.L","Tukey.U","Tukey.p"))
	outmat
}

