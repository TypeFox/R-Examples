estimateContrasts2 <- function (contrast.matrix, fit, alpha = 0.05, row = TRUE,L) 
{
    if (!inherits(fit, "lm")) 
        stop("Second input is not an \"lm\" object")
    if (length(dimnames(fit$model)[[2]]) != 3) 
        stop("estimateContrasts2 requires 2 factors in \"lm\" model!")
    if (!is.factor(fit$model[[2]])) 
        stop("1st explanatory variable in \"lm\" not a factor")
    if (!is.factor(fit$model[[3]])) 
        stop("2nd explanatory variable in \"lm\" not a factor")
	if (alpha < 0 || alpha > 1) 
	alpha <- 0.05    

 inter <- ifelse(  (nrow(anova(fit)) == 4)   ,TRUE,FALSE)
 dc <- dummy.coef(fit)
 n1 <- length(dc[[2]])
 n2 <- length(dc[[3]])
 if (!inter) dc[[4]] <- rep(0,n1*n2)
 	cellmeans <- (matrix(dc[[4]],n1,n2) + outer(dc[[2]],rep(1,n2)) + outer(rep(1,n1),dc[[3]])) + dc[[1]]
 	means <- if (row)		
			apply(cellmeans,1,mean)
           	 else   apply(cellmeans,2,mean)
	contrast.matrix <- if (!is.matrix(contrast.matrix)) 
		matrix(contrast.matrix, ncol = length(contrast.matrix))
	else contrast.matrix
	if (ncol(contrast.matrix) != length(means)) 
		stop("Contrast matrix has wrong number of columns")
	k <- length(fit$model[,1])/(n1 * n2)  # ONLY WORKS FOR BALANCED DATA
	groupsize <- ifelse(row,n2*k,n1*k) 
	contrasts <- contrast.matrix %*% means
      if (is.null(L))
    		L <- length(contrasts)
	stderr <- sqrt(contrast.matrix^2 %*% (1/rep(groupsize, length(means)))) * summary(fit)$sigma
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
		paste("C", 1:length(contrasts), sep = "")
	else dimnames(contrast.matrix)[[1]]   
	# drop unadjusted probs 2005          
	#dimnames(outmat) <- list(contrast.names,c("Estimate","Tukey.L","Tukey.U","Tukey.p","Unadj.p"))
	dimnames(outmat) <- list(contrast.names,c("Estimate","Tukey.L","Tukey.U","Tukey.p"))
	outmat
}

