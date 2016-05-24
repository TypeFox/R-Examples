ordGene <- function(xpr, lvs, type = c("RLRT", "LRT"), nsim = 1e6,
null.sample = NULL, progressBar = TRUE, ...){

  type <- match.arg(type)
  type <- switch(type, RLRT="RLRT", LRT="LRT")

  ## Check xpr and lvs
  if(!(is.matrix(xpr) | is.data.frame(xpr)))
    stop("xpr has to be a matrix or data frame")

  if(length(ncol(lvs)) > 0 | !is.numeric(lvs))
    stop("lvs has to be a numeric vector")

  if(any(is.na(xpr)))
    stop("Missing values in xpr are not allowed")

  if(any(is.na(lvs)))
    stop("Missing values in lvs are not allowed")

  # ordAOV null distribution
  x <- as.numeric(factor(lvs))
  y <- as.numeric(xpr[1,])
  if (length(null.sample)==0)
    RLRTnull <- ordAOV(x,y, type=type, nsim=nsim, ...)$sample
  else
    RLRTnull <- null.sample

  # ordAOV
  n <- nrow(xpr)
  pRLRT <- numeric(n)
  pANOVA <- numeric(n)
  pttest <- numeric(n)
  if (progressBar == TRUE) 
    {
      pb <- tkProgressBar(title = "progress bar", min = 0, max = n, width = 300)
    }

  for (j in 1:n)
    {
      y <- as.numeric(xpr[j,])
      pRLRT[j] <- ordAOV(x,y, type=type, null.sample=RLRTnull, ...)$p.value
      pANOVA[j] <- anova(lm(y ~ factor(x)))$"Pr(>F)"[1]
      pttest[j] <- anova(lm(y ~ x))$"Pr(>F)"[1]
      if (progressBar == TRUE) 
        {
          setTkProgressBar(pb, j, title = paste(round(j/n * 100, 0), "% done"))
        }
    }

  if (progressBar == TRUE) 
    close(pb)
    
  pvals <- cbind(pRLRT,pANOVA,pttest)
  rownames(pvals) <- row.names(xpr)
  if (type=="RLRT")
    colnames(pvals) <- c("RLRT", "ANOVA", "t-test")
  else
    colnames(pvals) <- c("LRT", "ANOVA", "t-test")

  # out
  return(pvals)
}
