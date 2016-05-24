#' @title Adaptation by Sanchez-Aluja of Lebart test
#' 
#' @description
#' Internal function. \code{test_sanchez_aluja} is called by \code{pathmox}
#' 
#' @param Y.lvs matrix with scores of latent variables
#' @param IDM inner design matrix
#' @param split.a split a (left side)
#' @param split.b split b (right side)
#' @export
#' @keywords internal
test_sanchez_aluja <- function(Y.lvs, IDM, split.a, split.b)
{    
  lvs <- ncol(Y.lvs)  
  endo <- rowSums(IDM)  	
  endo[endo!=0] <- 1  
  n.endo <- sum(endo)
  n <- nrow(Y.lvs)
  p.values <- rep(0,n.endo)		
  f.values <- rep(0,n.endo)
  df1.values <- rep(0,n.endo)
  df2.values <- rep(0,n.endo)
  s1s2.values <- rep(0,n.endo)
  who.endo <- which(endo == 1)
  SCE0 <- 0
  SCE1 <- 0		
  Sa <- 0
  Sb <- 0
  p.sum <- 0
  
  for (i in 1:n.endo)
  {
    aux <- who.endo[i]
    indep <- which(IDM[aux,1:aux] == 1)
    reg <- cbind(Y.lvs[,aux], Y.lvs[,indep])
    # SS0 calculation
    lm.res <- lm(Y.lvs[,aux] ~ Y.lvs[,indep])$residuals
    SS0 <- sum(lm.res^2)
    SCE0 <- SCE0 + SS0
    # regression matrices "Xa" and "Xb"
    Xa <- reg[split.a,]
    Xb <- reg[split.b,]
    lm.resa <- lm(Xa[,1] ~ Xa[,2:ncol(Xa)])$residuals
    lm.resb <- lm(Xb[,1] ~ Xb[,2:ncol(Xb)])$residuals
    SS1 <- sum(lm.resa^2) + sum(lm.resb^2)
    SCE1 <- SCE1 + SS1
    Sa <- Sa + lm.resa
    Sb <- Sb + lm.resb
    # test
    p <- length(indep)        
    p.sum <- p.sum + p
    df1 <- p
    df2 <- n - 2*p
    f <- ((n-2*p)/p) * ((SS0-SS1)/SS1)
    p.values[i] <- 1 - pf(f, df1, df2)
    f.values[i] <- f
    df1.values[i] <- df1
    df2.values[i] <- df2
    s1s2.values[i] <- round(sd(lm.resa) / sd(lm.resb), 4)
  }
  
  n.sum <- n.endo * (length(split.a) + length(split.b))
  # f global
  f.sum <- ((n.sum-2*p.sum)/p.sum)  *  ((SCE0 - SCE1)/SCE1) 
  df1.sum <- p.sum   						# gl num
  df2.sum <- n.sum - 2*p.sum					# gl den
  pval.sum <- 1 - pf(f.sum, df1.sum, df2.sum)	   		# pv
  s1s2.sum <- sd(Sa) / sd(Sb)
  F.global <- c(f.sum, df1.sum, df2.sum, pval.sum, s1s2.sum)   		
  
  F.partial <- cbind(f.values, df1.values, df2.values, p.values, s1s2.values) 
  rownames(F.partial) <- rownames(IDM)[who.endo]
  colnames(F.partial) <- c("f.stat", "df.num", "df.den", "p.val", "sa/sb")
  res <- list(F.global, F.partial)
  return(res)
}
