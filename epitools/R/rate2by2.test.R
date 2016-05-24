"rate2by2.test" <-
function(x, y = NULL, rr = 1, 
         rev = c("neither", "rows", "columns", "both")
         ){
  if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
  if(is.null(y)){
    x <- ratetable(x, rev = rev)
  } else {
    xn <- substitute(x)
    yn <- substitute(y)
    x <- ratetable(x, y, rev = rev)
    colnames(x) <- c(xn, yn)
  }
  tmx <- table.margins(x)[,-3]
  nr <- nrow(x)
  p.value <- matrix(NA, nr, 2)
  for(i in 2:nr){
    aa <- x[i,1]; bb <- x[1,1]; pt1 <- x[i,2]; pt0 <- x[1,2]
    pt <- pt0 + pt1
    mm <- aa + bb
    s <- rr*pt1/(rr*pt1 + pt0)
    p.lower <- dbinom(aa, mm, s)/2 + pbinom(aa-1, mm, s)
    p.upper <- 1 - p.lower
    pval1 <- min(p.lower, p.upper)
    pval2 <- 2*pval1
    ##Score p value
    num <- aa - (pt1/pt)*mm
    dem <- sqrt(mm*(pt1/pt)*(pt0/pt))
    zval <- num/dem
    chi2 <- (num/dem)^2
    pv <- 1-pnorm(abs(zval))
    pv2 <- 1-pchisq(chi2, df=1)
    p.value[i,] <- c(pval2, pv2)
  }
  colnames(p.value) <- c("midp.exact", "wald")
  rownames(p.value) <- rownames(x)  
  if(is.null(names(dimnames(x)))){
    names(dimnames(p.value)) <- c("Predictor", "Outcome")
  }
  if(!is.null(names(dimnames(x)))){
    names(dimnames(p.value)) <- c(names(dimnames(x))[1], "two-sided")
  }
  if(rr!=1) {p.value <- p.value[,"midp.exact"]}
  rrl <- list(x = x,
              p.value = p.value
              )
  rrl
}
