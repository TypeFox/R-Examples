"tab2by2.test" <-
function(x, y = NULL,
         correction = FALSE,
         rev = c("neither", "rows", "columns", "both")
         ){
  if(is.matrix(x) && !is.null(y)){stop("y argument should be NULL")}
  if(is.null(y)){
    x <- epitable(x, rev = rev)
  } else {
    x <- epitable(x, y, rev = rev)
  }
  nr <- nrow(x); nc <- ncol(x)
  fish <- chi2 <- midp <- rep(NA, nr)  
  for(i in 2:nr){
    xx <- x[c(1,i),]
    a0<-x[1,2]; b0<-x[1,1]; a1<-x[i,2]; b1<-x[i,1]
    fish[i] <- fisher.test(xx)$p.value    
    chi2[i] <- chisq.test(xx, correct = correction)$p.value
    midp[i] <- ormidp.test(a1, a0, b1, b0)$two.sided
  }
  pv <- cbind(midp, fish, chi2)
  colnames(pv) <- c("midp.exact", "fisher.exact", "chi.square")
  rownames(pv) <- rownames(x)
  names(dimnames(pv)) <- c(names(dimnames(x))[1], "two-sided")
  list(x = x,
       p.value = pv,
       correction = correction
       )
}
