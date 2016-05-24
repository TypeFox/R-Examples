score.calc <- function(marks,y,Z,X,K,M,Hinv,ploidy,model,min.MAF,max.geno.freq) {
  # needs to be modified to work if Hinv is not provided, not in a hurry, normally wouldn't be possible
  #
  m <- length(marks)
  n <- nrow(Z)
  scores <- numeric(m)*NA
  beta.out <- scores
  general <- length(grep("general",model,fixed=T))>0 
  ####################
  ## initialize the progress bar
  count <- 0
  tot <- m
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  #####################
  for (i in 1:m) {
    ###################
    count <- count + 1
    ###################
    if(ploidy == 2){
      S <- as.matrix(M[,marks[i]])
    }else{
      S <- design.score(M[,marks[i]],model,ploidy,min.MAF,max.geno.freq)
    }
    
    if (!is.null(S)) {
      v1 <- ncol(S)
      X2 <- cbind(X,Z%*%S)
      p <- ncol(X2)
      v2 <- n - p                 
      if (is.null(Hinv)) {			
        out <- try(mmer(y=y,X=X2,Z=list(list(Z=Z,K=K))),silent=TRUE)
        if (class(out)!="try-error") { 
          Hinv <- out$Hinv 
        }
      } 
      W <- crossprod(X2, Hinv %*% X2) 
      Winv <- try(solve(W),silent=TRUE)
      if (class(Winv) != "try-error") {
        beta <- Winv %*% crossprod(X2, Hinv %*% y)
        resid <- y - X2 %*% beta
        s2 <- as.double(crossprod(resid, Hinv %*% resid))/v2
        Q <- s2 * Winv[(p-v1+1):p,(p-v1+1):p]
        Tt <- solve(Q, silent= TRUE)
        if (class(Tt) != "try-error") {
          V <- beta[(p+1-v1):p]
          Fstat <- crossprod(V,Tt%*%V)/v1
          x <- v2/(v2+v1*Fstat)
          scores[i] <- -log10(pbeta(x, v2/2, v1/2)) 
          if (!general) {beta.out[i] <- beta[p]}                    
        }
      }
    }
    ################################
    setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    ################################
  }
  return(list(score=scores,beta=beta.out))			
}
