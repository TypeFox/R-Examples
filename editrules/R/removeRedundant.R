# 
removeRedundant <- function(E){
  
  m <- as.mip(E)
  A <- -getA(m$E)
  b <- getb(m$E)
  keep <- rep(TRUE, nrow(A))
  names(keep) <- rownames(A)
  
  sapply(seq_len(nrow(A)),function(r){
    m1 <- m
    #keep[r] <<- FALSE
    m1$E <- m1$E[keep,]
    
    m1$objfn <- A[r,]
    
    lps <- as.lp.mip(m1)
    statuscode <- solve(lps)
    
    o <- -1*get.objective(lps)
    
    keep[r] <<- (b[r] <= o) 
    c(o=o, t=b[r])
  })
  structure(E[names(keep)[keep],], removed=E[names(keep)[!keep],])
}

# E <- editmatrix(c(A="x>2", B="y > x", C="y>1" , D="x>1"))
# Er <- removeRedundant(E)
# Er
