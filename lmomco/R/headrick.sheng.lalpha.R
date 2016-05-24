"headrick.sheng.lalpha" <- function(x, ...) {
   ifelse(is.matrix(x), lco2 <- x,
                        lco2 <- Lcomoment.matrix(x, k=2)$matrix)
   dd <- dim(lco2)
   if(dd[1] != dd[2]) {
      warning("not a square matrix, returning NA")
      return(NA)
   }
   d <- dd[1]
   A <- sum(diag(as.matrix(lco2)))
   B <- sum(sapply(1:d, function(i) {
        sum(sapply(1:d, function(j) {
            ifelse(i == j, 0, lco2[i,j])
        } ))
   } ))
   alpha <- (d/(d-1)) * (1 - A/(A+B))
   zz <- list(alpha=alpha,
              title="Headrick and Shengs' L-alpha",
              source="headrick.sheng.alpha")
   return(zz)
}

"lalpha" <- function(x, ...) {
   return(headrick.sheng.lalpha(x))
}




