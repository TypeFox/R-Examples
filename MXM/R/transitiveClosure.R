## transitive closure

transitiveClosure <- function(amat) {
  
  n <- nrow(amat)  
  w <- which(amat != 0, arr.ind = TRUE)
  R <- relations::relation( graph = data.frame( A = w[, 1], B = w[, 2] ) )
  
  rowsInv <- unlist(R$.Data$domain$A)
  rowsToAppend <- setdiff(1:n, rowsInv)
  colsInv <- unlist(R$.Data$domain$B)
  colsToAppend <- setdiff(1:n, colsInv)
  
  ta <- as.matrix( relations::relation_incidence(R) )  
  
  for ( i in rowsToAppend ) {
    if ( i < max(rowsInv) ) {
      ta <- rbind( ta[1:i-1, ], rep( 0, ncol(ta) ), ta[i:nrow(ta), ] )
    } else {
      ta <- rbind( ta[1:nrow(ta),], rep( 0, ncol(ta) ) )
    }
  }
  
  for(j in colsToAppend) {
    if ( j < max(colsInv) ) {
      ta <- cbind( ta[, 1:j-1], rep( 0, nrow(ta) ), ta[, j:ncol(ta)] )
    } else {
      ta <- cbind( ta[, 1:ncol(ta)], rep( 0, nrow(ta) ) )
    }
  }
  
  r <- relations::relation( incidence = ta )
  
  closure <- relations::transitive_closure(r)
  relations::relation_incidence(closure)
  
  # the closure needs to be sorted
  rel <- as.matrix( relations::relation_incidence(closure) )
  sorted <- sort( as.numeric(rownames(rel) ), index.return = TRUE )
  rel <- rel[sorted$ix, ]
  rel <- rel[, sorted$ix]
  
  rel
}