Qhat <- function(y, 
                 a, 
                 g, 
                 wgt = NULL){  

  pai <- table(a) / length(a)

  ind <- a == g

  if( is(wgt, "NULL") ) wgt <- numeric(length(ind)) + 1.0

  multiplier <- wgt * ind / pai[as.character(a)]

  Q <- sum( y * multiplier ) / sum( multiplier )

  return(Q)

}
