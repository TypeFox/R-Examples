parse.constraints <- function(constr, k=2, allsame=FALSE) {
  if (!is.null(constr)) {
    alpha <- category <- rep(NA, k)
    u <- sum(is.na(constr)) # Number of unconstrained means
    category[is.na(constr)] <- 1:u
    numcon <- as.numeric(constr) # All non-numbers become NA here.
    category[!is.na(numcon)] <- 0  # These means will remain fixed constants
    alpha[!is.na(numcon)] <- numcon[!is.na(numcon)]  
    w <- which(is.na(numcon) & !is.na(constr)) # w indexes only character constraints
    if (length(w)>0) {  # Process all the character (equality) constraints
      constr <- as.character(constr)
      alpha[w] <- constr[w]
      nc <- nchar(constr)[w] # nc:  number of characters (in each constraint)
      ss <- substr(constr[w], nc, nc) # extract final character
      category[w] <- u + match(ss, unique(ss)) # make categories for final characters
      alpha[w] <- substr(alpha[w], 1, nc-1) # Now extract number before final character
      alpha[w][alpha[w]=="-"] <- "-1"
      alpha[w][alpha[w]==""] <- "1"
    }
    alpha <- as.numeric(alpha) # Entire alpha vector should be numeric or NA now.
  } else {
    alpha <- rep(1,k)
    if(allsame) 
      category <- alpha
    else 
      category <- 1:k
  }
  list(category=category, alpha=alpha)
}

