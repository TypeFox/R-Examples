## calculate the respective association measure between two variables of arbitrary types
association <- function(x, y){
  
  types <- c(data.class(x), data.class(y))
  
  # get binary variables (can be any of numeric, factor, ordered, logic)
  types[1] <- ifelse(length(na.omit(unique(x))) == 2, "binary", types[1])
  types[2] <- ifelse(length(na.omit(unique(y))) == 2, "binary", types[2])
  
  ## if (at least) one variable is continuous
  if(any(types == "numeric")){
    first <- which(types == "numeric")[1]
    second <- types[-first]
    if(first == 2){
      tmp <- x
      x <- y
      y <- tmp
    }
    
    # continuous - continuous/ordinal
    if(second == "numeric" | second == "ordered")
      res <- abs(cor(x, as.numeric(y), method="spearman", use="complete.obs"))
    
    # continuous - categorical
    else if(second == "factor")
      res <- assoc.rank.cat(x, y)

    # continuous - binary
    else if(second == "binary")
        res <- abs(myGKgamma(x, y))
  }
  
  ## otherwise, if (at least) one variable is ordinal
  else if(any(types == "ordered")){
    first <- which(types == "ordered")[1]
    second <- types[-first]
    if(first == 2){
      tmp <- x
      x <- y
      y <- tmp
    }
    
    # ordinal - ordinal
    if(second == "ordered")
      res <- abs(DescTools::GoodmanKruskalGamma(x, y))

    # ordinal - categorical
    else if(second == "factor")
      res <- assoc.rank.cat(x, y)
      
    # ordinal - binary
    else if(second == "binary")
      res <- abs(DescTools::GoodmanKruskalGamma(x, y))
  }
  
  ## if both variables are categorical/binary
  else
    res <- assoc.cat.cat(x, y)

  return(res)
}

