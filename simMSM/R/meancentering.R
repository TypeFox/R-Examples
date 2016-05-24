meancentering <- function(d, x, q, x.name = NULL, q.name = NULL){
  if(length(q) > 1.5){
    hi <- which(d$trans %in% q)
    x.q.mean <- mean(d[hi, x])
    ho <- d[, x] - x.q.mean
    ho[which(!(d$trans %in% q))] <- 0
    if(is.null(q.name)){
      q.name <- paste(q, collapse = ".")
    }
    if(is.null(x.name)){
      x.name <- x
    }
    name <- paste(x.name, q.name, sep = ".")
  }else{
    hi <- which(d$trans == q)
    x.q.mean <- mean(d[hi, x])
    ho <- d[, x] - x.q.mean
    ho[which(!(d$trans == q))] <- 0
    if(is.null(q.name)){
      q.name <- q
    }
    if(is.null(x.name)){
      x.name <- x
    }
    name <- paste(x.name, q.name, sep = ".")
  }
  return(list(x.q = ho, center = x.q.mean, name = name))
}