cldeviance <- function(object, groups.gr = "rawscore", pi.hat)
{
# computes the collapsed deviance of
# object of class ppar


  k <- dim(object$X)[2]                             #number of items
  N <- dim(object$X)[1]                             #number of persons (full)

  #----------- define group vector ---------------
  if (groups.gr == "rawscore") indvec.full <- rowSums(object$X, na.rm = TRUE)    #person raw scores
  #if (groups.gr == "pattern") {                                                  #pattern-wise
  #  X.string <- apply(object$X, 1, paste, collapse = "")
  #  indvec.full <- rank(X.string, ties.method = "min")
  #} 
  if (is.numeric(groups.gr)) {
    if (length(groups.gr) != dim(object$X)[1]) stop("Group vector must be of length N (number of subjects in object$X)!")
    indvec.full <- groups.gr
  }
  #---------- end define group vector -----------
  
  #---- reduce group vector (pers.ex)------
  if (length(object$pers.ex) > 0) {                 #persons eliminated
    y <- object$X[-object$pers.ex,]                #observed values
    indvec.red <- indvec.full[-object$pers.ex]
  } else {
    y <- (object$X)
    indvec.red <- indvec.full
  }

  #pi.hat <- pmat(object)
  #gmemb.ext <- rep(object$gmemb, each = k)          #gmemb extended to response vector
  #pi.hat <- as.vector(t(pmat(object)))              #fitted values
  
  dev.g <- tapply(1:length(indvec.red), indvec.red, function(ii) {     #D component for each group
                  n.g <- length(ii)                                    #number of group subjects
                  y.g <- colSums(rbind(y[ii,]))                        #group responses
                  pi.g <- rbind(pi.hat[ii,])[1,]                   #vector with fitted values
                  devvec <- mapply(function(yy, pp) {                  #compute deviance for each item
                             if ((yy > 0) && (yy < n.g)) {
                               term1 <- yy*log(yy/(n.g*pp))
                               term2 <- (n.g-yy)*log((n.g-yy)/(n.g*(1-pp)))
                               dev <- sign(yy-n.g*pp)*sqrt(2*(term1+term2))
                             }
                             if (yy == 0) dev <- -sqrt(2*n.g*abs(log(1-pp)))
                             if (yy == n.g) dev <- sqrt(2*n.g*abs(log(pp)))
                             return(dev)
                            },y.g, pi.g)
                  return(sum(devvec^2))                                #item-wise sum of squared devres
                })

  value <- sum(dev.g)
  df <- (length(unique(indvec.red)))*k
  p.value <- 1-pchisq(value, df = df)

  result <- list(value = value, df = df, p.value = p.value)
  return(result)
}
