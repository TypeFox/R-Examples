splitf <- function(set,fac,k) {
  prop <- table(fac)/length(fac)
  wl <- list()
  for (i in 1:nlevels(fac)) {wl[[i]] <- which(as.numeric(fac)==i)}
  res <- list()
  for (i in 1:(k-1)) {
    per.k <- length(unlist(wl))/(k-i+1)
    to.sample <- round(prop*per.k)
    ind <- NULL
    for (j in 1:length(wl)) {
	to.take <- sample(length(wl[[j]]),to.sample[j])
	ind <- c(ind,wl[[j]][to.take])
	wl[[j]] <- wl[[j]][-to.take]
    }
    res[[i]] <- set[ind,]
  }
  res[[k]] <- set[unlist(wl),]
  return(res)
}

