initConf <- function(init, diss, n, p, inddiff = FALSE) {
  
  if (inddiff) diss <- as.dist(apply(simplify2array(lapply(diss, as.matrix)), c(1,2), sum, na.rm = TRUE))
  
  if (length(init) == 1) {
   if (init == "torgerson") {
     meandiss <- mean(diss, na.rm = TRUE)   ## mean dissimilarity for imputation
     diss1 <- as.matrix(diss)
     diss1[is.na(diss1)] <- meandiss
     x <- torgerson(diss1, p = p) 
     init <- "dummy"
   }
   if (init == "random") {
     x <- matrix(runif(n*p, min = -1), ncol = p)
   }
  }
  if (is.data.frame(init)) init <- as.matrix(init)
  if (is.matrix(init)) {
    x <- as.matrix(init)
    if (any(dim(x) != c(n,p))) stop(paste0("Dimension of the starting configuration matrix needs to be ", n, " times ", p, "!"))
  }
  return(x)
}