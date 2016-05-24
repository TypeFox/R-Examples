epi.cp <- function(dat){
   obs <- as.data.frame(cbind(id = 1:nrow(dat), cp = rep(0, times = nrow(dat))))
   nvar <- dim(dat)[2]
   cp <- unique(dat)
   cp <- cbind(id = 1:nrow(cp),cp)

   for(i in 1:nrow(cp)){
      tmp <- rep(0, times = nrow(dat))
      for(j in 1:nvar){
        tmp. <- as.numeric(dat[,j] == cp[i,(j+1)])
        tmp <- cbind(tmp, tmp.)
        }
      tmp <- apply(tmp, MARGIN = 1, FUN = sum)
      id <- tmp == nvar
      obs$cp[id] <- cp$id[i]
   }
   n <- hist(obs$cp, breaks = seq(from = 0, to = nrow(cp), by = 1), plot = FALSE)
   n <- n$counts
   end <- nvar + 1
   cov.pattern <- as.data.frame(cbind(id = cp[,1], n, cp[,2:end]))
   rval <- list(cov.pattern = cov.pattern, id = obs$cp)
   rval
}
