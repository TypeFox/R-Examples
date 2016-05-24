"dccm.pca" <-
  function(x, pc = NULL, ncore = NULL, ...) {
   if (missing(x) || !"pca" %in% class(x))
     stop("dccm.pca: must supply a 'pca' object, i.e. from 'pca.xyz'")

   modes = pc

   ## Check for multiple cores
   ncore = setup.ncore(ncore)

   if(ncore > 1) {
     mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
     mccollect <- get("mccollect", envir = getNamespace("parallel"))
   }
   
   ## Set modes to be included
   if(is.null(modes))
      modes <- 1:length(x$L)
  
   ## If modes are negative, take modes complementary to them
   if( any(!is.numeric(modes)) || 
       any(!(abs(modes) %in% c(1:length(x$L)))) ||
       !(all(modes>0) || all(modes<0)) )
      stop("Incorrect mode index")
   if(all(modes < 0)) {
      modes <- setdiff(c(1:length(x$L)), abs(modes))
      if(length(modes) == 0)
         stop("No mode is selected")
   } 

   modes <- unique(modes)
   nmodes <- length(modes)
    
   ## Calc variance-covariance matrix over a subset of modes
   vcovmat <- function(r.inds, pca, vcov.mat = 0) {
     for ( i in seq_along(r.inds) ) {
       vcov.mat <- vcov.mat + (pca$U[, r.inds[i]] %o% pca$U[, r.inds[i]]) * pca$L[r.inds[i]]
       if(ncore > 1) writeBin(1, fpb)
       else setTxtProgressBar(pb, i)
     }
     return(vcov.mat)
   }

   ## Calculate variance-covariance matrix first ## 
   ## If contain $z, straightforward
   if(!is.null(x$z)) {

      q = x$z[, modes] %*% t(x$U[, modes])
      vcov = cov(q)

   } else {

      ## Initialize progress bar
      pb <- txtProgressBar(min=1, max=nmodes, style=3)

      if(ncore > 1) {   # Parallel

         # For progress bar
         fpb <- fifo(tempfile(), open = "w+b", blocking = T)

         # spawn a child process for message printing
         child <- mcparallel({ 
            progress <- 0.0
            while(progress < nmodes && !isIncomplete(fpb)) {
               msg <- readBin(fpb, "double")
               progress <- progress + as.numeric(msg)
               setTxtProgressBar(pb, progress)
            }
         } )
         ###################
 
         jobid <- rep(1:ncore, ceiling(nmodes/ncore))
         jobid <- jobid[1:nmodes]

         ltv <- mclapply(1:ncore, function(i) {
            j <- which(jobid %in% i)
            if(length(j) > 0) {
               m <- vcovmat(modes[j], x)
               m <- m[lower.tri(m, diag = TRUE)]
            } else {
               m = 0
            }
            return(m)
         } )
         
         ltv <- colSums(do.call(rbind, ltv))
         vcov <- matrix(0, nrow(x$U), nrow(x$U))
         vcov[lower.tri(vcov, diag = TRUE)] <- ltv
         vcov <- vcov + t(vcov)
         diag(vcov) <- diag(vcov) / 2

         close(fpb)
         mccollect(child) # End the child for message printing

      } else {       # Serial

         vcov <- vcovmat(modes, x)

      }
      close(pb)
   }
   
   corr.mat <- cov2dccm(vcov, ncore = ncore, ...)
   return(corr.mat)
}
