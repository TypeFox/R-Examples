 balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
         # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
         # nice we to get the ids in a list, split by class
###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
cat(i)
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
### Now do a clever sort to mix up the NAs
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }
