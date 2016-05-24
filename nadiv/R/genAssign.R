genAssign <- function(pedigree)
{ 
   n <- dim(pedigree)[1]
   generation <- vector("integer", length = n)
   if(any(apply(pedigree[, 1:3], MARGIN = 2, FUN = function(x){min(x, na.rm = TRUE) < 0}))){
      warning("Negative values in pedigree interpreted as missing values")
      pedigree[pedigree < 0] <- -998
   }
   if(!all(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.numeric)) | any(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.na))){
      pedigree[, 1:3] <- numPed(pedigree[, 1:3])
   }

   Cout <- .C("ga",
	as.integer(pedigree[, 2] - 1),
	as.integer(pedigree[, 3] - 1),
        generation,
	as.integer(n))

  Cout[[3]]
}
