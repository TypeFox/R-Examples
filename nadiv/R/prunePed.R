################################################
#Borrowed heavily from the 'prunePed' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################
prunePed <- function (pedigree, phenotyped) {
   nPed <- numPed(pedigree[, 1:3])
   ikeep <- keep <- match(phenotyped, pedigree[, 1])
   nind <- length(ikeep) + 1
   while(length(ikeep) != nind){
      nind <- length(ikeep)
      ikeep <- union(c(nPed[ikeep, 2:3]), ikeep)
      ikeep <- ikeep[which(ikeep > 0)]
   }

   pedigree <- pedigree[sort(ikeep), ]
   pedigree[, 1] <- as.factor(as.character(pedigree[, 1]))
   pedigree[, 2] <- as.factor(as.character(pedigree[, 2]))
   pedigree[, 3] <- as.factor(as.character(pedigree[, 3]))
 pedigree
}

