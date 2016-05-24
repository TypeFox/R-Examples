founderLine <- function(pedigree, sex){
   colsel <- match(sex, names(pedigree))
   if(!colsel %in% seq(ncol(pedigree))){
      stop("character argument to 'sex' must exactly match a column name in 'pedigree'")
   }
   nPed <- numPed(pedigree[, 1:3])
   line <- par <- nPed[, colsel]
   parKnown <- par > 0
   while(any(parKnown)){
      par[parKnown] <- nPed[line[parKnown], colsel]
      parKnown <- par > 0
      line[parKnown] <- par[parKnown]
   }
   line[which(line < 0 & pedigree[, 4] == pedigree[line[line > 0][1], 4])] <- which(line < 0 & pedigree[, 4] == pedigree[line[line > 0][1], 4]) 
   line[line < 0] <- NA 
 if(is.factor(pedigree[, 1])) as.character(pedigree[line, 1]) else pedigree[line, 1]
}

