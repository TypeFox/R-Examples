# merge fusion1 and pheno keeping only id's that are in both
fusion1m <- merge(fusion1, pheno, by='id', all.x=FALSE, all.y=FALSE)
