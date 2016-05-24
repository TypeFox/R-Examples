Fstat <-
function (genotypes, npop, pop.mbrship, ploidy = 2) 
{
    if (ploidy == 1) 
        stop("Fstat not implemented for haploid data")
    format <- FormatGenotypes(genotypes, ploidy = ploidy)
    genotypes <- format$genotypes
    allele.numbers <- format$allele.numbers
    if (sum(is.na(genotypes)) != 0) {
        warning("Genotypes contain missing values which might bias computations")
        genotypes[is.na(genotypes)] <- -999
    }
    Fis = rep(-999, npop)
    Fst = matrix(nrow = npop, ncol = npop, -999)
    if (npop > 1) {
        for (iclass1 in 1:(npop - 1)) {
            for (iclass2 in (iclass1 + 1):npop) {
                sub1 <- pop.mbrship == iclass1
                sub2 <- pop.mbrship == iclass2
                if ((sum(sub1) != 0) & (sum(sub2) != 0)) {
                  ztmp <- genotypes[sub1 | sub2, ]
                  nindivtmp <- nrow(ztmp)
                  pop.mbrshiptmp <- pop.mbrship[sub1 | sub2]
                  pop.mbrshiptmp[pop.mbrshiptmp == iclass1] <- 1
                  pop.mbrshiptmp[pop.mbrshiptmp == iclass2] <- 2
                  tabindiv <- matrix(nrow = nindivtmp, ncol = 2, 
                    data = -999)
                  kk <- numeric(2)
                  effcl <- table(pop.mbrshiptmp)
                  nloc <- length(allele.numbers)
                  nloc2 <- 2 * nloc
                  Fistmp <- Fsttmp <- Fittmp <- -999
                  res <- .Fortran("ggfst", PACKAGE = "Geneland", 
                    as.integer(nindivtmp), as.integer(nloc), 
                    as.integer(nloc2), as.integer(allele.numbers), 
                    as.integer(2), as.integer(effcl), as.integer(ztmp), 
                    as.integer(pop.mbrshiptmp), as.integer(tabindiv), 
                    as.integer(kk), as.double(Fistmp), as.double(Fsttmp), 
                    as.double(Fittmp))
                  Fst[iclass1, iclass2] <- res[[12]][1]
                }
            }
        }
    }
    for (iclass1 in 1:npop) {
        sub1 <- pop.mbrship == iclass1
        if ((sum(sub1) != 0)) {
            ztmp <- rbind(genotypes[sub1, ], genotypes[sub1, 
                ])
            nindivtmp <- nrow(ztmp)
            pop.mbrshiptmp <- c(rep(1, sum(sub1)), rep(2, sum(sub1)))
            tabindiv <- matrix(nrow = nindivtmp, ncol = 2, data = -999)
            kk <- numeric(2)
            effcl <- table(pop.mbrshiptmp)
            nloc <- length(allele.numbers)
            nloc2 <- 2 * nloc
            Fistmp <- Fsttmp <- Fittmp <- -999
            res <- .Fortran("ggfst", PACKAGE = "Geneland", as.integer(nindivtmp), 
                as.integer(nloc), as.integer(nloc2), as.integer(allele.numbers), 
                as.integer(2), as.integer(effcl), as.integer(ztmp), 
                as.integer(pop.mbrshiptmp), as.integer(tabindiv), 
                as.integer(kk), as.double(Fistmp), as.double(Fsttmp), 
                as.double(Fittmp))
            Fis[iclass1] <- res[[11]][1]
        }
    }
    Fst[lower.tri(Fst, diag = TRUE)] <- 0
    Fst <- Fst + t(Fst)
    Fis[Fis == -999] <- NA
    Fst[Fst == -999] <- NA
    list(Fis = Fis, Fst = Fst)
}
