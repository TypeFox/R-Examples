FamiliasPosterior <- function (pedigrees, loci, datamatrix, prior, ref = 1, kinship = 0, 
    simplifyMutations = FALSE) 
{
    if (missing(pedigrees) || length(pedigrees) < 1) 
        stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.")
    if (class(pedigrees) == "pedigree" | class(pedigrees) == 
        "FamiliasPedigree") 
        pedigrees <- list(pedigrees)
    if (class(pedigrees) != "list") 
        stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.")
    for (i in pedigrees) {
        if (class(i) != "pedigree" && class(i) != "FamiliasPedigree") 
            stop("The pedigrees parameter must be an object of type 'pedigree' or 'FamiliasPedigree', or a list of such.")
    }
    npeds <- length(pedigrees)
    if (ref < 1 | ref > npeds) 
        stop("Impossible reference pedigree index.")
    if (missing(prior)) 
        prior <- rep(1/npeds, npeds)
    if (length(prior) != npeds) 
        stop("The prior argument must be a vector of the same length as the number of pedigrees, if it is not missing.")
    if (any(prior < 0) || round(sum(prior), 6) != 1) 
        stop("The prior must consist of non-negative numbers summing to 1.")
    if (missing(datamatrix)) 
        stop("The datamatrix must be supplied.")
    persons <- rownames(datamatrix)
    if (length(persons) < 1) 
        stop("The row names of the datamatrix must contain the names of the persons you have data for.")
    npers <- length(persons)
    indextable <- matrix(0, npers, npeds)
    firstped <- pedigrees[[1]]
    for (i in 1:npeds) for (j in 1:npers) {
        indextable[j, i] <- match(persons[j], pedigrees[[i]]$id, 
            nomatch = 0)
        if (indextable[j, i] == 0) 
            stop(paste("Error: Person ", persons[j], "of the data matrix does not occur in pedigree", 
                i))
        if (i > 1) {
            if (pedigrees[[i]]$sex[indextable[j, i]] != firstped$sex[indextable[j, 
                1]]) 
                stop("Persons common to all pedigrees must have the same sex in all pedigrees!")
        }
    }
    .C("NewFamilias")
    for (j in 1:npers) {
        result <- .C("AddPerson", as.integer(!(firstped$sex[indextable[j, 
            1]] == "female")), as.integer(-1), as.integer(FALSE), 
            index = integer(1), error = integer(1))
        if (result$error > 0) 
            stop("ERROR: Problems with common persons in pedigree.")
    }
    for (i in pedigrees) {
        nPersons <- length(i$sex)
        neworder <- rep(0, nPersons)
        nExMales <- nExFemales <- 0
        for (j in 1:nPersons) {
            mm <- match(i$id[j], persons, nomatch = 0)
            if (mm > 0) 
                neworder[j] <- mm
            else if (i$sex[j] == "female") {
                nExFemales <- nExFemales + 1
                neworder[j] <- nExFemales
            }
            else {
                nExMales <- nExMales + 1
                neworder[j] <- nExMales
            }
        }
        for (j in 1:nPersons) {
            if (!(i$id[j] %in% persons)) {
                if (i$sex[j] == "female") 
                  neworder[j] <- neworder[j] + npers
                else neworder[j] <- neworder[j] + npers + nExFemales
            }
        }
        result <- .C("AddPedigree", as.integer(nExFemales), as.integer(nExMales), 
            index = integer(1), error = integer(1))
        if (result$error > 0) 
            stop("ERROR: Wrong input in pedigrees.")
        index <- result$index + 1
        for (j in 1:nPersons) {
            if (i$findex[j] > 0) {
                result <- .C("AddRelation", as.integer(neworder[i$findex[j]] - 
                  1), as.integer(neworder[j] - 1), as.integer(index - 
                  1), error = integer(1))
                if (result$error == 1) 
                  stop("ERROR: Wrong input.")
                else if (result$error == 2) 
                  stop("ERROR: Illegal relation based on Year-of-birth or is-Child data.")
                else if (result$error == 3) 
                  stop("ERROR: Cycle in the pedigree or duplicate parent.")
            }
            if (i$mindex[j] > 0) {
                result <- .C("AddRelation", as.integer(neworder[i$mindex[j]] - 
                  1), as.integer(neworder[j] - 1), as.integer(index - 
                  1), error = integer(1))
                if (result$error == 1) 
                  stop("ERROR: Wrong input.")
                else if (result$error == 2) 
                  stop("ERROR: Illegal relation based on Year-of-birth or is-Child data.")
                else if (result$error == 3) 
                  stop("ERROR: Cycle in the pedigree or duplicate parent.")
            }
        }
    }
    if (missing(loci) || length(loci) < 1) 
        stop("The loci argument must be a FamiliasLocus object or a list of such.")
    if (class(loci) == "FamiliasLocus") 
        loci <- list(loci)
    if (class(loci) != "list") 
        stop("The loci argument must be a FamiliasLocus object or a list of such.")
    nloci <- length(loci)
    nms <- rep("", nloci)
    for (i in 1:nloci) {
        if (class(loci[[i]]) != "FamiliasLocus") 
            stop("The loci argument must be a FamiliasLocus object or a list of such.")
        nms[i] <- loci[[i]]$locusname
    }
    if (anyDuplicated(nms)) 
        stop("There can be no duplicated names of the loci.")
    for (i in loci) {
        if (any(i$alleles <= 0)) 
            stop(paste("ERROR: Problems with allele frequencies in locus", 
                i$locusname))
        if (round(sum(i$alleles), 6) != 1) 
            stop(paste("ERROR: Allele frequencies must sum to 1 in locus", 
                i$locusname))
        nAlleles <- length(i$alleles)
        if (!is.matrix(i$femaleMutationMatrix) | dim(i$femaleMutationMatrix)[1] != 
            nAlleles | dim(i$femaleMutationMatrix)[2] != nAlleles) 
            stop(paste("The female mutation matrix must be of a dimension corresponding to the vector of frequencies in locus", 
                i$locusname))
        if (any(as.vector(i$femaleMutationMatrix) < 0)) 
            stop(paste("The female mutation matrix cannot have negative entries in locus", 
                i$locusname))
        if (any(round(apply(i$femaleMutationMatrix, 1, sum), 
            6) != 1)) 
            stop(paste("The rows in the female mutation matrix must sum to 1 in locus", 
                i$locusname))
        if (!is.matrix(i$maleMutationMatrix) | dim(i$maleMutationMatrix)[1] != 
            nAlleles | dim(i$maleMutationMatrix)[2] != nAlleles) 
            stop(paste("The male mutation matrix must be of a dimension corresponding to the vector of frequencies in locus", 
                i$locusname))
        if (any(as.vector(i$maleMutationMatrix) < 0)) 
            stop(paste("The male mutation matrix cannot have negative entries in locus", 
                i$locusname))
        if (any(round(apply(i$maleMutationMatrix, 1, sum), 6) != 
            1)) 
            stop(paste("The rows in the male mutation matrix must sum to 1 in locus", 
                i$locusname))
        lOfArrays <- nAlleles * nAlleles
        simplifyMatrix <- i$simpleMutationMatrices | simplifyMutations
        hasSilentAllele <- (names(i$alleles)[nAlleles] == "silent" | 
            names(i$alleles)[nAlleles] == "Silent")
        result <- .C("AddAlleleSystem", as.integer(nAlleles), 
            as.integer(lOfArrays), as.double(i$femaleMutationMatrix), 
            as.double(i$maleMutationMatrix), as.integer(simplifyMatrix), 
            as.integer(nAlleles), as.double(i$alleles), as.integer(hasSilentAllele), 
            index = integer(1), error = integer(1))
        if (result$error > 0) 
            stop(paste("ERROR: Problems with input of allele system.", 
                i$locusname))
    }
    if (dim(datamatrix)[2] != 2 * nloci) 
        stop("The datamatrix must have two columns for each locus.")
    for (i in 1:nloci) {
        for (j in 1:npers) {
            A1 <- datamatrix[j, 2 * i - 1]
            A2 <- datamatrix[j, 2 * i]
            if (!(is.na(A1) && is.na(A2))) {
                if (is.na(A1)) 
                  A1 <- A2
                if (is.na(A2)) 
                  A2 <- A1
                M1 <- match(A1, names(loci[[i]]$alleles), nomatch = 0)
                if (M1 == 0) 
                  stop(paste("Allele", A1, "is not found in locus", 
                    loci[[i]]$locusname))
                M2 <- match(A2, names(loci[[i]]$alleles), nomatch = 0)
                if (M2 == 0) 
                  stop(paste("Allele", A2, "is not found in locus", 
                    loci[[i]]$locusname))
                result <- .C("AddDNAObservation", as.integer(j - 
                  1), as.integer(i - 1), as.integer(M1 - 1), 
                  as.integer(M2 - 1), error = integer(1))
                if (result$error > 0) 
                  stop("ERROR: Problems with input of marker data.")
            }
        }
    }
    if (kinship < 0) 
        stop("ERROR: Kinship cannot be negative.")
    result <- .C("GetProbabilities", as.double(1), as.integer(-1), 
        as.double(1), as.double(1), as.integer(TRUE), as.double(kinship), 
        redundant = integer(npeds), probabilities = double(npeds), 
        likelihoods = double(nloci * npeds), error = integer(1))
    if (result$error == 1) 
        stop("ERROR: Problems computing probabilities.")
    if (result$error == 2) 
        stop("ERROR: All pedigrees have probability zero.")
    pedigreeDuplicated <- as.logical(result$redundant)
    if (any(pedigreeDuplicated)) 
        stop(paste("ERROR: Some of the listed pedigrees were equivalent. Duplicated pedigrees are numbers", 
            (1:npeds)[pedigreeDuplicated]))
    #prior <- result$probabilities
    likelihoodsPerSystem <- matrix(result$likelihoods, nloci, 
        npeds)
    likelihoods <- apply(likelihoodsPerSystem, 2, prod)
    posterior <- prior * likelihoods
    posterior <- posterior/sum(posterior)
    LR <- likelihoods/likelihoods[ref]
    LRperMarker <- likelihoodsPerSystem/likelihoodsPerSystem[, 
        ref]
    .C("NewFamilias")
    names(posterior) <- names(pedigrees)
    names(prior) <- names(pedigrees)
    names(LR) <- names(pedigrees)
    colnames(LRperMarker) <- names(pedigrees) 
    locusnames <- rep("", nloci)   
    for (i in 1:nloci) locusnames[i] <- loci[[i]]$locusname
    rownames(LRperMarker) <- locusnames
    names(likelihoods) <- names(pedigrees)
    colnames(likelihoodsPerSystem) <- names(pedigrees)
    rownames(likelihoodsPerSystem) <- locusnames
    list(posterior = posterior, prior = prior, LR = LR, LRperMarker = LRperMarker, 
        likelihoods = likelihoods, likelihoodsPerSystem = likelihoodsPerSystem)
}