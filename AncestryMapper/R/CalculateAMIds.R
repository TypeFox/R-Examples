#' calculateAMids
#'
#' \code{calculateAMids} returns for each individual a set of distances to a given group of population references
#'
#' Calculates genetics distance between a sample population and population references.
#' 
#' @param pedtxtFile Name of the ped file to be used. The ped file should include all 51 HGDP references and the individuals for which the user wishes to calculate the genetic distance.
#' @param fileReferences File detailing the individuals in the ped file that correspond to the references, and the populations they refer to. A file that uses the 51 HGDP reference populations is provided with the package.
#' @return returns a dataframe containing the genetic distance of each sample individual to the population references. This dataframe can be subsequently used as input to \code{plotAMids}.
#' @examples
#'
#' \dontrun{
#' library(AncestryMapper)
#' HGDP.References <- system.file('extdata', 'HGDP.References.txt', package='AncestryMapper')
#' HGDP.500SNPs <- system.file('extdata', 'HGDP.500SNPs.ped', package='AncestryMapper')
#' HGDP.Phenotypes <- system.file('extdata', 'HGDP.Phenotypes.txt', package='AncestryMapper')
#' genetic.distance <- calculateAMids(pedtxtFile=HGDP.500SNPs, fileReferences=HGDP.References)
#' }

calculateAMids <- function(pedtxtFile, fileReferences){

    FUN.sumPedG <- function(row){
        row <- as.numeric(row)
        ## Treat 0s and NAs
        row.na <- which(is.na(row))
        row.0 <- which(row==0)
        row.change <- c(row.na, row.0)
        
        ## When row is NA or O - random substitution of 1 or 2
        row[row.change] <- sample(c(1,2), length(row.change), replace=T) 
        line.i <- matrix(unlist(row), length(row)/2, 2, byrow=T)
        line.i <- line.i[,1] + line.i[,2]
        return(line.i)
    }
    
    ## function to calculate distances; euclidean divided by number of snps
    ## maximum of distance 2, minimum 0
    FUN.dist2Ref <- function(row){
        row <- as.numeric(row)
        distOut <- sqrt(as.vector(rowMeans(referencesMatrix^2) + mean(row^2) - 2 * referencesMatrix %*% row/ncol(referencesMatrix)))
        distOut <- round(distOut, 4)
        return(distOut)
    }

    ## Check only 1,2 and NAs
    FUN.checkVectors <- function(row){
        pedUnique <- unique(row)
        if(!identical(setdiff(pedUnique, c(0,1,2,NA)), character(0)))
            stop("Input file can have 0, 1, 2 and NA values only.")
        return(NULL)
    }
    
    
    ## Calculate sumPed from scratch, ped, map
    ## Import ped, faster with readLines so that it is faster
    ped <- readLines(pedtxtFile)
    
    ## strplit, each vector by ' '
    ped <- lapply(ped, function(row) unlist(strsplit(row, split=' ')))

    ## Individual names
    indNames <- unlist(lapply(ped, function(row) row[2])) ### row 2 because ped file is Id

    ##Check duplicates; only looks at the individual ids not fam ids
    if(sum(duplicated(indNames))!=0) stop("There are duplicate individuals in the input file.")

    ## Extract genotypes, remove pedigree informations
    ped <- lapply(ped, function(row) row[7:length(row)])
    
    ## Check each vector has only 0, 1, 2 and NAs
    checkRows <- lapply(ped, FUN.checkVectors)

    ##Sum each genotype per snp; add snps; convert 2 cols into one
    sumPed <- lapply(ped, FUN.sumPedG)
    print('Summation of genotypes complete.')
    sumPed <- as.data.frame(do.call('rbind', sumPed))
    row.names(sumPed) <- indNames
    
    ## Calculate indices: Assign 100 to most similar reference, 0 to least similar
    FUN.indexes <- function(vec) round(100 - ((vec - min(vec))/(max(vec) - min(vec)) * 100), 1)
    
    ## Load the references file
    references <- read.table(fileReferences, header=T, as.is=T)

    ## References matrix: 51 references
    ## Error if references not in ped file format
    referencesMatrix <- as.matrix(sumPed[references[,'Reference'],])
    if(nrow(referencesMatrix)!=51)
        warning(paste('Number of references is usually 51; you are using ',nrow(referencesMatrix)), immediate.=T)

    ## Distances of every individual to references
    count <- 0
    dist2Ref <- apply(sumPed, 1, FUN.dist2Ref)
    print('Distances of Samples to References Computed.')
    dist2Ref <- as.data.frame(t(dist2Ref))

    ## Names; same as references
    names(dist2Ref) <- paste('C_', references$Pop, sep='')

    ## Indices; min:0 max:100
    Indexes <- data.frame(t(apply(dist2Ref, 1, FUN.indexes)))
    names(Indexes) <- gsub('C_', 'I_', names(Indexes))
    distAll <- data.frame(dist2Ref, Indexes)
    distAll <- data.frame(Id=row.names(distAll), distAll) 

    ## Return Distances
    return(distAll)
}
