#########################################
#' Read genomic data to perform QTL analyses.
#'
#'
#' This function reads genomic data and is similar to the read.cross
#' function from r/qtl package (Broman and Sen, 2009) but allows
#' importing data from a flapjack format (Milne et al., 2010).
#' The files required include a file containing phenotypic
#' information (P.data), a file containing genotypic information (G.data),
#' and a file containing map information (map.data) for all markers.
#'
#' @usage qtl.cross (P.data = NULL, G.data, map.data, cross, heterozygotes = TRUE,
#'                   sep = "\t" )
#'
#' @param P.data Name of the file containing phenotypic information.
#' Each row represents the individuals while each column represents
#' the phenotypic traits. The first column should be labeled as 'genotype'
#' and should contain identification name for each individual.
#' The name of each trait should also be included.
#'
#' @param G.data Name of the file containing genotypic (marker scores)
#' information. Each row represents the individuals
#' while each column represents the markers. Headers for markers should be
#' included, but not for genotypes.
#' The first column contains the names of the genotypes.
#' The first row contains the names of the markers.
#' The marker genotypes are coded by two characters corresponding
#' to the alleles using a separator between alleles (by default a slash /).
#' If a single character is given, the genotype is assumed to be homozygous.
#' Missing values are indicated by default with '-'.
#' In the example below, the two alleles have been called 1 and 2 because it
#' is useful to link alleles to their origin, i.e. parent 1 or parent 2.
#' Therefore, 1 corresponds to homozygous for allele 1 (synonymous to 1/1),
#' 1/2 corresponds to heterozygous, and 2 corresponds to homozygous
#' for allele 2 (synonymous to 2/2).
#' In the case of partially informative markers (e.g. dominant markers)
#' genotypes are coded as 1/- or 2/-, depending on whether the dominant
#' allele originated from parent 1 or parent 2.
#'
#' @param  map.data Name of the file containing marker map information
#' (i.e. linkage group and position within linkage group).
#' The file is a text tab delimited file. Each row represents markers.
#' The file consists of three columns.
#' Column 1 gives the marker names,
#' column 2 the chromosome on which the marker has been mapped,
#' and column 3 indicates the position of the marker within the chromosome.
#'
#' @param cross The type of population studied. The type of population studied.
#' Options are: F2 (f2), doubled haploids (dh), backcross (bc),
#' recombinant inbred lines from selfing (riself, ri4self, or ri8self
#' depending on the number of parents used), recombinant inbred lines from
#' sib-mating (risib, ri4sib, or ri8sib depending on the number of parents
#' used), segregating F1 cross-pollinated populations (cp),
#'
#' @param heterozygotes It indicates whether there are heterozygotes or
#' not in the association mapping population. FALSE is set as default.
#'
#' @param sep To define the espace between the data.
#'
#' @return Creates an object of class cross to be used in QTL analysis.
#' The components are the same as r/qtl (Broman and Sen, 2009):
#' geno  This is a list with elements corresponding to chromosomes.names (geno)
#' contains the names of the chromosomes. Each chromosome is itself a list,
#' and is given class A or X according to whether it is autosomal or
#' the X chromosome.There are two components for each chromosome:
#' data, a matrix whose rows are individuals and whose columns are markers,
#' and map, either a vector of marker positions (in cM) or a matrix of dim
#' (2 x n.mar) where the rows correspond to marker positions in female and male
#' genetic distance, respectively.The genotype data gets converted into numeric
#' codes, as follows.
#' The genotype data for a backcross is coded as NA = missing, 1 = AA, 2 = AB.
#' For an F2 intercross, the coding is NA = missing, 1 = AA, 2 = AB, 3 = BB,
#' 4 = not BB (i.e. AA or AB; D in Mapmaker/qtl), 5 = not AA (i.e. AB or BB;
#' C in Mapmaker/qtl).
#' For a 4-way cross, the mother and father are assumed to have genotypes AB
#' and CD, respectively. The genotype data for the progeny is assumed to be
#' phase-known, with the following coding scheme: NA = missing, 1 = AC, 2 = BC,
#' 3 = AD, 4 = BD, 5 = A = AC or AD, 6 = B = BC or BD, 7 = C = AC or BC,
#' 8 = D = AD or BD, 9 = AC or BD, 10 = AD or BC, 11 = not AC, 12 = not BC,
#' 13 = not AD, 14 = not BD.
#' pheno   a data.frame of size (n.ind x n.phe) containing the phenotypes.
#' If a phenotype with the name genotype  is included, these identifiers will
#' be used in top.errorlod, plotErrorlod, and plotGeno as identifiers for the
#' individual.
#'
#' @references Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl.
#'             Springer, NewYork
#'             Comadran J, Thomas W, van Eeuwijk F, Ceccarelli S, Grando S, Stanca A,
#'             Pecchioni N, Akar T, Al-Yassin A, Benbelkacem A, Ouabbou H, Bort J,
#'             Romagosa I, Hackett C, Russell J (2009) Patterns of genetic diversity
#'             and linkage disequilibrium in a highly structured Hordeum vulgare
#'             association-mapping population for the Mediterranean basin.
#'             Theor Appl Genet 119:175-187
#'             Milne et al., (2010) Flapjack - graphical genotype visualization.
#'             Bioinformatics 26(24), 3133-3134.
#'
#' @author Lucia Gutierrez.
#'
#' @details The function creates an intermediate file called 'temp.csv' and
#' then uses the read.cross from r/qtl to read it.
#' The output object is an object of class=cross, the same as the
#' one produced by the function read.cross in r/qtl (Broman and Sen, 2009)
#'
#' @note All functions in this package uses cross data style.
#'
#' @seealso qtl.analysis, qtl.memq
#'
#' @import qtl
#' @import stringr
#' @import graphics
#' @import utils
#' @import grDevices
#' @import stats
#'
#' @export
#'
#' @examples
#' data (SxM_geno)
#' data (SxM_map)
#' data (SxM_pheno)
#'
#' P.data <- SxM_pheno
#' G.data <- SxM_geno
#' map.data <- SxM_map
#'
#' cross.data <- qtl.cross (P.data, G.data, map.data,
#' cross='dh', heterozygotes = FALSE)
#' summary (cross.data)
#'
qtl.cross <- function (P.data = NULL, G.data, map.data, cross = "dh",
  heterozygotes = TRUE, sep = "\t") {
    cross.data <- NULL

    G.data <- as.matrix(G.data)
    G.data <- G.data[3:nrow(G.data),]
    G.data[G.data == " 0"] <- "0"
    G.data[G.data == " 1"] <- "1"

    G.data <- G.data[which(G.data[, 1] %in% as.character(P.data$genotype)), ]
    P.data <- P.data[which(as.character(P.data$genotype) %in% G.data[, 1]), ]
    rownames(G.data) <- G.data[, 1]
    G.data <- G.data[, 2:ncol(G.data)]

    G.data <- G.data[, which(colnames(G.data) %in% map.data[, 1])]
    map.data <- map.data[which(map.data[, 1] %in% colnames(G.data)), ]
    G.data <- G.data[, match(map.data[, 1], colnames(G.data))]

    P.data <- P.data[match(rownames(G.data), as.character(P.data$genotype)), ]

    a <- matrix(rep("", (ncol(P.data) * 2)), 2, ncol(P.data))

    colnames(a) <- colnames(P.data)

    names(a) <- names(P.data)

    b <- data.frame(t(map.data[, 2:3]))

    names(b) <- as.character(map.data[, 1])

    c <- data.frame(a, b, check.names = FALSE, stringsAsFactors = FALSE)

    rownames(c) <- c(rownames(a), rownames(b))

    # extract id information for geno and pheno
    G.id <- row.names(G.data)
    P.id <- as.character(P.data$genotype)

    if (sum(G.id != P.id) > 0) {
        simpleError("IDs don't match")
    }

    geno.data <- G.data

    gens <- c(unlist(geno.data))
    gens <- sort(unique(gens[gens != "-"]))

    als <- unique(c(unlist(strsplit(gens, "/"))))

    Het <- gens[which(gens == paste(als[1], als[2],
      sep = "/") | gens == paste(als[2], als[1],
        sep = "/"))]
    if (length(Het) > 0) {
        hetpos <- grep(Het, geno.data)
        Apos <- setdiff(grep(als[1], geno.data), hetpos)
        Bpos <- setdiff(grep(als[2], geno.data), hetpos)
        geno.data[hetpos] <- "1/2"
    }
    if (length(Het) == 0) {
        Apos <- grep(als[1], geno.data)
        Bpos <- grep(als[2], geno.data)
    }
    geno.data[Apos] <- "1/1"
    geno.data[Bpos] <- "2/2"
    geno.data[geno.data == "1/1"] <- "AA"
    geno.data[geno.data == "2/2"] <- "BB"
    geno.data[geno.data == "1/2"] <- "AB"
    geno.data[geno.data == "1/-"] <- "C"
    geno.data[geno.data == "2/-"] <- "D"
    geno.data[geno.data == "-"] <- NA

    e <- cbind(P.data, geno.data)

    f <- rbind(c, e)

    names(f)[1] <- "id"

    zz <- file("cross.data.file.csv", "w")

    write.csv(f, zz, quote = FALSE, row.names = FALSE)

    close(zz)

    # this is to import the cross

    if(cross == "dh"){
      data.prueba <- read.cross("csv", file = "cross.data.file.csv",
        genotypes = c("AA", "BB"))
    } #check if more

    if(cross == "bc"){
      data.prueba <- read.cross("csv", file = "cross.data.file.csv",
        genotypes = c("AA", "AB"))
    } #check if more

    if(cross == "f2"){
      data.prueba <- read.cross("csv", file = "cross.data.file.csv",
        genotypes = c("AA", "AB", "BB", "C", "D"))
    }

    if(cross == "riself" | cross == "ri4self" | cross == "ri8self"){
      data.prueba <- read.cross("csv", file="cross.data.file.csv",
        genotypes=c("AA", "BB", "C", "D"))
    }

    if(cross == "risib" | cross == "ri4sib" | cross == "ri8sib"){
      data.prueba <- read.cross("csv", file = "cross.data.file.csv",
        genotypes = c("AA", "BB", "C", "D"))
    }

    data.prueba$gwas <- "QTL"

    cross.data <- data.prueba

    class(cross.data)[1] <- paste(cross)

    cross.data
}
