# Class definitions for mutation classes

#' @include shazam.R
NULL

#### Classes ####

#' S4 class defining replacement and silent mutation definitions
#' 
#' \code{MutationDefinition} defines a common data structure for defining the whether
#' a mutation is annotated as a replacement or silent mutation.
#' 
#' @slot    name            name of the MutationDefinition.
#' @slot    description     description of the model and its source.
#' @slot    classes         named character vectors with single-letter amino acid codes as names
#'                          and amino acid classes as values, with \code{NA} assigned to set of 
#'                          characters \code{c("X", "*", "-", ".")}. Replacement (R) is be 
#'                          defined as a change in amino acid class and silent (S) as no 
#'                          change in class. 
#' @slot    codonTable      matrix of codons (columns) and substitutions (rows).
#' @slot    citation        publication source.
#' 
#' @seealso
#' See \link{MUTATION_SCHEMES} for a set of predefined \code{MutationDefinition} objects.
#'    
#' @name         MutationDefinition-class
#' @rdname       MutationDefinition-class
#' @aliases      MutationDefinition
#' @exportClass  MutationDefinition
setClass("MutationDefinition", 
         slots=c(name="character",
                 description="character",
                 classes="character",
                 codonTable="matrix",
                 citation="character"))


#### Builder functions ####

# Create all codons one mutation away from input codon.
#
# All codons one mutation away from the input codon are generated.
#
# @param codon   starting codon to which mutations are added
# @return a vector of codons.
allCodonMuts <- function(codon) {
    codon_char <- seqinr::s2c(codon)
    matCodons <- t(array(codon_char, dim=c(3,12)))
    matCodons[1:4, 1] <- NUCLEOTIDES[1:4]
    matCodons[5:8, 2] <- NUCLEOTIDES[1:4]
    matCodons[9:12,3] <- NUCLEOTIDES[1:4]
    return(apply(matCodons, 1, seqinr::c2s))
}


# Generate codon table
#
# First generates all informative codons and determines types of mutations.
# Next generates uninformative codons (having either an N or a gap "-"
# character) and sets the mutation type as NA.
#
# @param   aminoAcidClasses  vector of amino acid trait classes
#                            if NULL then R or S is determined by amino acid identity
# @return  matrix with all codons as row and column names and the type of mutation as 
#           the corresponding value in the matrix.
computeCodonTable <- function(aminoAcidClasses=NULL) {
    # Initialize empty data.frame
    codon_table <- as.data.frame(matrix(NA, ncol=64, nrow=12))
    
    # Pre-compute every codon
    counter <- 1
    for(pOne in NUCLEOTIDES[1:4]) {
        for(pTwo in NUCLEOTIDES[1:4]) {
            for(pThree in NUCLEOTIDES[1:4]) {
                codon <- paste0(pOne, pTwo, pThree)
                colnames(codon_table)[counter] <- codon
                counter <- counter + 1
                all_muts <- allCodonMuts(codon)
                codon_table[, codon] <- sapply(all_muts, function(x) { mutationType(x, codon, aminoAcidClasses=aminoAcidClasses) })
            }
        }
    }
    
    # Set codons with N or . to be NA
    chars <- c("N","A","C","G","T", ".")
    for(n1 in chars) {
        for(n2 in chars) {
            for(n3 in chars) {
                if(n1=="N" | n2=="N" | n3=="N" | n1=="." | n2=="." | n3==".") {
                    codon_table[, paste0(n1, n2, n3)] <- rep(NA, 12)
                }
            }
        }
    }
    
    return(as.matrix(codon_table))
}


#' Creates a MutationDefinition
#' 
#' \code{createMutationDefinition} creates a \code{MutationDefinition}.
#'
#' @param    name           name of the mutation definition.
#' @param    classes        named character vectors with single-letter amino acid codes as names
#'                          and amino acid classes as values, with \code{NA} assigned to set of 
#'                          characters \code{c("X", "*", "-", ".")}. Replacement (R) is be 
#'                          defined as a change in amino acid class and silent (S) as no 
#'                          change in class. 
#' @param    description    description of the mutation definition and its source data.
#' @param    citation       publication source.
#' 
#' @return   A \code{MutationDefinition} object.
#' 
#' @seealso  See \code{\link{MutationDefinition}} for the return object.
#' 
#' @examples
#' # Define hydropathy classes
#' library(alakazam)
#' hydropathy <- list(hydrophobic=c("A", "I", "L", "M", "F", "W", "V"),
#'                    hydrophilic=c("R", "N", "D", "C", "Q", "E", "K"),
#'                    neutral=c("G", "H", "P", "S", "T", "Y"))
#' chars <- unlist(hydropathy, use.names=FALSE)
#' classes <- setNames(translateStrings(chars, hydropathy), chars)
#'
#' # Create hydropathy mutation definition
#' md <- createMutationDefinition("Hydropathy", classes)
#' 
#' @export
createMutationDefinition <- function(name,
                                     classes,
                                     description="",
                                     citation="") {
    # Build the codon table
    codonTable <- computeCodonTable(aminoAcidClasses=classes)

    # Define MutationDefinition object
    md <- new("MutationDefinition",
              name=name,
              description=description,
              classes=classes,
              codonTable=codonTable,
              citation=citation)
    
    return(md)
}


#### Data ####

#' Amino acid mutation definitions
#'
#' Definitions of replacement (R) and silent (S) mutations for different amino acid
#' physicochemical classes.
#'
#' @format A \code{\link{MutationDefinition}} object defining:
#' \itemize{
#'   \item  \code{CHARGE_MUTATIONS}:      Amino acid mutations are defined by changes
#'                                        in charge classes.
#'   \item  \code{HYDROPATHY_MUTATIONS}:  Amino acid mutations are defined by changes
#'                                        in hydrophobicitity classes.
#'   \item  \code{POLARITY_MUTATIONS}:    Amino acid mutations are defined by changes
#'                                        in polarity classes.
#' }
#' 
#' @references
#' \enumerate{
#'   \item  \url{http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html} 
#' }
#' 
#' @name   MUTATION_SCHEMES
NULL

#' @name    CHARGE_MUTATIONS
#' @rdname  MUTATION_SCHEMES
NULL

#' @name    HYDROPATHY_MUTATIONS
#' @rdname  MUTATION_SCHEMES
NULL

#' @name    POLARITY_MUTATIONS
#' @rdname  MUTATION_SCHEMES
NULL
