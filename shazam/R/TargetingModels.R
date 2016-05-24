# Targeting models

#' @include shazam.R
NULL

#### Data ####

#' Mouse single nucleotide distance matrix.
#'
#' Single nucleotide distance matrix of somatic hypermutation targeting based on 
#' Mus musculus Ig sequence data.
#'
#' @format   A symmetric matrix of nucleotide substitution distance scores with 
#'           row names and column names definition the specific subsitution.
#' 
#' @references
#' \enumerate{
#'   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#'            mutagenesis in normal and autoreactive B cells. 
#'            J Immunol. 1996 156:2642-52. 
#' }
#'
#' @seealso  See \code{\link{HS1FDistance}} for the human 1-mer distance matrix.
"M1NDistance"


#' Human single nucleotide distance matrix.
#'
#' Single nucleotide distance matrix of somatic hypermutation targeting based on 
#' human Ig sequence data.
#'
#' @format   A symmetric matrix of nucleotide substitution distance scores with 
#'           row names and column names definition the specific subsitution.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#' }
#'
#' @seealso  See \code{\link{HS5FModel}} for the 5-mer model from same publication.
"HS1FDistance"


#' Uniform 5-mer targeting model.
#'
#' A null 5-mer model of somatic hypermutation targeting where all substitution, mutability
#' and targeting rates are uniformly distributed.
#'
#' @format \code{\link{TargetingModel}} object.
#' 
#' @seealso  See \code{\link{HS5FModel}} the human 5-mer model.
"U5NModel"


#' Human 5-mer targeting model.
#'
#' 5-mer model of somatic hypermutation targeting based on analysis of silent mutations
#' in functional Ig sequences from Homo sapiens.
#'
#' @format \code{\link{TargetingModel}} object.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'  
#' @seealso  See \code{\link{HS1FDistance}} for the 1-mer distance matrix from the same 
#'           publication.
"HS5FModel"


#### Classes ####

#' S4 class defining a targeting model
#' 
#' \code{TargetingModel} defines a common data structure for mutability, substitution and
#' targeting of immunoglobulin (Ig) sequencing data in a 5-mer microsequence context.
#' 
#' @slot     name          Name of the model.
#' @slot     description   Description of the model and its source data.
#' @slot     species       Genus and species of the source sequencing data.
#' @slot     date          Date the model was built.
#' @slot     citation      Publication source.
#' @slot     substitution  Normalized rates of the center nucleotide of a given 5-mer 
#'                         mutating to a different nucleotide. The substitution model 
#'                         is stored as a 5x3125 matrix of rates. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         \code{c("A", "C", "G", "T", "N")}, and columns define the 
#'                         complete 5-mer of the unmutated nucleotide sequence.
#' @slot     mutability    Normalized rates of a given 5-mer being mutated. The 
#'                         mutability model is stored as a numeric vector of length 3125 
#'                         with mutability rates for each 5-mer.
#' @slot     targeting     Rate matrix of a given mutation ocurring, defined as 
#'                         \eqn{mutability * substitution}. The targeting model 
#'                         is stored as a 5x3125 matrix. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         \code{c("A", "C", "G", "T", "N")}, and columns define the complete 5-mer 
#'                         of the unmutated nucleotide sequence.
#' 
#' @seealso  See \link{createTargetingModel} building models from sequencing data.
#'           
#' @name         TargetingModel-class
#' @rdname       TargetingModel-class
#' @aliases      TargetingModel
#' @exportClass  TargetingModel
setClass("TargetingModel", 
         slots=c(name="character",
                 description="character",
                 species="character",
                 date="character",
                 citation="character",
                 mutability="numeric",
                 substitution="matrix",
                 targeting="matrix"),
         prototype=c(name="name",
                     description="description",
                     species="species",
                     date="2000-01-01",
                     citation="citation",
                     mutability=numeric(3125),
                     substitution=matrix(0, 5, 3125),
                     targeting=matrix(0, 5, 3125)))


#### Model building functions #####

#' Builds a substitution model
#'
#' \code{createSubstitutionMatrix} builds a 5-mer nucleotide substitution model by counting 
#' the number of substitution mutations occuring in the center position for all 5-mer 
#' motifs.
#'
#' @param    db                data.frame containing sequence data.
#' @param    model             type of model to create. The default model, "RS", creates 
#'                             a model by counting both replacement and silent mutations.
#'                             The "S" specification builds a model by counting only 
#'                             silent mutations.
#' @param    sequenceColumn    name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn    name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn       name of the column containing the V-segment allele call.
#' @param    multipleMutation  string specifying how to handle multiple mutations occuring 
#'                             within the same 5-mer. If \code{"independent"} then multiple 
#'                             mutations within the same 5-mer are counted indepedently. 
#'                             If \code{"ignore"} then 5-mers with multiple mutations are 
#'                             excluded from the total mutation tally.
#' @param    returnModel       string specifying what type of model to return; one of
#'                             \code{c("5mer", "1mer", "1mer_raw")}. If \code{"5mer"} 
#'                             (the default) then a 5-mer nucleotide context model is 
#'                             returned. If \code{"1mer"} or \code{"1mer_raw"} then a single 
#'                             nucleotide substitution matrix (no context) is returned;
#'                             where \code{"1mer_raw"} is the unnormalized version of the 
#'                             \code{"1mer"} model. Note, neither 1-mer model may be used
#'                             as input to \link{createMutabilityMatrix}.
#' @param    minNumMutations   minimum number of mutations required to compute the 5-mer 
#'                             substitution rates. If the number of mutations for a 5-mer
#'                             is below this threshold, its substitution rates will be 
#'                             estimated from neighboring 5-mers. Default is 50.                            
#' 
#' @return   A 4x1024 matrix of column normalized substitution rates for each 5-mer motif with 
#'           row names defining the center nucleotide, one of \code{c("A", "C", "G", "T")}, 
#'           and column names defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#' 
#' # Create model using only silent mutations and ignore multiple mutations
#' sub <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
#' 
#' @export
createSubstitutionMatrix <- function(db, model=c("RS", "S"), sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     vCallColumn="V_CALL",
                                     multipleMutation=c("independent", "ignore"),
                                     returnModel=c("5mer", "1mer", "1mer_raw"),
                                     minNumMutations=50)  {
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    returnModel <- match.arg(returnModel)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    db[[sequenceColumn]] = toupper(db[[sequenceColumn]])
    db[[germlineColumn]] = toupper(db[[germlineColumn]])
    
    # Setup
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_words <- seqinr::words(4, nuc_chars)
    
    # Define v_families (heavy or light chain) to only those found in the data
    v_families <- getFamily(db[[vCallColumn]])
    
    # Define empty return list of lists
    substitutionMatrix <- matrix(0, ncol=4, nrow=4, dimnames=list(nuc_chars, nuc_chars))
    substitutionList <- list()    
    for(v_fam in unique(v_families)) {
        substitutionList[[v_fam]] = list()
        for(word in nuc_words){
            substitutionList[[v_fam]][[word]] = substitutionMatrix
        }
    }
    
    # Remove IMGT gaps in the germline & input sequences
    matInputCollapsed <- removeCodonGaps(db[, c(sequenceColumn, germlineColumn)])
    # TODO: Unnecessary conversion
    db[[sequenceColumn]] <- matInputCollapsed[, 1]
    db[[germlineColumn]] <- matInputCollapsed[, 2]
    
    # Get mutations
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn,
                                       multipleMutation=multipleMutation,
                                       model=model)
    
    if (model == "S") { # Silent model
        for(index in 1:length(mutations)) {
            cSeq <-  s2c(db[index,sequenceColumn])
            cGL  <-  s2c(db[index,germlineColumn])
            indexMutation <- mutations[[index]]
            v_fam <- v_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL=="N") & !any(codonSeq=="N") ){
                    codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
                    codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
                    codonPermutate <- apply(codonPermutate,1,paste,collapse="")
                    codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
                    muType <- mutationTypeOptimized(codonPermutate)
                    if(!length(grep("N",wrd))){
                        if( sum(muType=="S") == length(muType) ){
                            substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] + 1)
                        }
                    }
                }
            }
        }
    } else if (model == "RS") { # RS model (All mutations)
        for (index in 1:length(mutations)) {
            cSeq <-  s2c(db[index,sequenceColumn])
            cGL  <-  s2c(db[index,germlineColumn])
            indexMutation <- mutations[[index]]
            v_fam <- v_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL=="N") & !any(codonSeq=="N") ){
                    if(!length(grep("N",wrd))){
                        substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] <- substitutionList[[v_fam]][[wrd]][glAtMutation,seqAtMutation] + 1
                    }
                }
            }
        }
    }
    
    
    # Convert substitutionList to listSubstitution to facilitate the aggregation of mutations
    arrNames <- c(outer(unique(v_families), nuc_words, paste, sep = "_"))
    listSubstitution <- array(0, dim=c(length(arrNames),4,4), dimnames=list(arrNames, nuc_chars, nuc_chars))
    
    for(v_fam in unique(v_families)){
        listSubstitution[paste(v_fam, nuc_words, sep="_"),,]<-t(sapply(nuc_words,function(word){substitutionList[[v_fam]][[word]]}))
    }
    
    # Aggregate mutations from all V families
    M<-list()
    subMat1mer <- matrix(0,4,4) # a single substitution matrix for all fivemers
    listSubNames<-sapply(dimnames(listSubstitution)[[1]],function(x)strsplit(x,"_",fixed=TRUE)[[1]])
    
    for (nuc_word in nuc_words) {
        M[[nuc_word]] <- t(sapply(1:4,function(i)apply(listSubstitution[listSubNames[2,] == nuc_word,i,],2,sum))) # sums mutations from all families
        rownames(M[[nuc_word]]) <- nuc_chars
        subMat1mer <- subMat1mer + M[[nuc_word]]
    }
    
    # Return 1-mer substitution model; this output cannot be used for createMutabilityMatrix
    if (returnModel == "1mer") {
        subMat1merNorm = t(apply(subMat1mer, 1, function(x){x/sum(x)}))
        return (subMat1merNorm)
    } else if (returnModel == "1mer_raw") {
        return (subMat1mer)
    } 
    
    # Aggregate mutations from neighboring bases for low frequency fivemers
    # fivemer=M; FIVEMER="CCATT"
    .simplifivemer <- function(fivemer, FIVEMER, Thresh=50) {
        Nuc=substr(FIVEMER,3,3)
        Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")
        
        # If the total number of mutation is greater than Thresh, and there are mutations to every other base
        if(sum(fivemer[[Nei]][Nuc,])>Thresh && sum(fivemer[[Nei]][Nuc,]==0)==1){
            return(fivemer[[Nei]][Nuc,]);
        }
        else{ # Otherwise aggregate mutations from 5-mers with the same inner 3-mer
            FIVE=fivemer[[Nei]][Nuc,]
            for(i in 1:4){
                for(j in 1:4){
                    MutatedNeighbor=paste(nuc_chars[i],substring(Nei,2,3),nuc_chars[j],collapse="",sep="")
                    FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                }
            }
            
            # If the total number of mutations is still not enough, aggregate mutations from all 5-mers 
            # i.e., use 1-mer model
            if(sum(FIVE) <= Thresh || sum(FIVE==0)!=1 ){
                FIVE=fivemer[[Nei]][Nuc,]
                MutatedNeighbors = seqinr::words(4, nuc_chars)
                for (MutatedNeighbor in MutatedNeighbors) {
                    FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                }
            } 
            
            return(FIVE)
        }
    }
    
    
    substitutionModel <- sapply(seqinr::words(5, nuc_chars), function(x) { .simplifivemer(M, x, Thresh=minNumMutations) })
    
    # Assign A->A, C->C, G->G, T->T to NA
    center_nuc <- gsub("..([ACGT])..", "\\1", colnames(substitutionModel))
    for (i in 1:length(center_nuc)) {
        substitutionModel[center_nuc[i], i] <- NA
    }
    
    # Normalize by column
    substitutionModel <- apply(substitutionModel, 2, function(x) { x / sum(x, na.rm=TRUE) })
    substitutionModel[!is.finite(substitutionModel)] <- NA
    
    return(substitutionModel)
}


#' Builds a mutability model
#'
#' \code{createMutabilityMatrix} builds a 5-mer nucleotide mutability model by counting 
#' the number of mutations occuring in the center position for all 5-mer motifs.
#'
#' @param    db                  data.frame containing sequence data.
#' @param    substitutionModel   matrix of 5-mer substitution rates built by 
#'                               \code{\link{createSubstitutionMatrix}}.
#' @param    model               type of model to create. The default model, "RS", creates 
#'                               a model by counting both replacement and silent mutations.
#'                               The "S" specification builds a model by counting only 
#'                               silent mutations.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn         name of the column containing the V-segment allele call.
#' @param    multipleMutation    string specifying how to handle multiple mutations occuring 
#'                               within the same 5-mer. If \code{"independent"} then multiple 
#'                               mutations within the same 5-mer are counted indepedently. 
#'                               If \code{"ignore"} then 5-mers with multiple mutations are 
#'                               excluded from the total mutation tally.
#' @param    minNumSeqMutations  minimum number of mutations in sequences containing each 5-mer
#'                               to compute the mutability rates. If the number is smaller 
#'                               than this threshold, the mutability for the 5-mer will be 
#'                               inferred. Default is 500.    
#' @param    returnSource        return the sources of 5-mer mutabilities (measured vs.
#'                               inferred). Default is false.                          
#'
#' @return   A named numeric vector of 1024 normalized mutability rates for each 5-mer 
#'           motif with names defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#'
#' # Create model using only silent mutations and ignore multiple mutations
#' sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="S", multipleMutation="ignore",
#'                                     minNumSeqMutations=10)
#' 
#' @export
createMutabilityMatrix <- function(db, substitutionModel, model=c("RS", "S"),
                                   sequenceColumn="SEQUENCE_IMGT", 
                                   germlineColumn="GERMLINE_IMGT_D_MASK",
                                   vCallColumn="V_CALL",
                                   multipleMutation=c("independent", "ignore"),
                                   minNumSeqMutations=500, 
                                   returnSource=FALSE) {
    # substitutionModel=sub_model; model="S"; sequenceColumn="SEQUENCE_IMGT"; germlineColumn="GERMLINE_IMGT_D_MASK"
    # vCallColumn="V_CALL"; multipleMutation="ignore"; minNumSeqMutations=10; returnSource=FALSE
    
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    db[[sequenceColumn]] = toupper(db[[sequenceColumn]])
    db[[germlineColumn]] = toupper(db[[germlineColumn]])
    
    # Check that the substitution model is valid
    if (any(dim(substitutionModel) != c(4, 1024))) {
        stop ("Please supply a valid 5-mer substitutionModel.")
    }

    # Set constants for function
    nuc_chars <- NUCLEOTIDES[1:4]
    
    # Remove IMGT gaps in the germline & input sequences
    matInputCollapsed <- removeCodonGaps(db[, c(sequenceColumn, germlineColumn)])
    # TODO: Unnecessary conversion
    db[[sequenceColumn]] <- matInputCollapsed[, 1]
    db[[germlineColumn]] <- matInputCollapsed[, 2]
    
    # Count mutations
    # TODO: this could be listMutations() instead, and skip the conversion from matInputCollapsed back to a data.frame
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn,
                                       multipleMutation=multipleMutation,
                                       model=model)
    
    # Foreground Count: Count the number of observed mutations for each 5-mer
    template <- rep(0, 1024)
    names(template) <- seqinr::words(5, nuc_chars)
    COUNT <- list()
    for(index in 1:length(mutations)){
        COUNT[[index]] <- template
        indexMutation <- mutations[[index]]
        if(!sum(is.na(indexMutation))){
            cSeq <-  s2c(db[index, sequenceColumn])
            cGL  <-  s2c(db[index, germlineColumn])
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions <= VLENGTH]
            positions <- positions[!is.na(positions)]
            for (position in  positions){
                wrd5 <- substr(db[index, germlineColumn], position - 2, position + 2)
                if(!grepl("[^ACGT]", wrd5) & nchar(wrd5) == 5){
                    codonNucs = getCodonPos(position)
                    codonGL = cGL[codonNucs]
                    codonSeq = cSeq[codonNucs]
                    muCodonPos = {position - 1} %% 3 + 1
                    #seqAtMutation <- codonSeq[muCodonPos]
                    glAtMutation <- codonGL[muCodonPos]
                    if (!any(codonGL %in% c("N", "-", ".")) & !any(codonSeq %in% c("N", "-", "."))) {
                        if (!length(grep("N", wrd5))) {
                            COUNT[[index]][wrd5]<- COUNT[[index]][wrd5] + 1;
                        }
                    }
                }
            }
        }
    }

    # Define sum of rates for nucleotide sets from substitution model
    # Two character sets
    wrd2Index <- combn(1:4, 2)
    wrd2Sums <- t(apply(wrd2Index, 2, function(x) colSums(substitutionModel[x, ], na.rm=TRUE)))
    rownames(wrd2Sums) <- apply(wrd2Index, 2, function(x) paste(nuc_chars[x], collapse=""))
    # Three character sets
    wrd3Index <- combn(1:4, 3)
    wrd3Sums <- t(apply(wrd3Index, 2, function(x) colSums(substitutionModel[x, ], na.rm=TRUE)))
    rownames(wrd3Sums) <- apply(wrd3Index, 2, function(x) paste(nuc_chars[x], collapse=""))
    # Merge single character, two character and three character sets
    substitutionSums <- rbind(substitutionModel, wrd2Sums, wrd3Sums)
        
    # Replace dots with Ns
    sSeqVec <- gsub("\\.", "N", db[[sequenceColumn]])
    sGermVec <- gsub("\\.", "N", db[[germlineColumn]])
    
    # Define template for 5-mer sums by position
    countTemplate <- matrix(0, VLENGTH, 1024, dimnames=list(1:VLENGTH, names(template)))
    
    # Background Count: Count the number of occurrences of each 5-mer
    BG_COUNT <- list()
    for (index in 1:length(mutations)) {
        tmpCounts <- countTemplate
        sGL <- sGermVec[index]
        cSeq <-  s2c(sSeqVec[index])
        cGL  <-  s2c(sGL)[1:VLENGTH]
        positions <- 3:(length(cGL) - 2)
        for (pos in  positions) {
            wrd5 <- substr(sGL, pos - 2, pos + 2)
            if (!grepl("[^ACGT]", wrd5) & nchar(wrd5) == 5) {
                codonNucs <- getCodonPos(pos)
                codonGL <- cGL[codonNucs]
                codonSeq <- cSeq[codonNucs]
                muCodonPos <- (pos - 1) %% 3 + 1
                glAtMutation <- codonGL[muCodonPos]
                if (!any(codonGL %in% c("N", "-")) & !any(codonSeq %in% c("N", "-"))) {
                    # Determine mutation types for NUCLEOTIDES[1:4]
                    muType <- CODON_TABLE[1:4 + 4*(muCodonPos - 1), stri_flatten(codonGL)]

                    # Set characters that meet mutation criteria
                    if (model == "S") {
                        muChars <- nuc_chars[1:4][nuc_chars[1:4] != glAtMutation & muType == "S"]
                    } else { 
                        muChars <- nuc_chars[1:4][nuc_chars[1:4] != glAtMutation]
                    }
                    
                    # Update counts
                    if (length(muChars) > 0) {
                        tmpCounts[pos, wrd5] <- substitutionSums[stri_flatten(muChars), wrd5]
                    }
                }
            }
        }
        BG_COUNT[[index]] <- colSums(tmpCounts)
        BG_COUNT[[index]][BG_COUNT[[index]] == 0] <- NA
    }
    
    Mutability <- list()
    for(i in 1:length(mutations)) {
        mut_mat <- COUNT[[i]] / BG_COUNT[[i]]
        mut_mat <- mut_mat / sum(mut_mat, na.rm=TRUE)
        mut_mat[!is.finite(mut_mat)] <- NA
        wgt_mat <- length(mutations[[i]])
        Mutability[[i]] <- list(mut_mat, wgt_mat)
    }
    
    
    # Aggregate mutability
    MutabilityMatrix <- sapply(Mutability, function(x) x[[1]])
    MutabilityWeights <- sapply(Mutability, function(x) x[[2]])
    Mutability_Mean <- apply(MutabilityMatrix, 1, weighted.mean, w=MutabilityWeights, na.rm=TRUE)
    Mutability_Mean[!is.finite(Mutability_Mean)] <- NA
    Mutability_Mean[Mutability_Mean == 0] <- NA
    
    # Filter out 5-mers with low number of observed mutations in the sequences
    NumSeqMutations <- sapply(1:1024,function(i)sum(MutabilityWeights[!is.na(MutabilityMatrix[i,])])) 
    Mutability_Mean[NumSeqMutations < minNumSeqMutations] <- NA
    
    
    # Infer mutability for missing 5-mers
    .fillHot <-function(FIVEMER,mutability){
        if(FIVEMER%in%names(mutability))if(!is.na(mutability[[FIVEMER]]))if(mutability[[FIVEMER]]>=0.0)return(mutability[[FIVEMER]])
        Nuc=substr(FIVEMER,3,3)
        #Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")
        FIVE=0
        COUNT=0
        
        # For A/T, infer mutability using the 3-mer model. 
        if(Nuc%in%c("A","T")){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(canMutateTo(substr(FIVEMER,1,1))[i],substr(FIVEMER,2,4),canMutateTo(substr(FIVEMER,5,5))[j],collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
        # For G, infer using 5-mers with the same downstream nucleotides 
        if(Nuc=="G"){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(canMutateTo(substr(FIVEMER,1,1))[i],canMutateTo(substr(FIVEMER,2,2))[j],substr(FIVEMER,3,5),collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
        
        # For C, infer using 5-mers with the same upstream nucleotides 
        if(Nuc=="C"){
            for(i in 1:3){
                for(j in 1:3){
                    MutatedNeighbor=paste(substr(FIVEMER,1,3),canMutateTo(substr(FIVEMER,4,4))[i],canMutateTo(substr(FIVEMER,5,5))[j],collapse="",sep="")
                    if(!is.na(mutability[[MutatedNeighbor]])){
                        FIVE=FIVE+mutability[[MutatedNeighbor]]
                        COUNT=COUNT+1
                    }
                }
            }
            return(FIVE/COUNT)
        }
    }
    
    Mutability_Mean_Complete <-sapply(words(5, nuc_chars), .fillHot, mutability = Mutability_Mean)
    
    for(i in names(which(is.na(Mutability_Mean_Complete)))){
        Mutability_Mean_Complete[i]<- .fillHot(i,mutability=Mutability_Mean_Complete)
    }
    for(i in names(which((Mutability_Mean_Complete)<1e-6))){
        Mutability_Mean_Complete[i]<- .fillHot(i,mutability=Mutability_Mean_Complete)
    }
    # If the neighboring 5-mers still don't have enough mutations, use 0 instead. 
    if (length(is.na(Mutability_Mean_Complete)) > 0) {
        warning("Insufficient number of mutations to infer some 5-mers. Filled with 0. ")
        Mutability_Mean_Complete[is.na(Mutability_Mean_Complete)] = 0 
    }
    
    
    # Normalize
    Mutability_Mean_Complete <- Mutability_Mean_Complete / sum(Mutability_Mean_Complete, na.rm=TRUE)
    
    # Return whether the 5-mer mutability is measured or inferred
    if (returnSource) {
        Mutability_Mean_Complete_Source = data.frame(Fivemer = names(Mutability_Mean_Complete),
                                                     Mutability = Mutability_Mean_Complete)
        Mutability_Mean_Complete_Source$Source = "Measured"
        Mutability_Mean_Complete_Source[Mutability_Mean_Complete_Source$Fivemer %in% 
                                            names(which(is.na(Mutability_Mean))), "Source"] = "Inferred"
        return(Mutability_Mean_Complete_Source)
    }
    
    return(Mutability_Mean_Complete)
}


#' Extends a substitution model to include Ns.
#' 
#' \code{extendSubstitutionMatrix} extends a 5-mer nucleotide substitution model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    substitutionModel  matrix of 5-mers substitution counts built by 
#'                              \code{\link{createSubstitutionMatrix}}.
#' 
#' @return   A 5x3125 matrix of normalized substitution rate for each 5-mer motif with 
#'           rows names defining the center nucleotide, one of \code{c("A", "C", "G", "T", "N")}, 
#'           and column names defining the 5-mer nucleotide sequence.
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#'
#' # Create model using only silent mutations and ignore multiple mutations
#' sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
#' ext_model <- extendSubstitutionMatrix(sub_model)
#' 
#' @export
extendSubstitutionMatrix <- function(substitutionModel) {
    # TODO: fix order so Ns are at the end? (c(input_names, words not in input_names))
    
    # Define old and new column/row names
    input_names <- colnames(substitutionModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- matrix(NA, nrow=length(nuc_chars), ncol=length(nuc_5mers), 
                         dimnames=list(nuc_chars, nuc_5mers))
    
    # Extend matrix with Ns
    for (mer in nuc_5mers) {
        if (mer %in% input_names) {
            extend_mat[, mer] <- c(substitutionModel[, mer], "N"=NA)
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[, mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[, mer] <- c(apply(substitutionModel[, mer_index], 1, mean, na.rm=TRUE), "N"=NA)
            }
        }
    }
    
    # Normalize
    #extend_mat <- apply(extend_mat, 2, function(x) { x/sum(x, na.rm=TRUE) })
    extend_mat[!is.finite(extend_mat)] <- NA
    
    return (extend_mat)
}


#' Extends a mutability model to include Ns.
#' 
#' \code{extendMutabilityMatrix} extends a 5-mer nucleotide mutability model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    mutabilityModel  vector of 5-mer mutability rates built by 
#'                            \code{\link{createMutabilityMatrix}}.
#' 
#' @return   A 3125 vector of normalized mutability rates for each 5-mer motif with 
#'           names defining the 5-mer nucleotide sequence.
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#'
#' # Create model using only silent mutations and ignore multiple mutations
#' sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="S", multipleMutation="ignore",
#'                                     minNumSeqMutations=10)
#' ext_model <- extendMutabilityMatrix(mut_model)
#' 
#' @export
extendMutabilityMatrix <- function(mutabilityModel) {
    # TODO: fix order so Ns are at the end? (c(input_names, words not in input_names))
    
    # Define old and new column/row names
    input_names <- names(mutabilityModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- array(NA, dim=length(nuc_5mers), dimnames=list(nuc_5mers))
    
    # Extend matrix with Ns
    for(mer in nuc_5mers) {
        #cat(mer,"\n")
        if (mer %in% input_names) {
            extend_mat[mer] <- mutabilityModel[mer]
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[mer] <- mean(mutabilityModel[mer_index], na.rm=TRUE)
            }
        }
    }
    
    # Normalize    
    #extend_mat <- extend_mat / sum(extend_mat, na.rm=TRUE)
    extend_mat[!is.finite(extend_mat)] <- NA
    
    return(extend_mat)
}
 

#' Calculates a targeting rate matrix
#' 
#' \code{createTargetingMatrix} calculates the targeting model matrix as the
#' combined probability of mutability and substitution.
#'
#' @param    substitutionModel  matrix of 5-mers substitution rates built by 
#'                              \code{\link{createSubstitutionMatrix}} or 
#'                              \code{\link{extendSubstitutionMatrix}}.
#' @param    mutabilityModel    vector of 5-mers mutability rates built by 
#'                              \code{\link{createMutabilityMatrix}} or 
#'                              \code{\link{extendMutabilityMatrix}}.
#' 
#' @return   A matrix with the same dimensions as the input \code{substitutionModel} 
#'           containing normalized targeting probabilities for each 5-mer motif with 
#'           row names defining the center nucleotide and column names defining the 
#'           5-mer nucleotide sequence.
#' 
#' @details
#' Targeting rates are calculated by multiplying the normalized mutability rate by the 
#' normalized substitution rates for each individual 5-mer.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#'
#' # Create 4x1024 model using only silent mutations and ignore multiple mutations
#' sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="S", multipleMutation="ignore",
#'                                     minNumSeqMutations=10)
#' tar_model <- createTargetingMatrix(sub_model, mut_model)
#' 
#' # Create 5x3125 model including Ns
#' sub_model <- extendSubstitutionMatrix(sub_model)
#' mut_model <- extendMutabilityMatrix(mut_model)
#' tar_model <- createTargetingMatrix(sub_model, mut_model)
#' 
#' @export
createTargetingMatrix <- function(substitutionModel, mutabilityModel) {
    # Calculate targeting
    tar_mat <- sweep(substitutionModel, 2, mutabilityModel, `*`)
    
    # Normalize    
    #tar_mat <- tar_mat / sum(tar_mat, na.rm=TRUE)
    tar_mat[!is.finite(tar_mat)] <- NA
    
    return(tar_mat)
}


#' Creates a TargetingModel
#' 
#' \code{createTargetingModel} creates a 5-mer \code{TargetingModel}.
#'
#' @param    db                  data.frame containing sequence data.
#' @param    model               type of model to create. The default model, "RS", creates 
#'                               a model by counting both replacement and silent mutations.
#'                               The "S" specification builds a model by counting only 
#'                               silent mutations.
#' @param    sequenceColumn      name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn      name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn         name of the column containing the V-segment allele calls.
#' @param    multipleMutation    string specifying how to handle multiple mutations occuring 
#'                               within the same 5-mer. If \code{"independent"} then multiple 
#'                               mutations within the same 5-mer are counted indepedently. 
#'                               If \code{"ignore"} then 5-mers with multiple mutations are 
#'                               excluded from the otal mutation tally.
#' @param    minNumMutations     minimum number of mutations required to compute the 5-mer 
#'                               substitution rates. If the number of mutations for a 5-mer
#'                               is below this threshold, its substitution rates will be 
#'                               estimated from neighboring 5-mers. Default is 50.   
#' @param    minNumSeqMutations  minimum number of mutations in sequences containing each 5-mer
#'                               to compute the mutability rates. If the number is smaller 
#'                               than this threshold, the mutability for the 5-mer will be 
#'                               inferred. Default is 500.   
#' @param    modelName           name of the model.
#' @param    modelDescription    description of the model and its source data.
#' @param    modelSpecies        genus and species of the source sequencing data.
#' @param    modelDate           date the model was built. If \code{NULL} the current date
#'                               will be used.
#' @param    modelCitation       publication source.
#' 
#' @return   A \code{\link{TargetingModel}} object.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data.
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  See \link{TargetingModel} for the return object.
#' @family   targeting model functions
#' 
#' @examples
#' # Subset example data to one isotype and sample as a demo
#' db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#'
#' # Create model using only silent mutations and ignore multiple mutations
#' model <- createTargetingModel(db, model="S", multipleMutation="ignore")
#' 
#' @export
createTargetingModel <- function(db, model=c("RS", "S"), sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 vCallColumn="V_CALL",
                                 multipleMutation=c("independent", "ignore"),
                                 minNumMutations=50, minNumSeqMutations=500,
                                 modelName="", modelDescription="", modelSpecies="", 
                                 modelCitation="", modelDate=NULL) {
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn, vCallColumn))
    if (check != TRUE) { stop(check) }
    
    # Set date
    if (is.null(modelDate)) { modelDate <- format(Sys.time(), "%Y-%m-%d") }

    # Create models
    sub_mat<- createSubstitutionMatrix(db, model=model, 
                                       sequenceColumn=sequenceColumn,
                                       germlineColumn=germlineColumn,
                                       vCallColumn=vCallColumn,
                                       multipleMutation=multipleMutation,
                                       minNumMutations=minNumMutations,
                                       returnModel="5mer")
    mut_mat <- createMutabilityMatrix(db, sub_mat, model=model,
                                      sequenceColumn=sequenceColumn,
                                      germlineColumn=germlineColumn,
                                      vCallColumn=vCallColumn,
                                      multipleMutation=multipleMutation,
                                      minNumSeqMutations=minNumSeqMutations,
                                      returnSource=FALSE)

    # Extend 5-mers with Ns
    sub_mat <- extendSubstitutionMatrix(sub_mat)
    mut_mat <- extendMutabilityMatrix(mut_mat)
    
    tar_mat <- createTargetingMatrix(sub_mat, mut_mat) 
    
    # Define TargetingModel object
    model_obj <- new("TargetingModel",
                     name=modelName,
                     description=modelDescription,
                     species=modelSpecies,
                     date=modelDate,
                     citation=modelCitation,
                     substitution=sub_mat,
                     mutability=mut_mat,
                     targeting=tar_mat)

    return(model_obj)
}


#' Calculates a 5-mer distance matrix from a TargetingModel object
#' 
#' \code{calcTargetingDistance} converts the targeting rates in a TargetingModel model 
#' to a matrix of 5-mer to single-nucleotide mutation distances.
#' 
#' @param    model     \link{TargetingModel} object with mutation likelihood information.
#'                                                
#' @return   A matrix of distances for each 5-mer motif with rows names defining 
#'           the center nucleotide and column names defining the 5-mer nucleotide 
#'           sequence.
#'           
#' @details
#' The targeting model is transformed into a distance matrix by:
#' \enumerate{
#'   \item  Converting the likelihood of being mutated \eqn{p=mutability*substitution} to 
#'          distance \eqn{d=-log10(p)}.
#'   \item  Dividing this distance by the mean of the distances
#'   \item  Converting all infinite, no change (e.g., A->A), and NA distances to 
#'          zero.
#' }
#' 
#' @seealso  Takes as input a \link{TargetingModel} object.
#' @family   targeting model functions
#' 
#' @examples
#' # Calculate targeting distance of HS5FModel
#' dist <- calcTargetingDistance(HS5FModel)
#' 
#' @export
calcTargetingDistance <- function(model) {
    if (is(model, "TargetingModel")) {
        model <- model@targeting
    } else if (!is(model, "matrix")) {
        stop("Input must be either a targeting matrix or TargetingModel object.")
    }
    
    # Take negative log10 of all probabilities
    model_dist <- -log10(model)
    # Center distances on 1 (divide by mean)
    model_dist <- model_dist/mean(model_dist, na.rm=T)
    # Set infinite values to NA
    model_dist[!is.finite(model_dist)] <- NA
    
    # TODO: the indexing of A-->A etc positions can probably be done smarter/faster
    # Assign no-change entries to distance 0
    center_nuc <- gsub("..([ACGTN])..", "\\1", colnames(model_dist))
    for (i in 1:length(center_nuc)) {
        model_dist[center_nuc[i], i] <- 0
    }

    return(model_dist)
}


# Rescales mutability probabilities from a TargetingModel
# 
# \code{rescaleMutability} renormalizes the mutability probabilities
# in a TargetingModel model and returns a rescaled matrix of mutability scores.
# 
# @param    model     \link{TargetingModel} object with mutation likelihood information.
# @param    mean      the mean value for the rescaled mutability scores.
#                                                
# @return   A named vector of mutability scores for each 5-mer motif with mean
#           equal to \code{mean}.
#           
# @seealso  Takes as input a \link{TargetingModel} object.
# @family   targeting model functions
# 
# @examples
# # Subset example data to one isotype and sample as a demo
# db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")
#
# # Create model and rescale mutabilities
# model <- createTargetingModel(db, model="S", multipleMutation="ignore")
# mut <- rescaleMutability(model)
#
rescaleMutability <- function(model, mean=1.0) {
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    
    # TODO:  perhaps this is not so useful.  or perhaps it should be relative to max(model).
    # Renormalize
    rescaled <- model / sum(model, na.rm=T) * sum(!is.na(model)) * mean
    rescaled[!is.finite(rescaled)] <- NA
    
    return(rescaled)
}


# Remove in-frame IMGT gaps
# 
# @param    matInput  Nx2 matrix with input and germline sequences in each column
# @return   A two column matrix with "..." codons removed.
removeCodonGaps <- function(matInput) {
    # Function to return valid codon sets
    # i = position, x = codon list 1, y = codon list 2
    .f1 <- function(i, x, y) {
        xcod <- x[i]
        ycod <- y[i]
        if (xcod != "..." & ycod != "...") {
            c(xcod, ycod)
        } else {
            c("", "")
        }
    }

    # Function to parse sequences
    # z = vector of 2 sequences
    .f2 <- function(z) {
        # Split strings into 3-mers
        cods <- stri_extract_all_regex(c(z[1], z[2]), ".{3}")
        cmat <- sapply(1:length(cods[[1]]), .f1, x=cods[[1]], y=cods[[2]])
        apply(cmat, 1, paste, collapse="")
    }
    
    matCollapsed <- t(apply(matInput, 1, .f2))
    
    return(matCollapsed)
}


#### I/O Functions ####

#' Write targeting model distances to a file
#' 
#' \code{writeTargetingDistance} writes a 5-mer targeting distance matrix 
#' to a tab-delimited file.
#' 
#' @param    model     \code{\link{TargetingModel}} object with 
#'                     mutation likelihood information.
#' @param    file      name of file to write.
#'                                                
#' @return   NULL
#'           
#' @details
#' The targeting distance write as a tab-delimited 5x3125 matrix. Rows define the mutated 
#' nucleotide at the center of each 5-mer, one of \code{c("A", "C", "G", "T", "N")}, 
#' and columns define the complete 5-mer of the unmutated nucleotide sequence. 
#' \code{NA} values in the distance matrix are replaced with distance 0.
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object and calculates  
#'           distances using \link{calcTargetingDistance}.
#' @family   targeting model functions
#' 
#' @examples
#' \dontrun{
#' # Write HS5F targeting model to working directory as hs5f.tab
#' writeTargetingModel(HS5FModel, "hs5f.tab") 
#' }
#' 
#' @export
writeTargetingDistance <- function(model, file) {
    to_write <- as.data.frame(calcTargetingDistance(model))
    to_write[is.na(to_write)] <- 0
    write.table(to_write, file, quote=FALSE, sep="\t")
}


#### Plotting functions ####

#' Plot mutability probabilities
#' 
#' \code{plotMutability} plots the mutability rates of a \code{TargetingModel}.
#' 
#' @param    model        \link{TargetingModel} object or matrix containing normalized 
#'                        mutability rates.
#' @param    nucleotides  vector of center nucleotide characters to plot.
#' @param    style        type of plot to draw. One of:
#'                        \itemize{
#'                          \item \code{"hedgehog"}:  circular plot showing higher mutability
#'                                                    scores further from the circle. The 5-mer
#'                                                    is denoted by the values of the inner 
#'                                                    circle. The 5-mer is read from the most interior 
#'                                                    position of the 5-mer (5') to most exterior position 
#'                                                    (3'), with the center nucleotide in the center ring.
#'                                                    Note, the order in which the 5-mers are plotted is
#'                                                    different for nucleotides \code{c("A", "C")} and 
#'                                                    \code{c("G", "T")}.
#'                          \item \code{"bar"}:       bar plot of mutability similar to the 
#'                                                    \code{hedgehog} style with the most 5' positions
#'                                                    of each 5-mer at the base of the plot.
#'                        }
#' @param    size         numeric scaling factor for lines and text in the plot.
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                        objects; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#'                                                
#' @return   A named list of ggplot objects defining the plots, with names defined by the 
#'           center nucleotide for the plot object.
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' @family   targeting model functions
#' 
#' @examples
#' # Plot one nucleotide in circular style
#' plotMutability(HS5FModel, "C")
#' 
#' # Plot two nucleotides in barchart style
#' plotMutability(HS5FModel, c("G","T"), style="bar")
#' 
#' @export
plotMutability <- function(model, nucleotides=c("A", "C", "G", "T"), 
                           style=c("hedgehog", "bar"), size=1, silent=FALSE, 
                           ...) {
    # model=HS5FModel
    # nucleotides=c("C")
    # nucleotides=c("A", "C", "G", "T")
    # style="hedgehog"
    # size=1
    # silent=FALSE
    
    # Check input
    nucleotides <- toupper(nucleotides)
    style <- match.arg(style)
    
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    
    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.margin=grid::unit(0, "lines"),
              panel.background=element_blank()) +
        theme(axis.text=element_text(margin=grid::unit(0, "lines"))) +
        theme(text=element_text(size=10*size),
              title=element_text(size=10*size),
              legend.margin=grid::unit(0, "lines"),
              legend.background=element_blank())

    # Scaling and layout parameters
    score_offset <- 0
    score_scale <- 15
    text_offset <- -5.6
    
    # Set guide colors
    motif_colors <- setNames(c("#4daf4a", "#e41a1c", "#094d85", "#999999"),
                             c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    dna_colors <- setNames(c("#7bce77", "#ff9b39", "#f04949", "#5796ca", "#c4c4c4"),
                           c("A", "C", "G", "T", "N"))
    
    # Build data.frame of mutability scores
    mut_scores <- model[!grepl("N", names(model))]
    mut_scores[!is.finite(mut_scores)] <- 0
    mut_words <- names(mut_scores)
    mut_positions <- as.data.frame(t(sapply(mut_words, seqinr::s2c)))
    colnames(mut_positions) <- paste0("pos", 1:ncol(mut_positions))
    mut_df <- data.frame(word=mut_words, 
                         score=mut_scores, 
                         mut_positions)

    # Add hot and cold-spot motif information
    mut_df$motif <- "Neutral"
    mut_df$motif[grepl("(.[AT]A..)|(..T[AT].)", mut_df$word, perl=TRUE)] <- "WA/TW"
    mut_df$motif[grepl("([AT][GA]C..)|(..G[CT][AT])", mut_df$word, perl=TRUE)] <- "WRC/GYW"
    mut_df$motif[grepl("([CG][CT]C..)|(..G[GA][CG])", mut_df$word, perl=TRUE)] <- "SYC/GRS"
    mut_df$motif <- factor(mut_df$motif, levels=c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    
    # Subset to nucleotides of interest
    mut_df <- mut_df[mut_df$pos3 %in% nucleotides, ]

    # Functions to transform and revert score for plotting
    score_max <- max(mut_df$score, na.rm=TRUE)
    .transform_score <- function(x) { x / score_max * score_scale + score_offset }
    .invert_score <- function(x) { (x - score_offset) / score_scale * score_max }
    
    # Rescale scores for plotting
    mut_df$score <- .transform_score(mut_df$score)
    
    # Build plots for each center nucleotide
    plot_list <- list()
    for (center_nuc in nucleotides) {
        # center_nuc <- "C"
        # Subset data to center nucleotide
        sub_df <- mut_df[mut_df$pos3 == center_nuc, ]
        
        # Order 5-mers by positions, with reversed order if center nucleotide is G or T
        if (center_nuc %in% c("A", "C")) {
            sub_df <- dplyr::arrange_(sub_df, .dots=c("pos1", "pos2", "pos4", "pos5"))
            sub_df$x <- 1:nrow(sub_df)            
        } else if (center_nuc %in% c("G", "T")) {
            sub_df <- dplyr::arrange_(sub_df, .dots=c("pos5", "pos4", "pos2", "pos1"))
            sub_df$x <- 1:nrow(sub_df)
        } else {
            stop("Invalid nucleotide choice")
        }
        
        # Melt 5-mer position data
        sub_melt <- sub_df %>% 
            tidyr::gather_("pos", "char", colnames(mut_positions)) %>% 
            select_(.dots=c("x", "pos", "char"))
        #sub_melt$pos <- factor(sub_melt$pos, levels=mut_names)
        #sub_melt$pos <- as.numeric(sub_melt$pos)
        sub_melt$pos <- as.numeric(gsub("pos", "", sub_melt$pos))

        # Define nucleotide text and rectangle position data
        sub_text <- list()
        for (i in 1:5) {
            nuc_rle <- rle(sub_melt$char[sub_melt$pos == i])

            # Set rectangle x limits
            rect_max <- cumsum(nuc_rle$lengths)
            rect_min <- rect_max - diff(c(0, rect_max))
            
            # Set text position
            if (length(rect_max) > 1) { 
                text_x <- rect_max - diff(c(0, rect_max)) / 2
            } else { 
                text_x <- rect_max / 2 
            }
        
            tmp_df <- data.frame(text_x=text_x, 
                                 text_y=i,
                                 text_label=factor(nuc_rle$values, levels=names(dna_colors)),
                                 rect_min=rect_min,
                                 rect_max=rect_max)
            
            sub_text[[i]] <- tmp_df
        }

        # Define text and rectangle positions for inner circle
        sub_melt$pos <- sub_melt$pos + text_offset
        sub_text <- lapply(sub_text, function(x) { dplyr::mutate_(x, text_y=interp(~ y + text_offset, y=as.name("text_y"))) })
        sub_rect <- dplyr::bind_rows(sub_text) %>%
            mutate_(rect_width=interp(~ y - x, x=as.name("rect_min"), y=as.name("rect_max")),
                    ymin=interp(~ y - 0.5, y=as.name("text_y")),
                    ymax=interp(~ y + 0.5, y=as.name("text_y")))
        

        # Define base plot object
        p1 <- ggplot(sub_df) + 
            base_theme + 
            #ggtitle(paste0("NN", center_nuc, "NN")) +
            xlab("") +
            ylab("") + 
            scale_color_manual(name="Motif", values=c(motif_colors, dna_colors), 
                               breaks=names(motif_colors)) +
            scale_fill_manual(name="", values=c(motif_colors, dna_colors), guide=FALSE) +
            geom_rect(data=sub_rect, 
                      mapping=aes_string(xmin="rect_min", xmax="rect_max", ymin="ymin", ymax="ymax", 
                                         fill="text_label", color="text_label"), 
                      size=0.5*size, alpha=1, show.legend=FALSE) +
            #geom_tile(data=sub_rect, 
            #          mapping=aes_string(x="text_x", y="text_y", width="rect_width", fill="text_label"), 
            #          size=0, alpha=0.7, show.legend=FALSE) +
            #geom_tile(data=sub_melt, mapping=aes_string(x="x", y="pos", fill="char"), size=0, alpha=0.7,
            #          show.legend=FALSE) +
            geom_text(data=sub_text[[3]], mapping=aes_string(x="text_x", y="text_y", label="text_label"), 
                      color="black", hjust=0.5, vjust=0.5, size=3*size, fontface=2)
        
        # Add 5-mer nucleotide labels
        if (center_nuc %in% c("A", "C")) {
            p1 <- p1 + geom_text(data=sub_text[[1]], mapping=aes_string(x="text_x", y="text_y", label="text_label"), 
                                 color="black", hjust=0.5, vjust=0.5, size=2*size) +
                geom_text(data=sub_text[[2]], mapping=aes_string(x="text_x", y="text_y", label="text_label"), 
                          color="black", hjust=0.5, vjust=0.5, size=2*size)
        } else if (center_nuc %in% c("G", "T")) {
            p1 <- p1 + geom_text(data=sub_text[[4]], mapping=aes_string(x="text_x", y="text_y", label="text_label"), 
                                 color="black", hjust=0.5, vjust=0.5, size=2*size) +
                geom_text(data=sub_text[[5]], mapping=aes_string(x="text_x", y="text_y", label="text_label"), 
                          color="black", hjust=0.5, vjust=0.5, size=2*size)
        }

        # Add style options and mutability scores
        if (style == "hedgehog") {
            y_limits <- c(text_offset - 1, score_scale + score_offset)
            #orient_x <- sub_text[[3]]$text_x[1]
            #orient_y <- text_offset - 1
            p1 <- p1 + theme(plot.margin=grid::unit(c(0, 0, 0, 0), "lines"),
                             panel.grid=element_blank(), 
                             panel.border=element_blank(),
                             axis.title=element_blank(),
                             axis.text=element_blank(), 
                             axis.ticks=element_blank(),
                             legend.direction="horizontal",
                             legend.justification=c(0.5, 1),
                             legend.position=c(0.5, 1)) +
                guides(color=guide_legend(override.aes=list(linetype=1, size=2*size))) +
                scale_x_continuous(expand=c(0, 0)) +
                scale_y_continuous(limits=y_limits, expand=c(0, 0)) +
                coord_polar(theta="x") +
                geom_segment(data=sub_df, mapping=aes_string(x="x", xend="x", yend="score", color="motif"), 
                             y=score_offset, size=0.75*size)
        } else if (style == "bar") {
            y_breaks <- seq(score_offset, score_scale + score_offset, 1)
            y_limits <- c(text_offset + 0.5, score_scale + score_offset)
            sub_colors <- motif_colors[names(motif_colors) %in% sub_df$motif]
            p1 <- p1 + theme(plot.margin=grid::unit(c(1, 1, 1, 1), "lines"),
                             panel.grid=element_blank(), 
                             panel.border=element_rect(color="black"),
                             axis.text.x=element_blank(), 
                             axis.ticks.x=element_blank(),
                             legend.position="top") +
                guides(color=guide_legend(override.aes=list(fill=sub_colors, linetype=0))) +
                ylab("Mutability") +
                scale_x_continuous(expand=c(0, 1)) +
                scale_y_continuous(limits=y_limits, breaks=y_breaks, expand=c(0, 0.5),
                                   labels=function(x) scales::scientific(.invert_score(x))) +
                geom_bar(data=sub_df, mapping=aes_string(x="x", y="score", fill="motif", color="motif"), 
                         stat="identity", position="identity", size=0, width=0.7)
        }

        # Add additional theme elements
        p1 <- p1 + do.call(theme, list(...))
        
        # Add plot to list
        plot_list[[center_nuc]] <- p1
    }

    
    # Plot
    if (!silent) { 
        do.call(multiggplot, args=c(plot_list, ncol=length(plot_list))) 
    }
    
    invisible(plot_list)
}


#### Original BASELINe functions ####

# Given a nuc, returns the other 3 nucs it can mutate to
canMutateTo <- function(nuc) {
    NUCLEOTIDES[1:4][NUCLEOTIDES[1:4] != nuc]
}


# Compute the mutations types
mutationTypeOptimized <- function(matOfCodons) {
    apply(matOfCodons, 1, function(x) { mutationType(x[2], x[1]) })
}


# row 1 = GL
# row 2 = Seq
# in_matrix <- matrix(c(c("A","A","A","C","C","C"), c("A","G","G","C","C","A")), 2 ,6, byrow=T)
# analyzeMutations2NucUri(in_matrix)
analyzeMutations2NucUri <- function(in_matrix) {
    if(ncol(in_matrix) > VLENGTH) {
        paramGL = in_matrix[2,1:VLENGTH]
        paramSeq = in_matrix[1,1:VLENGTH]
    } else {
        paramGL = in_matrix[2,]
        paramSeq = in_matrix[1,]
    }
    #mutations = apply(rbind(paramGL,paramSeq), 2, function(x){!x[1]==x[2]})
    mutations_val = paramGL != paramSeq
    if (any(mutations_val)) {
        mutationPos = {1:length(mutations_val)}[mutations_val]
        #mutationPos = mutationPos[sapply(mutationPos, function(x){!any(paramSeq[getCodonPos(x)]=="N")})]
        length_mutations =length(mutationPos)
        mutationInfo = rep(NA,length_mutations)
        if (any(mutationPos)) {
            pos<- mutationPos
            pos <- pos[!is.na(pos)]
            pos_array <- array(sapply(pos,getCodonPos))
            codonGL <- paramGL[pos_array]
            codonGL[is.na(codonGL)] <- "N"
            codonSeq = sapply(pos,function(x){
                seqP = paramGL[getCodonPos(x)]
                muCodonPos = {x-1}%%3+1
                seqP[muCodonPos] = paramSeq[x]
                return(seqP)
            })
            codonSeq[is.na(codonSeq)] <- "N"
            GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
            Seqcodons =   apply(codonSeq,2,c2s)
            mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
            names(mutationInfo) = mutationPos
        }
        if (any(!is.na(mutationInfo))) {
            return(mutationInfo[!is.na(mutationInfo)])
        } else {
            return(NA)
        }
        
        
    } else {
        return (NA)
    }
}



# List mutations
listMutations <- function(seqInput, seqGL, multipleMutation, model) {
    #if( is.na(c(seqInput, seqGL)) ) return(array(NA,4))
    if (is.na(seqInput) | is.na(seqGL)) { return(NA) }
    seqI = s2c(seqInput)
    seqG = s2c(seqGL)
    matIGL = matrix(c(seqI, seqG), ncol=length(seqI), nrow=2, byrow=T)
    mutations <- analyzeMutations2NucUri(matIGL)
    mutations <- mutations[!is.na(mutations)]
    positions <- as.numeric(names(mutations))
    # mutations <- mutations[positions <= VLENGTH]
    
    #remove the nucleotide mutations in the codons with multiple mutations
    if (multipleMutation == "ignore") {
       mutationCodons = getCodonNumb(as.numeric(names(mutations)))
       tableMutationCodons <- table(mutationCodons)
       codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
       mutations <- mutations[!(mutationCodons %in% codonsWithMultipleMutations)]
    }
    if (model == "S") {
       mutations <- mutations[mutations == "S"]
    }
    if (length(mutations) > 0) {
        return(mutations)
    } else {
        return(NA)
    }
}


# List the numbers of observed mutations
#
# This lists the observed number of mutation.
#
# @param   db  a data.frame of the DB file.
# @param   sequenceColumn  The name of the sequence column.
# @param   germlineColumn  The name of the germline column.
# 
# @return  list of mutations in each clone
listObservedMutations <- function(db, sequenceColumn="SEQUENCE_IMGT", 
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  multipleMutation=c("independent", "ignore"),
                                  model = c("RS", "S"))  {
    
    # Make sure the columns specified exist 
    if (!(sequenceColumn %in% names(db))) {
        stop("The sequence column", sequenceColumn, "was not found.")
    } 
    if (!(germlineColumn %in% names(db))) {
        stop("The germline column", germlineColumn, "was not found.")
    } 
    
    mutations <- mapply(listMutations, db[[sequenceColumn]], db[[germlineColumn]], 
                        multipleMutation, model, USE.NAMES=FALSE)
    return(mutations)
}


#### Testing functions ####

# Function to make dummy data for testing targetting functions
# 
# @param   nseq  number of sequences
# @param   nmut  number of mutations per sequence
# @param   nmer  number of 5-mers per sequence (sequence length = 5 * nmer)
#
# @return  a data.frame with columns SEQUENCE_ID, SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK, V_CALL.
# 
# @examples
# db <- makeTargetingTestDb(500)
makeTargetingTestDb <- function(nseq=10, nmut=40, nmers=50) {
    nuc_chars <- c("A", "C", "G", "T")
    
    .mut <- function(x, n) {
       i <- sample(1:nchar(x), n) 
       y <- seqinr::s2c(x)
       y[i] <- sapply(y[i], function(z) sample(nuc_chars[nuc_chars != z], 1))
       return(seqinr::c2s(y))
    }
    
    seq <- apply(replicate(nseq, sample(seqinr::words(5, nuc_chars), nmers)), 2, paste, collapse="")
    germ <- sapply(seq, .mut, n=nmut)
    db <- data.frame(SEQUENCE_ID=paste0("SEQ", 1:nseq),
                     SEQUENCE_IMGT=seq,
                     GERMLINE_IMGT_D_MASK=germ,
                     V_CALL="Homsap IGHV3-66*02 F", stringsAsFactors=FALSE)
    rownames(db) <- NULL
    
    return(db)
}
