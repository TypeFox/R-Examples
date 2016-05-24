# Mutation profiling

#' @include shazam.R
NULL

#### Clonal Consensus building functions ####

#' Identifies clonal consensus sequences
#'
#' Identifies effective/consensus sequences collapsed by clone
#'
#' \code{collapseByClone} identifies the consensus sequence of each clonal 
#' group and appends a column to the input \code{data.frame} containing the clonal 
#' consensus for each sequence.
#'
#' @param   db                  \code{data.frame} containing sequence data.
#' @param   cloneColumn         \code{character} name of the column containing clonal 
#'                              identifiers.
#' @param   sequenceColumn      \code{character} name of the column containing input 
#'                              sequences.
#' @param   germlineColumn      \code{character} name of the column containing germline 
#'                              sequences.
#' @param   expandedDb          \code{logical} indicating whether or not to return the 
#'                              expanded \code{db}, containing all the sequences (as opposed
#'                              to returning just one sequence per clone collapsed by )
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences.
#' @param   nonTerminalOnly     \code{logical} indicating whether to include mutations
#'                              at the leaves.
#' @param   nproc               Number of cores to distribute the operation over. If the 
#'                              \code{cluster} has already been set earlier, then pass the 
#'                              \code{cluster}. This will ensure that it is not reset.
#'                              
#' 
#' @return   A modified \code{db} data.frame with clonal consensus sequences in the
#'           CLONAL_SEQUENCE column.
#'
#' @details
#' For sequences identified to be part of the same clone, this function defines an 
#' effective sequence that will be representative for all mutations in the clone. Each 
#' position in this consensus (or effective) sequence is created by a weighted sampling 
#' of each mutated base (and non "N", "." or "-" characters) from all the sequences in 
#' the clone. 
#' 
#' For example, in a clone with 5 sequences that have a C at position 1, and 5 sequences
#' with a T at this same position, the consensus sequence will have a C 50\%  and T 50\% 
#' of the time it is called.
#' 
#' The function returns an updated \code{db} that collpases all the sequences by clones 
#' defined in the \code{cloneColumn} column argument.
#' 
#' Non-terminal branch mutations are defined as the set of mutations that occur on 
#' branches of the lineage tree that are not connected to a leaf. For computational 
#' efficiency, the set of non-terminal branch mutations is approximated as those that are
#' shared between more than one sequence in a clone. In this case the terminal branch 
#' mutations are filtered out.
#' 
#' This function can be parallelized if \code{db} contains thousands of sequences. 
#' Specify the number of cores available using the \code{nproc} parameter.
#' 
#' @seealso
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
#'                           BARCODE %in% c("RL016","RL018","RL019","RL021"))
#' 
#' # Run collapseByClone
#' db_new <- collapseByClone(db, cloneColumn="CLONE", 
#'                           sequenceColumn="SEQUENCE_IMGT",
#'                           germlineColumn="GERMLINE_IMGT_D_MASK",
#'                           expandedDb=FALSE,
#'                           regionDefinition=IMGT_V_NO_CDR3,
#'                           nproc=1)
#' 
#' @export
collapseByClone <- function(db, 
                            cloneColumn="CLONE", 
                            sequenceColumn="SEQUENCE_IMGT",
                            germlineColumn="GERMLINE_IMGT_D_MASK",
                            expandedDb=FALSE,
                            regionDefinition=NULL,
                            nonTerminalOnly=FALSE,
                            nproc=1) {
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL

    # Check for valid columns
    check <- checkColumns(db, c(cloneColumn, sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    if (class(expandedDb) != "logical") {
        stop ("expandedDb must be TRUE or FALSE.")
    }
    
    if (class(nonTerminalOnly) != "logical") {
        stop ("nonTerminalOnly must be TRUE or FALSE.")
    }
    
    # TODO: Why is this converted to numeric?  Won't that break?
    db[[cloneColumn]] <- as.numeric(db[[cloneColumn]])
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of CLONES, instead of doing dplyr
    dt <- data.table(db)
    setkeyv(dt, cloneColumn)
    # Get the group indexes
    #dt <- dt[, list(yidx=list(.I)), by=list(CLONE)]
    dt <- dt[, list(yidx=list(.I)), by=cloneColumn]
    groups <- dt[, yidx]
    lenGroups <- length(groups)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        if (nproc != 0) { 
            #cluster <- makeCluster(nproc, type="SOCK") 
            cluster <- parallel::makeCluster(nproc, type= "PSOCK")
        }
        parallel::clusterExport(cluster, list('db', 
                                              'sequenceColumn', 'germlineColumn', 'cloneColumn',
                                              'regionDefinition', 'calcClonalConsensus',
                                              'groups', 'c2s', 's2c', 'words', 'translate'), 
                                envir=environment() )
        registerDoParallel(cluster)
    }
    
    # Printing status to console
    cat("Collapsing clonal sequences...\n")
    
    list_ClonalConsensus <-
        foreach(idx=iterators::icount(lenGroups), .combine=c, .verbose=FALSE, .errorhandling='pass') %dopar% {
            calcClonalConsensus(inputSeq = db[[sequenceColumn]][groups[[idx]]],
                                germlineSeq = db[[germlineColumn]][groups[[idx]]],
                                regionDefinition = regionDefinition, 
                                nonTerminalOnly = nonTerminalOnly)
        }
    
    # Stop cluster
    if(nproc > 1) { #
        parallel::stopCluster(cluster) 
    }
    
    # If expandedDb is FALSE then collapse the db by clones
    if(expandedDb==FALSE){ 
        uniqueCloneIDs <-  unique(db[[cloneColumn]])
        indexOfFirstOccurenceOfClone <- match(uniqueCloneIDs, db[[cloneColumn]])
        db_ClonalConsensus <- db[indexOfFirstOccurenceOfClone, ]
        db_ClonalConsensus$CLONAL_SEQUENCE <- unlist(list_ClonalConsensus)
    }else{
        # Match the ClonalConsensus to all the sequences in the clone
        vec_ClonalConsensus <- unlist(list_ClonalConsensus)
        expanded_ClonalConsensus <- tapply(db[[cloneColumn]], 1:nrow(db),function(x){return(vec_ClonalConsensus[x])})
        db_ClonalConsensus <- db
        db_ClonalConsensus$CLONAL_SEQUENCE <- unlist(expanded_ClonalConsensus)
    }
    
    return(db_ClonalConsensus)
}



# Helper function for calcDBClonalConsensus
calcClonalConsensus <- function(inputSeq, germlineSeq, 
                                regionDefinition=NULL, 
                                nonTerminalOnly=FALSE) {
    # If only one sequence in clone, return it
    if(length(inputSeq) == 1) {
        return(inputSeq)
    }
    
    # Find length of shortest input sequence
    # This is used to trim all the sequences to that length
    # or, if a regionDefinition is passed, then only analyze till the end of the defined length
    len_inputSeq <- sapply(inputSeq, function(x){nchar(x)})
    len_shortest <- min(len_inputSeq, na.rm=TRUE)
    if(!is.null(regionDefinition)) {
        len_shortest <- min(len_shortest, regionDefinition@seqLength, na.rm=TRUE)
    }
    
    if (class(nonTerminalOnly) != "logical") {
        stop ("nonTerminalOnly must be TRUE or FALSE.")
    }
    
    #Find the length of the longest germline sequence
    len_germlineSeq <- sapply(germlineSeq, function(x){nchar(x)})
    len_longest <- max(len_germlineSeq, na.rm=TRUE)
    germlineSeq <- germlineSeq[(which(len_longest==len_germlineSeq))[1]]    
    
    # Identify the consensus sequence
    inputSeq <- unique(inputSeq)
    charInputSeqs <- sapply(inputSeq, function(x){ s2c(x)[1:len_shortest]})
    charGLSeq <- s2c(germlineSeq)
    matClone <- sapply(1:len_shortest, function(i){
        
        # Identify the nucleotides (in seqs and germline) at the current position
        posNucs = unique(charInputSeqs[i,])
        posGL = charGLSeq[i]
        error = FALSE
        
        # If the current position is a gap in both germline and the sequence,
        # return a gap
        if(posGL %in% c("-", ".") & sum(!(posNucs%in%c("-", ".", "N", "n")))==0 ){
            return(c(".",error))
        }
        
        # If all the sequences in the clone have the same nucleotide at the current
        # position, return the value at the current positions
        if(length(posNucs)==1)
            return(c(posNucs[1],error))
        else{         
            if("N"%in%posNucs){
                error=TRUE
            }
            
            # If the current nucleotide matches germline, return germline 
            if(sum(!posNucs[!posNucs%in%c("N", "n")]%in%posGL)==0){
                return( c(posGL,error) )
            }else{
                #return( c(sample(posNucs[posNucs!="N"],1),error) )
                
                # If we look at all nodes (including terminal nodes), sample a nucleotide from the possible
                # nucleotides in the clonal sequences at this position
                if(!nonTerminalOnly){
                    return( c(sample(charInputSeqs[i,!charInputSeqs[i,]%in% c("N", "n") & charInputSeqs[i,]!=posGL],1),error) )
                }else{
                    
                    # If we look at only non-terminal nodes, we only sample the nucleotides that appear more 
                    # than once (this is a quick approximation)
                    posNucs = charInputSeqs[i,!charInputSeqs[i,]%in% c("N", "n") & charInputSeqs[i,]!=posGL]
                    posNucsTable = table(posNucs)
                    if(sum(posNucsTable>1)==0){
                        return( c(posGL,error) )
                    }else{
                        return( c(sample( posNucs[posNucs%in%names(posNucsTable)[posNucsTable>1]],1),error) )
                    }
                }
                
            }
        }
        if(error==TRUE){warning("Error while attempting to collapse by clone!")}
    })
    
    return( c2s(matClone[1,]) )
}



#### Mutation counting functions ####

#' Calculate observed numbers of mutations
#'
#' \code{calcDBObservedMutations} calculates the observed number of mutations for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence.
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. Default is \code{FALSE}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. If NULL, mutations 
#'                               are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    nproc               number of cores to distribute the operation over. If the 
#'                               cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with observed mutation counts for each 
#'           sequence listed. The columns names are dynamically created based on the
#'           regions in the \code{regionDefinition}. For example, when using the
#'           \link{IMGT_V_NO_CDR3} definition, which defines positions for CDR and
#'           FWR, the following columns are added:
#'           \itemize{
#'             \item  \code{OBSERVED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{OBSERVED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{OBSERVED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{OBSERVED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' Mutation count are determined by comparing the input sequences (in the column specified 
#' by \code{sequenceColumn}) to the germline sequence (in the column specified by 
#' \code{germlineColumn}). 
#' 
#' The mutations are binned as either replacement (R) or silent (S) across the different 
#' regions of the sequences as defined by \code{regionDefinition}. Typically, this would 
#' be the framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
#' nucleotide sequences. Mutation counts are appended to the input \code{db} as 
#' additional columns.
#' 
#' @seealso  
#' \link{calcObservedMutations} is called by this function to get the list of mutations 
#' in each sequence grouped by the \link{RegionDefinition}. 
#' See \link{IMGT_SCHEMES} for a set of predefined \link{RegionDefinition} objects.
#' See \link{calcDBExpectedMutations} for calculating expected mutation frequencies.
#'           
#' 
#' @examples
#' # Subset example data
#' db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
#'                           BARCODE %in% c("RL016","RL018","RL019","RL021"))
#'
#' # Run calcDBObservedMutations()
#' db_new <- calcDBObservedMutations(db, sequenceColumn="SEQUENCE_IMGT",
#'                                   germlineColumn="GERMLINE_IMGT_D_MASK",
#'                                   frequency=TRUE,
#'                                   regionDefinition=IMGT_V_NO_CDR3,
#'                                   nproc=1)
#'                      
#' @export
calcDBObservedMutations <- function(db, 
                                    sequenceColumn="SEQUENCE_IMGT",
                                    germlineColumn="GERMLINE_IMGT_D_MASK",
                                    frequency=FALSE,
                                    regionDefinition=NULL,
                                    mutationDefinition=NULL,
                                    nproc=1) {
    # Hack for visibility of data.table and foreach index variables
    idx <- NULL
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition', 'frequency',
                                              'calcObservedMutations','s2c','c2s','NUCLEOTIDES',
                                              'getCodonPos','getContextInCodon','mutationType',
                                              'translateCodonToAminoAcid','AMINO_ACIDS','binMutationsByRegion',
                                              'collapseMatrixToVector'), 
                                envir=environment())
        registerDoParallel(cluster)
    } else if (nproc==1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating observed number of mutations...\n")
    
    # Identify all the mutations in the sequences
    # observedMutations helper function returns a list (1 element per sequence)
    # containing an array of mutations (s or R) and the array labels indicate
    # the nucleotide position of the mutations.
    numbOfSeqs <- nrow(db)
    observedMutations_list <-
        foreach(idx=iterators::icount(numbOfSeqs)) %dopar% {
            calcObservedMutations(db[idx, sequenceColumn], 
                                  db[idx, germlineColumn],
                                  frequency=frequency,
                                  regionDefinition=regionDefinition,
                                  mutationDefinition=mutationDefinition)
        }
    
    # Convert list of mutations to data.frame
    if (!is.null(regionDefinition)) {
        labels_length <- length(regionDefinition@labels)
    } else{
        #labels_length=1
        labels_length <- length(makeNullRegionDefinition()@labels)
    }
    observed_mutations <- do.call( rbind, lapply(observedMutations_list, function(x) { 
        length(x) <- labels_length 
        return(x)
    }))
    
    
    sep <- "_"
    if (ncol(observed_mutations) > 1) sep <- "_"
    observed_mutations[is.na(observed_mutations)] <- 0
    if (frequency == TRUE) {
        colnames(observed_mutations) <- paste("MU_FREQ", colnames(observed_mutations), sep=sep)
    } else {
        colnames(observed_mutations) <- paste("OBSERVED", colnames(observed_mutations), sep=sep)
    }
    
    # Properly shutting down the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, observed_mutations)
    return(db_new)    
}


#' Count the number of observed mutations in a sequence.
#'
#' \code{calcObservedMutations} determines all the mutations in a given input seqeunce compared
#' to its germline sequence.
#'
#' @param    inputSeq            input sequence.
#' @param    germlineSeq         germline sequence.
#' @param    frequency           \code{logical} indicating whether or not to calculate
#'                               mutation frequencies. Default is \code{FALSE}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences. Note, only the part of
#'                               sequences defined in \code{regionDefinition} are analyzed.
#'                               If NULL, mutations are counted for entire sequence.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#'                               
#' @return   An \code{array} of the mutations, replacement (R) or silent(S), with the 
#'           names indicating the nucleotide postion of the mutations in the sequence.
#'           
#' @details
#' Each mutation is considered independently in its codon context. Note, only the part of 
#' \code{inputSeq} defined in \code{regionDefinition} is analyzed. For example, when using 
#' the default \link{IMGT_V_NO_CDR3} definition, then mutations in positions beyond 
#' 312 will be ignored.
#' 
#' @seealso  See \link{calcDBObservedMutations} for counting the number of observed mutations.
#' 
#' @examples
#' # Extracting the first entry in the example data to use for input and germline sequences.
#' inputSeq <- InfluenzaDb[1, "SEQUENCE_IMGT"]
#' germlineSeq <-  InfluenzaDb[1, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' calcObservedMutations(inputSeq, germlineSeq)
#' 
#' # Identify only mutations the V segment minus CDR3
#' calcObservedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3)
#'  
#' # Identify mutations by change in hydropathy class
#' calcObservedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3,
#'                       mutationDefinition=HYDROPATHY_MUTATIONS, frequency=TRUE)
#' 
#' @export
calcObservedMutations <- function(inputSeq, germlineSeq, frequency=FALSE,
                                  regionDefinition=NULL, mutationDefinition=NULL) {
    # Assign mutation definition
    aminoAcidClasses <- if (is.null(mutationDefinition)) { NULL } else { mutationDefinition@classes }
        
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    germlineSeq <- gsub("\\.\\.\\.", "XXX", germlineSeq)
    #If there is a single gap left convert it to an N
    germlineSeq <- gsub("\\.", "N", germlineSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    germlineSeq <- gsub("XXX", "...", germlineSeq)
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    inputSeq <- gsub("\\.\\.\\.", "XXX", inputSeq)
    #If there is a single gap left convert it to an N
    inputSeq <- gsub("\\.", "N", inputSeq)
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    inputSeq <- gsub("XXX", "...", inputSeq)    
    
    # Trim the input and germline sequence to the shortest
    len_inputSeq <- nchar(inputSeq)
    len_germlineSeq <- nchar(germlineSeq)
    
    # If a regionDefinition is passed,
    # then only analyze till the end of the defined length
    if(!is.null(regionDefinition)) {
        rdLength  <- regionDefinition@seqLength
    } else {
        rdLength <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        # Create full sequence RegionDefinition object
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    len_shortest <- min(c(len_inputSeq, len_germlineSeq, rdLength), na.rm=TRUE)
    
    c_inputSeq <- s2c(inputSeq)[1:len_shortest]
    c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
    
    # If the sequence and germline (which now should be the same length) is shorter
    # than the rdLength, pad it with Ns
    if(len_shortest<rdLength){
        fillWithNs <- array("N",rdLength-len_shortest)
        c_inputSeq <- c( c_inputSeq, fillWithNs)
        c_germlineSeq <- c( c_germlineSeq, fillWithNs)
    }
    
    mutations_array <- NA
    mutations = (c_germlineSeq != c_inputSeq) & (c_germlineSeq%in%NUCLEOTIDES[1:5]) & (c_inputSeq%in%NUCLEOTIDES[1:5])
    if (sum(mutations) > 0) {
        # The nucleotide positions of the mutations
        mutations_pos <- which(mutations==TRUE)
        # For every mutations_pos, extract the entire codon from germline
        mutations_pos_codons <- array(sapply(mutations_pos,getCodonPos))
        c_germlineSeq_codons <- c_germlineSeq[mutations_pos_codons]
        # For every mutations_pos, extract the codon from germline (without other mutations 
        # at the same codon, if any).
        c_inputSeq_codons <- array(sapply(mutations_pos, function(x) {
            seqP = c_germlineSeq[getCodonPos(x)]
            seqP[getContextInCodon(x)] = c_inputSeq[x]
            return(seqP) }))
        # split the string of codons into vector of codons
        c_germlineSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_germlineSeq_codons)), " ")[[1]]
        c_inputSeq_codons <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", c2s(c_inputSeq_codons)), " ")[[1]]
        
        # Determine whether the mutations are R or S
        mutations_array <- apply(rbind(c_germlineSeq_codons, c_inputSeq_codons), 2, 
                                 function(x) { mutationType(c2s(x[1]), c2s(x[2]), aminoAcidClasses=aminoAcidClasses) })
        names(mutations_array) = mutations_pos
        mutations_array<- mutations_array[!is.na(mutations_array)]
        if(length(mutations_array)==sum(is.na(mutations_array))){
            mutations_array<-NA    
        } else {
            mutations_array <- binMutationsByRegion(mutations_array,regionDefinition)
            if (frequency==TRUE) {
                tempNames <- sapply(regionDefinition@labels, function(x) { substr(x, 1, nchar(x)-2) })
                # Subset boundaries to only non N bases (in both seq and gl)
                boundaries <- regionDefinition@boundaries[c_inputSeq%in%NUCLEOTIDES[1:4] &  
                                                              c_germlineSeq%in%NUCLEOTIDES[1:4]]
                # Freq = numb of mutations / numb of non N bases (in both seq and gl)
                denoms <- sapply(tempNames, function(x) { sum(boundaries==x) })
                mutations_array <- mutations_array/denoms
            }
        }        
    }    
    return(mutations_array)
}


# Aggregate mutations by region
#
# \code{binMutationsByRegion} takes an array of observed mutations (e.g. from 
# \code{\link{calcObservedMutations}}) and bins them by the different regions defined in the 
# \code{regionDefinition}.
#
# @param   mutationsArray     \code{array} containing the mutations (R/S) with the names
#                             indicating the nucleotide positions of the mutations.                             
# @param   regionDefinition   \code{\link{RegionDefinition}} object defining the regions
#                             and boundaries of the Ig sequences.
# 
# @return An \code{array} of R/S mutations binned across all the unique regions defined
# by \code{regionDefinition}.
# 
# @details
# Note, only the part of sequences defined in \code{regionDefinition} are analyzed.
# For example, if the default \link{IMGT_V_NO_CDR3} definition is used, then mutations
# in positions beyond 312 will be ignored.
# 
# @seealso  
# See \code{\link{calcDBObservedMutations}} for identifying and counting the 
# number of observed mutations.
# This function is also used in \code{\link{calcObservedMutations}}.
# 
# @examples
# # Generate a random mutation array
# numbOfMutations <- sample(3:10, 1) 
# posOfMutations <- sort(sample(330, numbOfMutations))
# mutation_types <- sample(c("R","S"), length(posOfMutations), replace=TRUE)
# mutations_array <- array(mutation_types, dimnames=list(posOfMutations))
# 
# # Random mutations
# binMutationsByRegion(mutations_array, regionDefinition=NULL)
# binMutationsByRegion(mutations_array, regionDefinition=IMGT_V_NO_CDR3)
binMutationsByRegion <- function(mutationsArray, 
                                 regionDefinition=NULL) {
    # Create full sequence RegionDefinition object 
    # The seqLength will be the largest index of a mutation
    if (is.null(regionDefinition)) {
        regionDefinition <- makeNullRegionDefinition(max(as.numeric(names(mutationsArray))))
    }
    
    # Make a factor of R/S
    mutatedPositions <- as.numeric(names(mutationsArray))
    mutations <- array(NA,  dim=regionDefinition@seqLength)
    mutations[mutatedPositions] <- mutationsArray
    mutations <- mutations[1:regionDefinition@seqLength]
    mutations <- factor(mutations,levels=c("R", "S"))
    
    mutations_region_counts <- collapseMatrixToVector(table(regionDefinition@boundaries, mutations))
    
    sortingOrder <- match(regionDefinition@labels, names(mutations_region_counts))
    mutations_region_counts <- mutations_region_counts[sortingOrder]
    return(mutations_region_counts)
}


#### Expected frequencies calculating functions ####

#' Calculate expected mutation frequencies
#'
#' \code{calcDBExpectedMutations} calculates the expected mutation frequencies for each 
#' sequence in the input \code{data.frame}.
#'
#' @param    db                  \code{data.frame} containing sequence data.
#' @param    sequenceColumn      \code{character} name of the column containing input 
#'                               sequences.
#' @param    germlineColumn      \code{character} name of the column containing 
#'                               the germline or reference sequence.
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#' @param    nproc               \code{numeric} number of cores to distribute the operation
#'                               over. If the cluster has already been set the call function with 
#'                               \code{nproc} = 0 to not reset or reinitialize. Default is 
#'                               \code{nproc} = 1.
#' 
#' @return   A modified \code{db} \code{data.frame} with expected mutation frequencies 
#'           for each region defined in \code{regionDefinition}.
#'          
#'           The columns names are dynamically created based on the regions in  
#'           \code{regionDefinition}. For example, when using the \link{IMGT_V_NO_CDR3}
#'           definition, which defines positions for CDR and FWR, the following columns are
#'           added:  
#'           \itemize{
#'             \item  \code{EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                            CDR2 of the V-segment.
#'             \item  \code{EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                            of the V-segment.
#'             \item  \code{EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                            FWR2 and FWR3 of the V-segment.
#'             \item  \code{EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                            FWR3 of the V-segment.
#'           }
#'           
#' @details
#' Only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the \link{IMGT_V_NO_CDR3} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  
#' \link{calcExpectedMutations} is called by this function to calculate the expected 
#' mutation frequencies. See \link{calcDBObservedMutations} for getting observed 
#' mutation counts. See \link{IMGT_SCHEMES} for a set of predefined 
#' \link{RegionDefinition} objects.
#' 
#' @examples
#' # Subset example data
#' db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
#'                           BARCODE %in% c("RL016","RL018","RL019","RL021"))
#'
#' # Calculate expected mutations over V region
#' db <- calcDBExpectedMutations(db,
#'                               sequenceColumn="SEQUENCE_IMGT",
#'                               germlineColumn="GERMLINE_IMGT_D_MASK",
#'                               regionDefinition=IMGT_V_NO_CDR3,
#'                               nproc=1)
#' 
#' # Calculate hydropathy expected mutations over V region
#' db <- calcDBExpectedMutations(db,
#'                               sequenceColumn="SEQUENCE_IMGT",
#'                               germlineColumn="GERMLINE_IMGT_D_MASK",
#'                               regionDefinition=IMGT_V_NO_CDR3,
#'                               mutationDefinition=HYDROPATHY_MUTATIONS,
#'                               nproc=1)
#'
#' @export
calcDBExpectedMutations <- function(db, 
                                    sequenceColumn="SEQUENCE_IMGT",
                                    germlineColumn="GERMLINE_IMGT_D_MASK",
                                    targetingModel=HS5FModel,
                                    regionDefinition=NULL,
                                    mutationDefinition=NULL,
                                    nproc=1) {
    # Hack for visibility of data.table and foreach index variables
    idx <- NULL
    
    # Check for valid columns
    check <- checkColumns(db, c(sequenceColumn, germlineColumn))
    if (check != TRUE) { stop(check) }
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){ 
        cluster = nproc 
        nproc = 0
    }
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc(), na.rm=T)
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if (nproc > 1) {        
        cluster <- parallel::makeCluster(nproc, type = "PSOCK")
        parallel::clusterExport(cluster, list('db', 'sequenceColumn', 'germlineColumn', 
                                              'regionDefinition','targetingModel',
                                              'calcExpectedMutations','calculateTargeting',
                                              's2c','c2s','NUCLEOTIDES','HS5FModel',
                                              'calculateMutationalPaths','CODON_TABLE'),
                                envir=environment() )
        registerDoParallel(cluster)
    } else if (nproc == 1) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    
    # Printing status to console
    cat("Calculating the expected frequencies of mutations...\n")
    
    # Calculate targeting for each sequence (based on the germline)
    # Should be a 5 x N matrix where N in the number of nucleotides defined by
    # the regionDefinition
    numbOfSeqs <- nrow(db)
    
    targeting_list <-
        foreach (idx=iterators::icount(numbOfSeqs)) %dopar% {
            calcExpectedMutations(germlineSeq=db[idx, germlineColumn],
                                  inputSeq=db[idx, sequenceColumn],
                                  targetingModel=targetingModel,
                                  regionDefinition=regionDefinition,
                                  mutationDefinition=mutationDefinition)
        }
    
    # Convert list of expected mutation freq to data.frame
    if (is.null(regionDefinition)) {
        labels_length <- length(makeNullRegionDefinition()@labels)
    } else {
        labels_length <- length(regionDefinition@labels)
    }
    expectedMutationFrequencies <- do.call(rbind, lapply(targeting_list, function(x) { 
        length(x) <- labels_length 
        return(x) })) 
    
    expectedMutationFrequencies[is.na(expectedMutationFrequencies)] <- 0
    colnames(expectedMutationFrequencies) <- paste0("EXPECTED_", colnames(expectedMutationFrequencies))
    
    # Properly shutting down the cluster
    if(nproc>1){ parallel::stopCluster(cluster) }
    
    # Bind the observed mutations to db
    db_new <- cbind(db, expectedMutationFrequencies)
    return(db_new)    
    
}


#' Calculate expected mutation frequencies of a sequence
#'
#' \code{calcExpectedMutations} calculates the expected mutation
#' frequencies of a given sequence. This is primarily a helper function for
#' \link{calcDBExpectedMutations}. 
#'
#' @param    germlineSeq         germline (reference) sequence.
#' @param    inputSeq            input (observed) sequence. If this is not \code{NULL}, 
#'                               then \code{germlineSeq} will be processed to be the same
#'                               same length as \code{inputSeq} and positions in 
#'                               \code{germlineSeq} corresponding to positions with Ns in 
#'                               \code{inputSeq} will also be assigned an N. 
#' @param    targetingModel      \link{TargetingModel} object. Default is \link{HS5FModel}.
#' @param    regionDefinition    \link{RegionDefinition} object defining the regions
#'                               and boundaries of the Ig sequences.
#' @param    mutationDefinition  \link{MutationDefinition} object defining replacement
#'                               and silent mutation criteria. If \code{NULL} then 
#'                               replacement and silent are determined by exact 
#'                               amino acid identity.
#'                               
#' @return   A \code{numeric} vector of the expected frequencies of mutations in the 
#'           regions in the \code{regionDefinition}. For example, when using the default 
#'           \link{IMGT_V_NO_CDR3} definition, which defines positions for CDR and 
#'           FWR, the following columns are calculated:
#'           \itemize{
#'              \item  \code{EXPECTED_CDR_R}:  number of replacement mutations in CDR1 and 
#'                                             CDR2 of the V-segment.
#'              \item  \code{EXPECTED_CDR_S}:  number of silent mutations in CDR1 and CDR2 
#'                                             of the V-segment.
#'              \item  \code{EXPECTED_FWR_R}:  number of replacement mutations in FWR1, 
#'                                             FWR2 and FWR3 of the V-segment.
#'              \item  \code{EXPECTED_FWR_S}:  number of silent mutations in FWR1, FWR2 and
#'                                             FWR3 of the V-segment.
#'            }
#'           
#' @details
#' \code{calcExpectedMutations} calculates the expected mutation frequencies of a 
#' given sequence and its germline. 
#' 
#' Note, only the part of the sequences defined in \code{regionDefinition} are analyzed. 
#' For example, when using the default \link{IMGT_V_NO_CDR3} definition, mutations in
#' positions beyond 312 will be ignored.
#' 
#' @seealso  \link{calcDBExpectedMutations} calls this function.
#' To create a custom \code{targetingModel} see \link{createTargetingModel}.
#' See \link{calcObservedMutations} for getting observed mutation counts.
#' 
#' @examples
#' # Extracting the first entry in the exampled data to use for input and germline sequences.
#' inputSeq <- InfluenzaDb[1, "SEQUENCE_IMGT"]
#' germlineSeq <-  InfluenzaDb[1, "GERMLINE_IMGT_D_MASK"]
#' 
#' # Identify all mutations in the sequence
#' calcExpectedMutations(inputSeq, germlineSeq)
#' 
#' # Identify only mutations the V segment minus CDR3
#' calcExpectedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3)
#' 
#' # Define mutations based on hydropathy
#' calcExpectedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3,
#'                       mutationDefinition=HYDROPATHY_MUTATIONS)
#' 
#' @export
calcExpectedMutations <- function(germlineSeq,
                                  inputSeq=NULL,
                                  targetingModel=HS5FModel,
                                  regionDefinition=NULL,
                                  mutationDefinition=NULL) {
    # Assign codon table
    codonTable <- if (is.null(mutationDefinition)) { CODON_TABLE } else { mutationDefinition@codonTable }
    
    # Get targeting
    targeting <- calculateTargeting(germlineSeq=germlineSeq, 
                                    inputSeq=inputSeq,
                                    targetingModel=targetingModel,
                                    regionDefinition=regionDefinition)
    
    # Determine the mutations paths (i.e. determine R and S mutation frequencies)
    mutationalPaths <- calculateMutationalPaths(germlineSeq=c2s(colnames(targeting)), 
                                                regionDefinition=regionDefinition,
                                                codonTable=codonTable)
    
    typesOfMutations <- c("R", "S")
    mutationalPaths[!(mutationalPaths %in% typesOfMutations)] <- NA
    
    if (is.null(regionDefinition)) {
        rdLength <- max(nchar(inputSeq), nchar(germlineSeq), na.rm=TRUE)
        regionDefinition <- makeNullRegionDefinition(rdLength)
    }
    listExpectedMutationFrequencies <- list()
    for(region in regionDefinition@regions){
        for(typeOfMutation in typesOfMutations){
            region_mutation <- paste(region, typeOfMutation, sep="_")    
            
            targeting_region <- targeting[1:4, regionDefinition@boundaries %in% region]
            mutationalPaths_region <- mutationalPaths[, regionDefinition@boundaries[1:ncol(mutationalPaths)] %in% region]
            targeting_typeOfMutation_region <- sum(targeting_region[mutationalPaths_region == typeOfMutation], 
                                                   na.rm=TRUE)
            
            listExpectedMutationFrequencies[[region_mutation]] <- targeting_typeOfMutation_region
            
        }
    }
    expectedMutationFrequencies <- unlist(listExpectedMutationFrequencies)
    expectedMutationFrequencies[!is.finite(expectedMutationFrequencies)] <- NA
    expectedMutationFrequencies <- expectedMutationFrequencies/sum(expectedMutationFrequencies, na.rm=TRUE)
    return(expectedMutationFrequencies)    
}


calculateTargeting <- function(germlineSeq,
                               inputSeq=NULL,
                               targetingModel=HS5FModel,
                               regionDefinition=NULL) {
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same legth, mask germlineSequence with Ns where inputSequence is also N
    # If not needed then  you may skip this step by passing in inputSequence=NULL 
    # (which is default). 
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- nchar(inputSeq)
        len_germlineSeq <- nchar(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest<length_regionDefinition){
            fillWithNs <- array("N",length_regionDefinition-len_shortest)
            c_inputSeq <- c( c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[ c_inputSeq=="N" |  !c_inputSeq%in%c(NUCLEOTIDES[1:5],".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    }else{
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    
    # Removing IMGT gaps (they should come in threes)
    # After converting ... to XXX any other . is not an IMGT gap & will be treated like N
    gaplessSeq <- gsub("\\.\\.\\.", "XXX", s_germlineSeq)
    #If there is a single gap left convert it to an N
    gaplessSeq <- gsub("\\.", "N", gaplessSeq)
    
    # Re-assigning s_germlineSeq (now has all "." that are not IMGT gaps converted to Ns)
    s_germlineSeq <- gsub("XXX", "...", gaplessSeq)
    c_germlineSeq <- s2c(s_germlineSeq)
    # Matrix to hold targeting values for each position in c_germlineSeq
    germlineSeqTargeting <- matrix(NA, 
                                   ncol=nchar(s_germlineSeq), 
                                   nrow=length(NUCLEOTIDES[1:5]),
                                   dimnames=list(NUCLEOTIDES[1:5], c_germlineSeq))
    
    # Now remove the IMGT gaps so that the correct 5mers can be made to calculate
    # targeting. e.g.
    # GAGAAA......TAG yields: "GAGAA" "AGAAA" "GAAAT" "AAATA" "AATAG"
    # (because the IMGT gaps are NOT real gaps in sequence!!!)
    gaplessSeq <- gsub("\\.\\.\\.", "", s_germlineSeq)
    gaplessSeqLen <- nchar(gaplessSeq)
    
    #Slide through 5-mers and look up targeting
    gaplessSeq <- paste("NN", gaplessSeq, "NN", sep="")
    gaplessSeqLen <- nchar(gaplessSeq)
    pos<- 3:(gaplessSeqLen - 2)
    subSeq =  substr(rep(gaplessSeq, gaplessSeqLen - 4), (pos - 2), (pos + 2))
    germlineSeqTargeting_gapless <- targetingModel@targeting[,subSeq]
#     germlineSeqTargeting_gapless <- sapply(subSeq, function(x) { 
#         targetingModel@targeting[, x] })
    
    germlineSeqTargeting[, c_germlineSeq != "."] <- germlineSeqTargeting_gapless  
    
    # Set self-mutating targeting values to be NA
    mutatingToSelf <- colnames(germlineSeqTargeting)
    mutatingToSelf[!(mutatingToSelf %in% NUCLEOTIDES[1:5])] <- "N"
#     # TODO: What's with this <<- business?
#     # TODO: I think this is assigning NA to all self-mutations, which are already NA
#     sapply(1:ncol(germlineSeqTargeting), function(pos) { germlineSeqTargeting[mutatingToSelf[pos], pos] <<- NA })
    
    germlineSeqTargeting[!is.finite(germlineSeqTargeting)] <- NA
    return(germlineSeqTargeting)
}

calculateMutationalPaths <- function(germlineSeq,
                                     inputSeq=NULL,
                                     regionDefinition=NULL,
                                     codonTable=NULL) {    
    # Set codon table if required
    if (is.null(codonTable)) { codonTable <- CODON_TABLE }
    
    # If an inputSequence is passed then process the germlineSequence
    # to be the same length, mask germlineSequence with Ns where inputSequence is also N
    # If this function is being called after running calculateTargeting you may skip
    # this step by passing in inputSequence=NULL (which is default). This way you save
    # some processing time.
    if(!is.null(inputSeq)){    
        # Trim the input and germline sequence to the shortest
        len_inputSeq <- nchar(inputSeq)
        len_germlineSeq <- nchar(germlineSeq)
        # If a regionDefinition is passed,
        # then only analyze till the end of the defined length
        if(!is.null(regionDefinition)){
            length_regionDefinition  <- regionDefinition@seqLength
        } else{
            length_regionDefinition <- max(len_inputSeq, len_germlineSeq, na.rm=TRUE)
        }
        len_shortest <- min( c(len_inputSeq,len_germlineSeq,length_regionDefinition),  na.rm=TRUE)
        
        c_inputSeq <- s2c(inputSeq)[1:len_shortest]
        c_germlineSeq <- s2c(germlineSeq)[1:len_shortest]
        
        # If the sequence and germline (which now should be the same length) is shorter
        # than the length_regionDefinition, pad it with Ns
        if(len_shortest<length_regionDefinition){
            fillWithNs <- array("N",length_regionDefinition-len_shortest)
            c_inputSeq <- c( c_inputSeq, fillWithNs)
            c_germlineSeq <- c( c_germlineSeq, fillWithNs)
        }
        
        # Mask germline with Ns where input sequence has Ns
        c_germlineSeq[c_inputSeq=="N" |  !c_inputSeq %in% c(NUCLEOTIDES[1:5], ".") ] = "N"    
        s_germlineSeq <- c2s(c_germlineSeq)
    } else {
        s_germlineSeq <- germlineSeq
        c_germlineSeq <- s2c(s_germlineSeq)
    }
    

    s_germlineSeq_len <- nchar(s_germlineSeq)    
    vecCodons <- sapply({1:(s_germlineSeq_len/3)}*3 - 2, function(x) { substr(s_germlineSeq, x, x + 2) })
    vecCodons[!vecCodons %in% colnames(codonTable)] <- "NNN"
    matMutationTypes = matrix(codonTable[, vecCodons], nrow=4, byrow=F,
                              dimnames=list(NUCLEOTIDES[1:4], c_germlineSeq[ 1:(3 * length(vecCodons))]))
    
    return(matMutationTypes)
}

#### Additional helper functions ####
# Given a nuclotide position, returns the codon number
# e.g. nuc 86  = codon 29
getCodonNumb <- function(nucPos){
  return( ceiling(nucPos/3) )
}

# Given a codon, returns all the nuc positions that make the codon
getCodonNucs <- function(codonNumb){
  getCodonPos(codonNumb*3)
}

# Given a nucleotide postions return the position in the codon
getContextInCodon <- function(nucPos){
  return( {nucPos-1}%%3+1 )
}

# Given a nuclotide position, returns the pos of the 3 nucs that made the codon
# e.g. nuc 86 is part of nucs 85,86,87
getCodonPos <- function(nucPos) {
  codonNum =  (ceiling(nucPos / 3)) * 3
  return ((codonNum - 2):codonNum)
}

# Translate codon to amino acid
translateCodonToAminoAcid <- function(Codon) {
  return (AMINO_ACIDS[Codon])
}

# Given two codons, tells you if the mutation is R or S (based on your definition)
#
# @param   codonFrom         starting codon
# @param   codonTo           ending codon
# @param   aminoAcidClasses  vector of amino acid trait classes
#                            if NULL then R or S is determined by amino acid identity
# @return  Mutation type as "R" (replacement), "S" (silent), "Stop" (stop) or NA (input is NA).
#
# @examples
# # Without classes
# shazam:::mutationType("TTT", "TTC")
# shazam:::mutationType("TTT", "TTA")
# shazam:::mutationType("TTT", "TGA")
#
# # With classes
# classes <- HYDROPATHY_MUTATIONS@classes
# shazam:::mutationType("TTT", "TTC", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TTA", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TCT", aminoAcidClasses=classes)
# shazam:::mutationType("TTT", "TGA", aminoAcidClasses=classes)
mutationType <- function(codonFrom, codonTo, aminoAcidClasses=NULL) {
    # codonFrom="TTT"; codonTo="TTA"
    # codonFrom="TTT"; codonTo="TGA"
    
    # Translate codons
    aaFrom <- translateCodonToAminoAcid(codonFrom)
    aaTo <- translateCodonToAminoAcid(codonTo)
    
    # If any codon is NA then return NA
    if (any(is.na(c(codonFrom, codonTo, aaFrom, aaTo)))) { 
        return(NA) 
    }
    
    # If any amino acid is Stop then return "Stop"
    if (any(c(aaFrom, aaTo) == "*")) { 
        return("Stop") 
    }
    
    if (is.null(aminoAcidClasses)) {
        # Check for exact identity if no amino acid classes are specified
        mutation <- if (aaFrom == aaTo) { "S" } else { "R" }
    } else {
        # Check for amino acid class identity if classes are specified
        mutation <- if (aminoAcidClasses[aaFrom] == aminoAcidClasses[aaTo]) { "S" } else { "R" }
    }
    return(mutation)
}