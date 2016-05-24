# Generates distance to nearest neighbor

#' @include shazam.R
NULL


# Returns a 5-mer sliding window of given sequence
#
# @param   strSequence   The sequence string
# @return  An array of 5-mer sliding windows
#
# @examples
# slidingArrayOf5mers("ACGTNACGTNACGTN")
slidingArrayOf5mers <- function(strSequence){
    seqLength <- nchar(strSequence)
    return( substr( rep(strSequence,seqLength-4), 1:(seqLength-4), 5:seqLength ) )
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1                first nucleotide sequence, broken into 5mers.
# @param    seq2                second nucleotide sequence, broken into 5mers.
# @param    targeting_model     targeting model.
# @param    normalize           The method of normalization. Default is "none".
#                               "length" = normalize distance by length of junction.
# @param    symmetry            if model is hs5f, distance between seq1 and seq2 is either the
#                               average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @return   distance between two sequences.
#
# @examples
# seq1 = c("NNACG", "NACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
# seq2 = c("NNACG", "NACGA", "ACGAA", "CGAAC", "GAACG", "AACGT", "ACGTA", "CGTAC", "GTACG", "TACGT", "ACGTN", "CGTNN")
#
# distSeq5mers(seq1, seq2, HS5FModel)
distSeq5mers <- function(seq1, seq2, targeting_model, 
                         normalize=c("none" ,"length", "mutations"),
                         symmetry=c("avg","min")) {
  # Evaluate choices
  normalize <- match.arg(normalize)
  symmetry <- match.arg(symmetry)
  
  # Get distance from targeting model
  targeting_dist <- calcTargetingDistance(targeting_model)
  
  # Compute length of sequence (for normalization, if specified)
  juncLength <- length(seq1)

  # Compute distance only on fivemers that have mutations
  fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
  #fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
  #fivemersWithMu <- fivemersWithMu & fivemersWithNonNuc
  seq1 <- seq1[fivemersWithMu]
  seq2 <- seq2[fivemersWithMu]
  
  # Number of mutations (for normalization, if specified)
  numbOfMutation <- sum(fivemersWithMu)

  dist <- NA
  tryCatch({
    if (length(seq1)==1){
      seq1_to_seq2 <- targeting_dist[substr(seq2,3,3),seq1]
      seq2_to_seq1 <- targeting_dist[substr(seq1,3,3),seq2]
    } else {
      seq1_to_seq2 <- sum( diag(targeting_dist[substr(seq2,3,3),seq1]) )
      seq2_to_seq1 <- sum( diag(targeting_dist[substr(seq1,3,3),seq2]) )
    }
    if (symmetry == "avg") {
      dist <- mean(c(seq1_to_seq2, seq2_to_seq1))
    } else if (symmetry == "min") {
      dist <- min(c(seq1_to_seq2, seq2_to_seq1))
    }
  },error = function(e){
    warning("Invalid sequence. Cannot compute distance.")
    return(NA)
  })

  # Normalize distances
  if (normalize == "length") { 
      dist <- dist/juncLength
  } else if (normalize == "mutations") { 
      dist <- dist/numbOfMutation 
  }
  
  return(dist)
}


# Get distance between two sequences of same length, broken by a sliding window of 5mers
#
# @param    seq1          first nucleotide sequence.
# @param    seq2          second nucleotide sequence.
# @param    model         DNA (ham) or amino acid (aa) hamming distance model or
#                         mouse (m1n) or human (hs1f) single nucleotide distance model
# @param    normalize     The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   distance between two sequences.
#
# @examples
# seq1 = "ATG-C"
# seq2 = "AT--C"
# 
# distSeqMat(seq1, seq2)
distSeqMat <- function(seq1, seq2, model=c("ham","aa","m1n","hs1f"),
                       normalize=c("none" ,"length", "mutations")) {
  # Evaluate choices
  model <- match.arg(model)
  normalize <- match.arg(normalize)
  
  # Get character distance matrix
  if (model == "ham") {
    dist_mat <- getDNAMatrix(gap=0)
  } else if (model == "m1n") {
    dist_mat <- M1NDistance
  } else if (model == "hs1f") {
    dist_mat <- HS1FDistance
  } else if (model == "aa") {
    
    # Translate sequences
    seq1 <- strsplit(tolower(gsub("[-.]","N",seq1)), "")[[1]]
    seq2 <- strsplit(tolower(gsub("[-.]","N",seq2)), "")[[1]]
    seq1 <- translate(seq1, ambiguous=T)
    seq2 <- translate(seq2, ambiguous=T)
    
    dist_mat <- getAAMatrix()
  }
  
  # Calculate distance
  dist <- tryCatch(getSeqDistance(seq1, seq2, dist_mat=dist_mat),
                   error=function(e) {
                       warning("Invalid character in sequence. Cannot compute distance.")
                       return(NA)
                   })
  
  # Normalize distances
  if (normalize == "length") { 
    dist <- dist/sum(nchar(seq1))
  } else if (normalize == "mutations") {
    dist <- dist/sum(strsplit(seq1,"")[[1]] != strsplit(seq2,"")[[1]])
  }
  
  return(dist)
}


# Given an array of nucleotide sequences, find the pairwise distances
# 
# @param   arrJunctions   character vector of nucleotide sequences.
# @param   model          SHM targeting model.
# @param   normalize      The method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
# @param    symmetry      if model is hs5f, distance between seq1 and seq2 is either the
#                         average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @return  A matrix of pairwise distances between junction sequences.
# 
# @details
# needs method details
# 
# @seealso needs links
# 
# @examples
# # working example
getPairwiseDistances <- function(arrJunctions, targeting_model, 
                                 normalize=c("none" ,"length", "mutations"),
                                 symmetry=c("avg","min")) {
  # Initial checks
  normalize <- match.arg(normalize)
  symmetry <- match.arg(symmetry)
  
  # Convert junctions to uppercase
  arrJunctions <- toupper(arrJunctions)
  # Convert gaps to Ns
  arrJunctions <- gsub('[-.]', 'N', arrJunctions, fixed=T)
  # Add 'NN' to front and end of each sequence for fivemers
  arrJunctions <- as.vector(sapply(arrJunctions, function(x){ paste("NN", x, "NN", sep="") }))

  numbOfJunctions<-length(arrJunctions)

  #Junctions are broken in to 5-mers based on a sliding window (of one) and placed in matrix
  #Each column is a junction
  #E.g. junctions 1234567, ABCDEFG, JKLMNOP becomes:
  # 12345   ABCDE   JKLMN
  # 23456   BCDEF   KLMNO
  # 34567   CDEFG   LMNOP
  .matSeqSlidingFiveMer <- sapply(arrJunctions, function(x) { slidingArrayOf5mers(x) }, simplify="matrix")

  # Compute pairwise distance between all sequences' fivemers (by column)
  matDistance <-
    sapply(1:numbOfJunctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfJunctions, function(j) {
      distSeq5mers(.matSeqSlidingFiveMer[,i],
                   .matSeqSlidingFiveMer[,j],
                   targeting_model,
                   normalize=normalize,
                   symmetry=symmetry)
    })))
  # Make distance matrix symmetric
  matDistance <- matDistance + t(matDistance)
  return(matDistance)
}


# Given an array of junctions, generate distance array for pairwise distances with
# 0 if junction is non-unique and return this array along with unique junctions that 
# are only observed once
findUniqueJunctions <- function(arrJunctions) {
  # Initialize array of distances
  arrJunctionsDist <- rep(NA,length(arrJunctions))
  
  # Filter unique junctions
  arrJunctionsUnique <- unique(arrJunctions)
  
  # Map indices of unique to its non-unique in the original arrJunctions
  indexJunctions <- match(arrJunctions, arrJunctionsUnique)
  
  # Identify junctions with multiple non-unique sequences and set its distances to 0
  indexJunctionsCounts <- table(indexJunctions)
  indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
  indexRepeated <- indexJunctions%in%indexRepeated
  # arrJunctionsDist[ indexRepeated ] <- 0
  names(arrJunctionsDist) <- arrJunctions
  
  # Subset unique junctions to those that are only observed once
  #arrJunctionsUnique <- arrJunctionsUnique[indexJunctionsCounts==1]
  
  return(list('arrJunctionsDist'=arrJunctionsDist, 
              'arrJunctionsUnique'=arrJunctionsUnique,
              'indexJunctionsCounts'=indexJunctionsCounts))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    targeting_model     targeting model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @param    symmetry      if model is hs5f, distance between seq1 and seq2 is either the
#                         average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
# @return   A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getClosestBy5mers(arrJunctions, HS5FModel, normalize="none" )
getClosestBy5mers <- function(arrJunctions, targeting_model, 
                              normalize=c("none" ,"length", "mutations"),
                              symmetry=c("avg","min")) {
  # Initial checks
  normalize <- match.arg(normalize)
  symmetry <- match.arg(symmetry)
  
  # Find unique sequences and return distance array with mapping
  l <- findUniqueJunctions(arrJunctions)
  arrJunctionsDist <- l$arrJunctionsDist
  arrJunctionsUnique <- l$arrJunctionsUnique
  # indexJunctionsCounts <- l$indexJunctionsCounts

  # Compute distances between junctions
  numbOfUniqueJunctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
  if (numbOfUniqueJunctions>1){
    # Calculate symmetric distance matrix
    matDistance <- getPairwiseDistances(arrJunctionsUnique, targeting_model, normalize, symmetry)
    # Find minimum distance for each sequence
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i){ min(matDistance[,i][matDistance[,i]>0]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }

  # Fill the distances for unique sequences
  # arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist[indexJunctionsCounts==1]
  arrJunctionsDist <- arrUniqueJunctionsDist[match(names(arrJunctionsDist), names(arrUniqueJunctionsDist))]
  return(round(arrJunctionsDist,4))
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param    arrJunctions  character vector of junction sequences.
# @param    model         DNA (ham) or amino acid (aa) hamming distance model or
#                         mouse (m1n) or human (hs1f) single nucleotide distance model
# @param    normalize     method of normalization. Default is "none".
#                         "length" = normalize distance by length of junction.
#                         "mutations" = normalize distance by number of mutations in junction.
# @return   A vector of distances to the closest sequence.
# @examples
# arrJunctions <- c( "ACGTACGTACGT","ACGAACGTACGT",
#                    "ACGAACGTATGT", "ACGAACGTATGC",
#                    "ACGAACGTATCC","AAAAAAAAAAAA")
# getClosestMat(arrJunctions, normalize="none" )
getClosestMat <- function(arrJunctions, model=c("ham","aa","m1n","hs1f"),
                          normalize=c("none" ,"length", "mutations")) {
  # Initial checks
  model <- match.arg(model)
  normalize <- match.arg(normalize)
  
  # Find unique sequences and return distance array with mapping
  l <- findUniqueJunctions(arrJunctions)
  arrJunctionsDist <- l$arrJunctionsDist
  arrJunctionsUnique <- l$arrJunctionsUnique
  # indexJunctionsCounts <- l$indexJunctionsCounts
  
  # Compute distances between junctions
  numbOfUniqueJunctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJunctions)
  if (numbOfUniqueJunctions>1){
    # Calculate symmetric distance matrix
    matDistance <-
      sapply(1:numbOfUniqueJunctions, function(i) c(rep.int(0,i-1), sapply(i:numbOfUniqueJunctions, function(j) {
        distSeqMat(arrJunctionsUnique[i], arrJunctionsUnique[j], model=model, normalize=normalize)
      })))
    matDistance <- matDistance + t(matDistance)
    # Find minimum distance for each sequence
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJunctions, function(i){ min(matDistance[,i][matDistance[,i]>0]) })
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }
  
  # Fill the distances for unique sequences
  # arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist[indexJunctionsCounts==1]
  arrJunctionsDist <- arrUniqueJunctionsDist[match(names(arrJunctionsDist), names(arrUniqueJunctionsDist))]
  return(round(arrJunctionsDist,4))
}

#' Distance to nearest neighbor
#'
#' Get distance of every sequence to its nearest sequence sharing same V gene, J gene, and
#' sequence length.
#'
#' @param    db              data.frame containing sequence data.
#' @param    sequenceColumn  name of the column containing nucleotide sequences to compare. 
#'                           Also used to determine sequence length for grouping.
#' @param    vCallColumn     name of the column containing the V-segment allele calls.
#' @param    jCallColumn     name of the column containing the J-segment allele calls.
#' @param    model           underlying SHM model, which must be one of 
#'                           \code{c("m1n", "ham", "aa", "hs5f")}.
#'                           See Details for further information.
#' @param    normalize       method of normalization. The default is \code{"length"}, which 
#'                           divides the distance by the length of the sequence group. If 
#'                           \code{"none"} then no normalization if performed
#' @param    symmetry        if model is hs5f, distance between seq1 and seq2 is either the
#'                           average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).
#' @param    first           if \code{TRUE} only the first call of the gene assignments is used.
#'                           If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                           group all sequences with any overlapping gene calls.
#' @param    nproc           number of cores to distribute the function over.
#' @param    fields          Additional fields to use for grouping
#'
#' @return   Returns a modified \code{db} data.frame with nearest neighbor distances in the 
#'           \code{DIST_NEAREST} column.
#'
#' @details
#' The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
#' sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
#' with the ideal threshold being a value that separates the two modes.
#' 
#' "hs5f" use distance derived from the \link{HS5FModel}
#' using \link{calcTargetingDistance}. "hs1f" and "m1n" use \link{HS1FDistance} and \link{M1NDistance}
#'  to calculate distances respectively. "ham" uses a nucleotide hamming distance matrix from 
#'  \link{getDNAMatrix}, with gaps being zero. "aa" uses an amino acid hamming distance matrix 
#'  from \link{getAAMatrix}.
#' 
#' @references
#' \enumerate{
#'   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#'            mutagenesis in normal and autoreactive B cells. 
#'            J Immunol. 1996 156:2642-52. 
#'   \item  Glanville J, Kuo TC, von Budingen H-C, et al. 
#'            Naive antibody gene-segment frequencies are heritable and unaltered by 
#'            chronic lymphocyte ablation. 
#'            Proc Natl Acad Sci USA. 2011 108(50):20066-71.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'  }
#'  
#' @seealso  See \link{calcTargetingDistance} for generating nucleotide distance matrices 
#'           from a \link{TargetingModel} object. See \link{M1NDistance}, 
#'           \link{HS5FModel}, \link{getDNAMatrix}, and \link{getAAMatrix}
#'           for individual model details.
#' 
#' @examples
#' # Subset data for demo purposes
#' db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
#'              BARCODE %in% c("RL016","RL018","RL019","RL021"))
#' 
#' # Use genotyped V assignments, HS1F model, and normalize by junction length
#' dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
#'                            model="hs1f", first=FALSE, normalize="length")
#'                            
#' # Plot histogram of non-NA distances
#' p1 <- ggplot(data=subset(dist_hs1f, !is.na(DIST_NEAREST))) + theme_bw() + 
#'     ggtitle("Distance to nearest: hs1f") + xlab("distance") +
#'     geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
#'                    fill="steelblue", color="white")
#' plot(p1)
#'
#' @export
distToNearest <- function(db, sequenceColumn="JUNCTION", vCallColumn="V_CALL", 
                          jCallColumn="J_CALL", model=c("hs1f", "m1n", "ham", "aa", "hs5f"), 
                          normalize=c("length", "none"), symmetry=c("avg","min"),
                          first=TRUE, nproc=1, fields=NULL) {
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL
    
    # Initial checks
    model <- match.arg(model)
    normalize <- match.arg(normalize)
    symmetry <- match.arg(symmetry)
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    
    # Check for valid columns
    columns <- c(sequenceColumn, vCallColumn, jCallColumn,fields)
    columns <- columns[!is.null(columns)]
    
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Convert case check for invalid characters
    db[, sequenceColumn] <- toupper(db[[sequenceColumn]])
    #check <- grepl("[^ACGTN]", db[[sequenceColumn]], perl=TRUE)
    #if (any(check)) {
    #  stop("Invalid sequence characters in the ", sequenceColumn, " column.")
    #}
    
    # Get targeting model
    if (model == "hs5f") {
        targeting_model <- HS5FModel
    }
    
    # Parse V and J columns to get gene
    # cat("V+J Column parsing\n")
    if (first) {
        db$V <- getGene(db[[vCallColumn]])
        db$J <- getGene(db[[jCallColumn]])
    } else {
        db$V1 <- getGene(db[[vCallColumn]], first=FALSE)
        db$J1 <- getGene(db[[jCallColumn]], first=FALSE)
        db$V <- db$V1
        db$J <- db$J1
        # Reassign V genes to most general group of genes
        for(ambig in unique(db$V1[grepl(',', db$V1)])) {
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$V[grepl(g, db$V1)] = ambig
            }
        }
        # Reassign J genes to most general group of genes
        for(ambig in unique(db$J1[grepl(',',db$J1)])) {
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$J[grepl(g, db$J1)] = ambig
            }
        }
    }
    
    # Create new column for distance to nearest neighbor
    db$DIST_NEAREST <- rep(NA, nrow(db))
    db$ROW_ID <- 1:nrow(db)
    db$L <- nchar(db[[sequenceColumn]])
    
    # Create cluster of nproc size and export namespaces
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # Calculate distance to nearest neighbor
    # cat("Calculating distance to nearest neighbor\n")
    
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of V J L, instead of doing dplyr
    dt <- data.table(db)
    # Get the group indexes
    group_cols <- c("V","J","L")
    if (!is.null(fields)) {
        group_cols <- append(group_cols,fields)
    }
    dt <- dt[, list( yidx = list(.I) ) , by = group_cols ]
    groups <- dt[,yidx]
    lenGroups <- length(groups)
    
    # Export groups to the clusters
    if (nproc>1) { parallel::clusterExport(cluster, list("db", 
                                               "groups", 
                                               "sequenceColumn","model",
                                               "normalize","symmetry","getClosestMat", 
                                               "HS1FDistance","distSeqMat",
                                               "calcTargetingDistance",
                                               "findUniqueJunctions","getPairwiseDistances"), envir=environment()) }
    
    if (model %in% c("hs5f")) {
        # Export targeting model to processes
        if (nproc>1) { parallel::clusterExport(cluster, list("targeting_model","getClosestBy5mers"), envir=environment()) }    
        list_db <-
            foreach(idx=iterators::icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[idx]],]
                db_group$DIST_NEAREST <-
                    getClosestBy5mers( db[groups[[idx]],sequenceColumn],
                                       targeting_model=targeting_model,
                                       normalize=normalize,
                                       symmetry=symmetry )
                return(db_group)
            }    
    } else if (model %in% c("ham", "aa", "m1n", "hs1f")) {    
        list_db <-
            foreach(idx=iterators::icount(lenGroups), .errorhandling='pass') %dopar% {
                db_group <- db[groups[[idx]],]
                db_group$DIST_NEAREST <-
                    getClosestMat( db[groups[[idx]],sequenceColumn],
                                   model=model,
                                   normalize=normalize )
                return(db_group)
            }        
    }
    
    # Convert list from foreach into a db data.frame
    db <- dplyr::bind_rows(list_db)
    db <- db[order(db$ROW_ID),]
    
    # Stop the cluster
    if( nproc>1) { parallel::stopCluster(cluster) }
    
    return(db[, !(names(db) %in% c("V", "J", "L", "ROW_ID", "V1", "J1"))])
}
