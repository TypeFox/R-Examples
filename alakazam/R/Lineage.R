# Ig lineage reconstruction via maximum parsimony

#' @include Alakazam.R
NULL

#### Classes ####

#' S4 class defining a clone
#' 
#' \code{ChangeoClone} defines a common data structure for perform lineage recontruction
#' from Change-O data.
#' 
#' @slot     data      data.frame containing sequences and annotations. Contains the
#'                     columns \code{SEQUENCE_ID} and \code{SEQUENCE}, as well as any additional 
#'                     sequence-specific annotation columns.
#' @slot     clone     string defining the clone identifier.
#' @slot     germline  string containing the germline sequence for the clone.
#' @slot     v_gene    string defining the V segment gene call.
#' @slot     j_gene    string defining the J segment gene call.
#' @slot     junc_len  numeric junction length (nucleotide count).
#' 
#' @seealso  See \code{\link{makeChangeoClone}} and \code{\link{buildPhylipLineage}} for use.
#'           
#' @name         ChangeoClone-class
#' @rdname       ChangeoClone-class
#' @aliases      ChangeoClone
#' @exportClass  ChangeoClone
setClass("ChangeoClone", 
         slots=c(data="data.frame",
                 clone="character",
                 germline="character", 
                 v_gene="character", 
                 j_gene="character", 
                 junc_len="numeric"))

#### Data ####

#' Example Ig lineage trees
#'
#' A set of Ig lineage trees generated from the \code{extdata/ExampleDb.gz} file.
#'
#' @format   A list of igraph objects output by \code{\link{buildPhylipLineage}}.
"ExampleTrees"


#### Preprocessing functions ####

#' Generate a ChangeoClone object for lineage construction
#' 
#' \code{makeChangeoClone} takes a data.frame with Change-O style columns as input and 
#' masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
#' annotations associated with duplicate sequences. It returns a \code{ChangeoClone} 
#' object which serves as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the Change-O data for a clone. See Details
#'                        for the list of required columns and their default values.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    germ         name of the column containing germline DNA sequences. All entries 
#'                        in this column should be identical for any given clone, and they
#'                        must be multiple aligned with the data in the \code{seq} column.
#' @param    vcall        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    jcall        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    max_mask     maximum number of characters to mask at the leading and trailing
#'                        sequence ends. If \code{NULL} then the upper masking bound will 
#'                        be automatically determined from the maximum number of observed 
#'                        leading or trailing Ns amongst all sequences. If set to \code{0} 
#'                        (default) then masking will not be performed.
#' @param    text_fields  text annotation columns to retain and merge during duplicate removal.
#' @param    num_fields   numeric annotation columns to retain and sum during duplicate removal.
#' @param    seq_fields   sequence annotation columns to retain and collapse during duplicate 
#'                        removal. Note, this is distinct from the \code{seq} and \code{germ} 
#'                        arguments, which contain the primary sequence data for the clone
#'                        and should not be repeated in this argument.
#' @param    add_count    if \code{TRUE} add an additional annotation column called 
#'                        COLLAPSE_COUNT during duplicate removal that indicates the 
#'                        number of sequences that were collapsed.
#'
#' @return   A \code{\link{ChangeoClone}} object containing the modified clone.
#'
#' @details
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{vcall}, \code{jcall}, 
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "SEQUENCE_ID"}:           unique sequence identifier.
#'   \item  \code{seq      = "SEQUENCE_IMGT"}:         IMGT-gapped sample sequence.
#'   \item  \code{germ     = "GERMLINE_IMGT_D_MASK"}:  IMGT-gapped germline sequence.
#'   \item  \code{vcall    = "V_CALL"}:                V-segment allele call.
#'   \item  \code{jcall    = "J_CALL"}:                J-segment allele call.
#'   \item  \code{junc_len = "JUNCTION_LENGTH"}:       junction sequence length.
#'   \item  \code{clone    = "CLONE"}:                 clone identifier.
#' }
#' Additional annotation columns specified in the \code{text_fields}, \code{num_fields} 
#' or \code{seq_fields} arguments will be retained in the \code{data} slot of the return 
#' object, but are not required. If the input data.frame \code{data} already contains a 
#' column named \code{SEQUENCE}, which is not used as the \code{seq} argument, then that 
#' column will not be retained.
#' 
#' The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
#' However, all sequences (both observed and germline) must be multiple aligned using
#' some scheme for both proper duplicate removal and lineage reconstruction. 
#'
#' The value for the germline sequence, V-segment gene call, J-segment gene call, 
#' junction length, and clone identifier are determined from the first entry in the 
#' \code{germ}, \code{vcall}, \code{jcall}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'  
#' @seealso  Executes in order \code{\link{maskSeqGaps}}, \code{\link{maskSeqEnds}}
#'           and \code{\link{collapseDuplicates}}. 
#'           Returns a \code{\link{ChangeoClone}} object which serves as input to
#'           \code{\link{buildPhylipLineage}}.
#' 
#' @examples
#' # Example Change-O data.frame
#' df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
#'                  SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
#'                  V_CALL="Homsap IGKV1-39*01 F",
#'                  J_CALL="Homsap IGKJ5*01 F",
#'                  JUNCTION_LENGTH=2,
#'                  GERMLINE_IMGT_D_MASK="CCCCAGGG",
#'                  CLONE=1,
#'                  TYPE=c("IgM", "IgG", "IgG", "IgA"),
#'                  COUNT=1:4,
#'                  stringsAsFactors=FALSE)
#' 
#' # Without end masking
#' makeChangeoClone(df, text_fields="TYPE", num_fields="COUNT")
#'
#' # With end masking
#' makeChangeoClone(df, max_mask=3, text_fields="TYPE", num_fields="COUNT")
#'
#' @export
makeChangeoClone <- function(data, id="SEQUENCE_ID", seq="SEQUENCE_IMGT", 
                             germ="GERMLINE_IMGT_D_MASK", vcall="V_CALL", jcall="J_CALL",
                             junc_len="JUNCTION_LENGTH", clone="CLONE",
                             max_mask=0, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
                             add_count=TRUE) {
    # Check for valid fields
    check <- checkColumns(data, c(id, seq, germ, vcall, jcall, junc_len, clone, 
                                  text_fields, num_fields, seq_fields))
    if (check != TRUE) { stop(check) }
    
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, c(id, seq, text_fields, num_fields, seq_fields)]
    tmp_df[[seq]] <- maskSeqGaps(tmp_df[[seq]], outer_only=FALSE)
    tmp_df[[seq]] <- maskSeqEnds(tmp_df[[seq]], max_mask=max_mask, trim=FALSE)
    
    # Remove duplicates
    tmp_df <- collapseDuplicates(tmp_df, id=id, seq=seq, text_fields=text_fields, 
                                 num_fields=num_fields, seq_fields=seq_fields,
                                 add_count=add_count)
    
    # Define return object
    tmp_names <- names(tmp_df)
    if ("SEQUENCE" %in% tmp_names & seq != "SEQUENCE") {
        tmp_df <- tmp_df[, tmp_names != "SEQUENCE"]
        tmp_names <- names(tmp_df)
    }
    names(tmp_df)[tmp_names == seq] <- "SEQUENCE"
    clone <- new("ChangeoClone", 
                 data=as.data.frame(tmp_df),
                 clone=as.character(data[[clone]][1]),
                 germline=maskSeqGaps(data[[germ]][1], outer_only=FALSE), 
                 v_gene=getGene(data[[vcall]][1]), 
                 j_gene=getGene(data[[jcall]][1]), 
                 junc_len=data[[junc_len]][1])
    
    return(clone)
}


#### PHYLIP functions ####

# Create PHYLIP input files in a temporary folder
#
# @param   clone  a ChangeoClone object
# @param   path   a directory to store the write the output files to
# @return  a named vector translating SEQUENCE_ID (names) to PHYLIP taxa (values)
writePhylipInput <- function(clone, path) {
    # Define PHYLIP columns
    nseq <- nrow(clone@data)
    v1 <- c(sprintf('%-9s', nseq + 1),
            sprintf("%-9s", "Germline"), 
            sprintf("SAM%-6s", 1:nseq))
    v2 <- c(stri_length(clone@germline),
            clone@germline, 
            clone@data[["SEQUENCE"]])
    phy_df <- data.frame(v1, v2, stringsAsFactors=F)
    
    # Define names vector mapping taxa names to original sequence identifiers
    id_map <- setNames(gsub("^\\s+|\\s+$", "", v1[-(1:2)]), clone@data[["SEQUENCE_ID"]])
    
    # Create PHYLIP input file
    write.table(phy_df, file=file.path(path, "infile"), 
                quote=F, sep=" ", col.names=F, row.names=F)    
    
    return(id_map)
}


# Run PHYLIP dnapars application
#
# @param   path          temporary directory containing infile.
# @param   dnapars_exec  path to the dnapars executable.
# @param   verbose       if TRUE suppress phylip console output.
# @return  TRUE if phylip ran successfully and FALSE otherwise
runPhylip <- function(path, dnapars_exec, verbose=FALSE) {
    # Expand shell variables
    dnapars_exec <- path.expand(dnapars_exec)
    
    # Remove old files
    if (file.exists(file.path(path, "outfile"))) { file.remove(file.path(path, "outfile")) }
    if (file.exists(file.path(path, "outtree"))) { file.remove(file.path(path, "outtree")) }    
    
    # Set platform specific options
    if (.Platform$OS.type == "windows") { 
        quiet_params <- list(ignore.stdout=TRUE, ignore.stderr=TRUE)
        invoke <- shell
    } else { 
        quiet_params <- list(stdout=FALSE, stderr=FALSE)
        invoke <- system2
    } 
    
    # Set dnapars options
    phy_options <- c("S", "Y", "I", "4", "5", ".")
    params <- list(dnapars_exec, input=c(phy_options, "Y"), wait=TRUE)
    if (!verbose) {
        params <- append(params, quiet_params)
    }
    
    # Call phylip
    wd <- getwd()
    setwd(path)
    status <- tryCatch(do.call(invoke, params), error=function(e) e)
    setwd(wd)
    
    # Return TRUE if phylip ran successfully
    invisible(status == 0)
}


# Reads in the PHYLIP outfile
#
# @param   path  the temporary folder containing the dnapars outfile
# @return  a character vector with each item as a line in the outfile
readPhylipOutput <- function(path) {
    phylip_out <- scan(file.path(path, "outfile"), what="character", sep="\n", 
                       blank.lines.skip=FALSE, strip.white=FALSE, quiet=TRUE)
    return(phylip_out)
}


# Test for successful PHYLIP dnapars run by checking the outfile
#
# @param   phylip_out  a character vector returned by readPhylipOut
# @return  TRUE if trees built 
#          FALSE if no trees built
checkPhylipOutput <- function(phylip_out) {
    # Check for failed tree build
    result <- !(any(grepl('-1 trees in all found', phylip_out)))
    
    return(result)
}


# Extracts inferred sequences from PHYLIP dnapars outfile
#
# @param   phylip_out   a character vector returned by readPhylipOutput
# @return  a list containing an id vector, a sequence vector and an annotation data.frame
getPhylipInferred <- function(phylip_out) {
    # Process dnapars output
    seq_start <- min(grep("From\\s+To\\s+Any Steps\\?\\s+State at upper node", 
                          phylip_out, perl=T, fixed=F))
    seq_empty <- grep("^\\s*$", phylip_out[seq_start:length(phylip_out)], perl=T, fixed=F)
    seq_len <- seq_empty[min(which(seq_empty[-1] == (seq_empty[-length(seq_empty)] + 1)))]
    seq_block <- paste(phylip_out[(seq_start + 2):(seq_start + seq_len - 2)], collapse="\n")
    seq_df <- read.table(textConnection(seq_block), as.is=T, fill=T, blank.lines.skip=F)
    
    # Correct first line of block and remove blank rows
    fix.row <- c(1, which(is.na(seq_df[,1])) + 1)
    end_col <-  ncol(seq_df) - 2
    #seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:5], stringsAsFactors=F)
    seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:end_col], stringsAsFactors=F)
    seq_df <- seq_df[-(fix.row[-1] - 1), ]
    
    # Create data.frame of inferred sequences
    inferred_num <- unique(grep("^[0-9]+$", seq_df[, 2], value=T))
    inferred_seq <- sapply(inferred_num, function(n) { paste(t(as.matrix(seq_df[seq_df[, 2] == n, -c(1:3)])), collapse="") })
    
    return(data.frame(SEQUENCE_ID=paste0("Inferred", inferred_num), SEQUENCE=inferred_seq))
}


# Extracts graph edge list from a PHYLIP dnapars outfile
#
# @param   phylip_out  character vector returned by readPhylipOutput
# @param   id_map      named vector of PHYLIP taxa names (values) to sequence 
#                      identifiers (names) that will be translated. If NULL
#                      no taxa name translation is performed
# @return  a data.frame of edges with columns (from, to, weight)
getPhylipEdges <- function(phylip_out, id_map=NULL) {
    # Process dnapars output
    edge_start <- min(grep('between\\s+and\\s+length', phylip_out, 
                           perl=TRUE, fixed=FALSE))
    edge_len <- min(grep('^\\s*$', phylip_out[edge_start:length(phylip_out)], 
                         perl=TRUE, fixed=FALSE))
    edge_block <- paste(phylip_out[(edge_start + 2):(edge_start + edge_len - 2)], collapse='\n')
    edge_df <- read.table(textConnection(edge_block), col.names=c('from', 'to', 'weight'), 
                          as.is=TRUE)

    # Modify inferred taxa names to include "Inferred"
    inf_map <- unique(grep("^[0-9]+$", c(edge_df$from, edge_df$to), value=T))
    names(inf_map) <- paste0("Inferred", inf_map)
    edge_df$from <- translateStrings(edge_df$from, inf_map)
    edge_df$to <- translateStrings(edge_df$to, inf_map)
    
    if (!is.null(id_map)) {
        # Reassign PHYLIP taxa names to sequence IDs
        edge_df$from <- translateStrings(edge_df$from, id_map)
        edge_df$to <- translateStrings(edge_df$to, id_map)
    }
    
    return(edge_df)
}


# Modify edges of phylip output
#
# @param   edges     data.frame of edges returned by getPhylipEdges
# @param   clone     a ChangeoClone object containg sequence data
# @param   dist_mat  DNA character distance matrix
# @return  a list of modified edges data.frame and clone object
modifyPhylipEdges <- function(edges, clone, dist_mat=getDNAMatrix(gap=0)) {
    # Move germline to root position
    germ_idx <- which(edges$to == "Germline")
    edges[germ_idx, c('from', 'to')] <- edges[germ_idx, c('to', 'from')]
    
    # Calculate edge mutations
    for (i in 1:nrow(edges)) {
        if (edges$from[i] == "Germline") {
            seq1 <- clone@germline
        } else {
            seq1 <- clone@data[["SEQUENCE"]][clone@data[["SEQUENCE_ID"]] == edges$from[i]]
        }
        seq2 <- clone@data[["SEQUENCE"]][clone@data[["SEQUENCE_ID"]] == edges$to[i]]
        edges$weight[i] <- getSeqDistance(seq1, seq2, dist_mat)        
    }
    
    # Find rows zero weight edges with inferred parent nodes
    remove_row <- which(edges$weight == 0 & 
                        edges$from != "Germline" & 
                        grepl('^Inferred\\d+$', edges$from))
    
    # Replace inferred parent nodes with child nodes when edge weight is zero
    while (length(remove_row) > 0) {
        # Remove first node with zero distance to parent
        r <- remove_row[1]
        r_idx <- which(edges[c('from', 'to')] == edges$from[r], arr.ind=T)
        edges[r_idx] <- edges$to[r]
        
        # Recalculate edge weights for modified rows
        r_mod <- r_idx[, 1][r_idx[, 1] != r]
        for (i in r_mod) {
            if (edges$from[i] == "Germline") {
                seq1 <- clone@germline
            } else {
                seq1 <- clone@data[["SEQUENCE"]][clone@data[["SEQUENCE_ID"]] == edges$from[i]]
            }
            seq2 <- clone@data[["SEQUENCE"]][clone@data[["SEQUENCE_ID"]] == edges$to[i]]
            edges$weight[i] <- getSeqDistance(seq1, seq2, dist_mat)      
        }
        
        # Remove row
        edges <- edges[-r, ]
        
        # Re-determine rows to remove
        remove_row <- which(edges$weight == 0 & 
                            edges$from != "Germline" & 
                            grepl('^Inferred\\d+$', edges$from))      
    }
    
    # Remove rows from clone
    keep_clone <- clone@data[["SEQUENCE_ID"]] %in% unique(c(edges$from, edges$to))
    clone@data <- as.data.frame(clone@data[keep_clone, ])
    
    return(list(edges=edges, clone=clone))
}

# Convert edge data.frame and clone object to igraph graph object
#
# @param   edges  data.frame of edges returned by getPhylipEdges
# @param   clone  a ChangeoClone object containg sequence data
# @return  an igraph graph object
phylipToGraph <- function(edges, clone) {
    # Create igraph object
    g <- igraph::graph_from_data_frame(edges, directed=T)
    
    # Add germline sequence
    germ_idx <- which(igraph::V(g)$name == "Germline")
    g <- igraph::set_vertex_attr(g, "sequence", index=germ_idx, clone@germline)
    
    # Add sample sequences and names
    clone_idx <- match(clone@data[["SEQUENCE_ID"]], igraph::V(g)$name) 
    g <- igraph::set_vertex_attr(g, "sequence", index=clone_idx, clone@data[["SEQUENCE"]])
    
    # Add annotations
    ann_fields <- names(clone@data)[!(names(clone@data) %in% c("SEQUENCE_ID", "SEQUENCE"))]
    for (n in ann_fields) {
        g <- igraph::set_vertex_attr(g, n, index=germ_idx, NA)
        g <- igraph::set_vertex_attr(g, n, index=clone_idx, clone@data[[n]])
    }
    
    # Add edge and vertex labels
    igraph::V(g)$label <- igraph::V(g)$name
    igraph::E(g)$label <- igraph::E(g)$weight
    
    # Add graph attributes
    g$clone <- clone@clone
    g$v_gene <- clone@v_gene
    g$j_gene <- clone@j_gene
    g$junc_len <- clone@junc_len
    
    return(g)
}


#' Infer an Ig lineage using PHYLIP
#' 
#' \code{buildPhylipLineage} reconstructs an Ig lineage via maximum parsimony using the 
#' dnapars application of the PHYLIP package.
#' 
#' @param    clone         \code{\link{ChangeoClone}} object containing clone data.
#' @param    dnapars_exec  path to the PHYLIP dnapars executable.
#' @param    rm_temp       if \code{TRUE} delete the temporary directory after running dnapars;
#'                         if \code{FALSE} keep the temporary directory.
#' @param    verbose       if \code{FALSE} suppress the output of dnapars; 
#'                         if \code{TRUE} STDOUT and STDERR of dnapars will be passed to 
#'                         the console.
#'                                                
#' @return   An igraph \code{graph} object defining the Ig lineage tree. Each unique input 
#'           sequence in \code{clone} is a vertex of the tree, with additional vertices being
#'           either the germline (root) sequences or inferred intermediates. The \code{graph} 
#'           object has the following attributes.
#'           
#'           Vertex attributes:
#'           \itemize{
#'             \item  \code{name}:      value in the \code{SEQUENCE_ID} column of the \code{data} 
#'                                      slot of the input \code{clone} for observed sequences. 
#'                                      The germline (root) vertex is assigned the name 
#'                                      "Germline" and inferred intermediates are assigned
#'                                      names with the format {"Inferred1", "Inferred2", ...}.
#'             \item  \code{sequence}:  value in the \code{SEQUENCE} column of the \code{data} 
#'                                      slot of the input \code{clone} for observed sequences.
#'                                      The germline (root) vertex is assigned the sequence
#'                                      in the \code{germline} slot of the input \code{clone}.
#'                                      The sequence of inferred intermediates are extracted
#'                                      from the dnapars output.
#'             \item  \code{label}:     same as the \code{name} attribute.
#'           }
#'           Additionally, each other column in the \code{data} slot of the input 
#'           \code{clone} is added as a vertex attribute with the attribute name set to 
#'           the source column name. For the germline and inferred intermediate vertices,
#'           these additional vertex attributes are all assigned a value of \code{NA}.
#'           
#'           Edge attributes:
#'           \itemize{
#'             \item  \code{weight}:    Hamming distance between the \code{sequence} attributes
#'                                      of the two vertices.
#'             \item  \code{label}:     same as the \code{weight} attribute.
#'           }
#'           Graph attributes:
#'           \itemize{
#'             \item  \code{clone}:     clone identifier from the \code{clone} slot of the
#'                                      input \code{ChangeoClone}.
#'             \item  \code{v_gene}:    V-segment gene call from the \code{v_gene} slot of 
#'                                      the input \code{ChangeoClone}.
#'             \item  \code{j_gene}:    J-segment gene call from the \code{j_gene} slot of 
#'                                      the input \code{ChangeoClone}.
#'             \item  \code{junc_len}:  junction length (nucleotide count) from the 
#'                                      \code{junc_len} slot of the input \code{ChangeoClone}.
#'           }
#'           
#' @details
#' \code{buildPhylipLineage} builds the lineage tree of a set of unique Ig sequences via
#' maximum parsimony through an external call to the dnapars application of the PHYLIP
#' package. dnapars is called with default algorithm options, except for the search option, 
#' which is set to "Rearrange on one best tree". The germline sequence of the clone is used 
#' for the outgroup. 
#' 
#' Following tree construction using dnapars, the dnapars output is modified to allow
#' input sequences to appear as internal nodes of the tree. Intermediate sequences 
#' inferred by dnapars are replaced by children within the tree having a Hamming distance 
#' of zero from their parent node. The distance calculation allows IUPAC ambiguous 
#' character matches, where an ambiguous character has distance zero to any character in 
#' the set of characters it represents. Distance calculation and movement of child nodes 
#' up the tree is repeated until all parent-child pairs have a distance greater than zero 
#' between them. The germline sequence (outgroup) is moved to the root of the tree and
#' excluded from the node replacement processes, which permits the trunk of the tree to be
#' the only edge with a distance of zero. Edge weights of the resultant tree are assigned 
#' as the distance between each sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Felsenstein J. PHYLIP - Phylogeny Inference Package (Version 3.2). 
#'            Cladistics. 1989 5:164-166.
#'   \item  Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
#'            sclerosis brain mature in the draining cervical lymph nodes. 
#'            Sci Transl Med. 2014 6(248):248ra107.
#' }
#'   
#' @seealso  Takes as input a \code{\link{ChangeoClone}}. 
#'           Temporary directories are created with \code{\link{makeTempDir}}.
#'           Distance is calculated using \code{\link{getSeqDistance}}. 
#'           See \code{\link{igraph}} and \code{\link{igraph.plotting}} for working 
#'           with igraph \code{graph} objects. 
#'
#' @examples
#' \dontrun{
#' # Load example data
#' file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Preprocess clone
#' clone <- subset(df, CLONE == 164)
#' clone <- makeChangeoClone(clone, text_fields=c("SAMPLE", "ISOTYPE"), num_fields="DUPCOUNT")
#' 
#' # Run PHYLIP and process output
#' dnapars_exec <- "~/apps/phylip-3.69/dnapars"
#' graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
#' 
#' # Plot graph with a tree layout
#' library(igraph)
#' ly <- layout_as_tree(graph, root="Germline", circular=F, flip.y=T)
#' plot(graph, layout=ly)
#' }
#' 
#' @export
buildPhylipLineage <- function(clone, dnapars_exec, rm_temp=FALSE, verbose=FALSE) {
    # Check clone size
    if (nrow(clone@data) < 2) {
        warning("Clone ", clone@clone, " was skipped as it does not contain at least 
                2 unique sequences")
        return(NULL)
    }
    
    # Check fields
    seq_len = unique(stri_length(clone@data[["SEQUENCE"]]))
    germ_len = ifelse(length(clone@germline) == 0, 0, stri_length(clone@germline))
    if(germ_len == 0) {
        stop("Clone ", clone@clone, "does not contain a germline sequence.")
    }
    if(length(seq_len) != 1) {
        stop("Clone ", clone@clone, "does not contain sequences of equal length.")
    }
    if(seq_len != germ_len) {
        stop("The germline and input sequences are not the same length for clone ", clone@clone)
    }
    
    # Check dnapars access
    if (file.access(dnapars_exec, mode=1) == -1) {
        stop("The file ", dnapars_exec, " cannot be executed.")
    }
    
    # Create temporary directory
    temp_path <- makeTempDir(paste0(clone@clone, "-phylip"))
    if (verbose) {
        cat("TEMP_DIR> ", temp_path, "\n", sep="")
    }
    
    # Run PHYLIP
    id_map <- writePhylipInput(clone, temp_path)
    runPhylip(temp_path, dnapars_exec, verbose=verbose)
    phylip_out <- readPhylipOutput(temp_path)
    
    # Remove temporary directory
    if (rm_temp) {
        unlink(temp_path, recursive=TRUE)
    }
    
    # Check output for trees
    if (!checkPhylipOutput(phylip_out)) {
        warning('PHYLIP failed to generate trees for clone ', clone)
        return(NULL)
    }
    
    # Extract inferred sequences from PHYLIP output
    inf_df <- getPhylipInferred(phylip_out)
    clone@data <- as.data.frame(bind_rows(clone@data, inf_df))

    # Extract edge table from PHYLIP output 
    edges <- getPhylipEdges(phylip_out, id_map=id_map)
    
    # Modify PHYLIP tree to remove 0 distance edges
    mod_list <- modifyPhylipEdges(edges, clone)

    # Convert edges and clone data to igraph graph object
    graph <- phylipToGraph(mod_list$edges, mod_list$clone)
    
    return(graph)
}