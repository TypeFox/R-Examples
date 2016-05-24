# Alakazam package documentation and import directives

#' The alakazam package
#' 
#' \code{alakazam} in a member of the Change-O suite of tools and serves five main 
#' purposes:
#' \itemize{
#'   \item  Providing core functionality for other R packages in the Change-O suite. This
#'          includes common tasks such as file I/O, basic DNA sequence manipulation, and
#'          interacting with V(D)J segment and gene annotations.
#'   \item  Providing an R interface for interacting with the output of the pRESTO 
#'          tool suite.
#'   \item  Performing lineage reconstruction on clonal populations of immunoglobulin 
#'          (Ig) sequences. 
#'   \item  Performing clonal abundance and diversity analysis on lymphocyte repertoires.
#'   \item  Performing physicochemical property analyses of lymphocyte receptor sequences.
#' }
#' For additional details regarding the use of the \code{alakazam} package see the 
#' vignettes:\cr
#' \code{browseVignettes("alakazam")}
#' 
#' @section  File I/O:
#' \itemize{
#'   \item  \code{\link{readChangeoDb}}:        Input Change-O style files.
#'   \item  \code{\link{writeChangeoDb}}:       Output Change-O style files.
#' }
#' 
#' @section  Sequence cleaning:
#' \itemize{
#'   \item  \code{\link{maskSeqEnds}}:          Mask ragged ends.
#'   \item  \code{\link{maskSeqGaps}}:          Mask gap characters.
#'   \item  \code{\link{collapseDuplicates}}:   Remove duplicate sequences.
#' }
#' 
#' @section  Lineage reconstruction:
#' \itemize{
#'   \item  \code{\link{makeChangeoClone}}:     Clean sequences for lineage reconstruction.
#'   \item  \code{\link{buildPhylipLineage}}:   Perform lineage reconstruction of Ig sequences.
#' }
#' 
#' @section  Diversity analysis:
#' \itemize{
#'   \item  \code{\link{countClones}}:          Calculate clonal abundance.
#'   \item  \code{\link{estimateAbundance}}:    Infer complete clonal abundance distribution with
#'                                              confidence intervals.
#'   \item  \code{\link{rarefyDiversity}}:      Generate clonal diversity curves.
#'   \item  \code{\link{testDiversity}}:        Test significance of clonal diversity scores.
#'   \item  \code{\link{plotAbundance}}:        Plot clone size distribution as a rank-abundance 
#'                                              curve.
#'   \item  \code{\link{plotDiversityCurve}}:   Plot clonal diversity curves.
#' }
#' 
#' @section  Ig and TCR sequence annotation:
#' \itemize{
#'   \item  \code{\link{countGenes}}:           Calculate Ig and TCR allele, gene and family usage.
#'   \item  \code{\link{extractVRegion}}:       Extract CDRs and FWRs sub-sequences.
#'   \item  \code{\link{getAllele}}:            Get V(D)J allele names.
#'   \item  \code{\link{getGene}}:              Get V(D)J gene names.
#'   \item  \code{\link{getFamily}}:            Get V(D)J family names.
#' }
#' 
#' @section  Sequence distance calculation:
#' \itemize{
#'   \item  \code{\link{getSeqDistance}}:       Calculate Hamming distance between sequences.
#'   \item  \code{\link{getSeqMatrix}}:         Calculate a matrix of pairwise Hamming distances 
#'                                              for a sequence set.
#'   \item  \code{\link{testSeqEqual}}:         Test sequences for equality.
#' }
#' 
#' @section  Amino acid propertes:
#' \itemize{
#'   \item  \code{\link{translateDNA}}:         Translate DNA sequences to amino acid sequences.
#'   \item  \code{\link{aminoAcidProperties}}:  Calculate various physicochemical properties of amino acid 
#'                                              sequences.
#'   \item  \code{\link{countPatterns}}:        Count patterns in sequences.
#'                                              
#' }
#' 
#' @section  General data manipulation:
#' \itemize{
#'   \item  \code{\link{translateStrings}}:     Perform multiple string replacement operations.
#' } 
#' 
#' @name     alakazam
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Vander Heiden JA, Yaari G, et al. pRESTO: a toolkit for processing 
#'            high-throughput sequencing raw reads of lymphocyte receptor repertoires. 
#'            Bioinformatics. 2014 30(13):1930-2.
#'   \item  Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
#'            sclerosis brain mature in the draining cervical lymph nodes. 
#'            Sci Transl Med. 2014 6(248):248ra107.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#'   \item  Gupta NT, Vander Heiden JA, et al. Change-O: a toolkit for analyzing 
#'            large-scale B cell immunoglobulin repertoire sequencing data.
#'            Under review.
#' }
#' 
#' @import      ggplot2
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  dplyr     do n desc %>%
#'                        as_data_frame data_frame data_frame_
#'                        bind_cols bind_rows combine
#'                        filter filter_ select select_ arrange arrange_
#'                        group_by group_by_ ungroup
#'                        mutate mutate_ transmute transmute_
#'                        rename rename_ summarize summarize_
#' @importFrom  igraph    V E graph_from_data_frame
#'                        vertex_attr set_vertex_attr
#' @importFrom  lazyeval  interp
#' @importFrom  scales    log2_trans log10_trans trans_breaks trans_format
#'                        math_format percent scientific
#' @importFrom  seqinr    translate
#' @importFrom  stats     na.omit setNames ecdf sd cor cov median mad
#'                        dbinom pbinom qbinom rbinom
#'                        dnorm pnorm qnorm rnorm
#'                        dmultinom rmultinom
#' @importFrom  stringi   stri_dup stri_flatten stri_join stri_length 
#'                        stri_count_boundaries stri_count_regex 
#'                        stri_extract_all_regex stri_extract_first_regex  
#'                        stri_replace_all_regex stri_replace_first_regex
NULL


#### Sysdata ####

# 1x20 vector of default amino acid hydropathy scores
# HYDROPATHY_KYTJ82

# 1x20 vector of default amino acid bulkiness scores
# BULKINESS_ZIMJ68

# 1x20 vector of default amino acid polarity scores
# POLARITY_GRAR74

# 1x7 vector of default amino acid pK values
# PK_EMBOSS
