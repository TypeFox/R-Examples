# shazam package documentation and import directives

#' The shazam package
#'
#' Provides tools for advanced anaylisis of immunoglobulin (Ig) somatic hypermutation 
#' (SHM), including BASELINe, a novel method for quantifying antigen-driven selection in 
#' high-throughput Ig sequencing data.
#' 
#' Dramatic improvements in high-throughput sequencing technologies now enable 
#' large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
#' antigen-receptor proteins located on the surface of T and B lymphocytes.
#' 
#' The \code{shazam} package provides tools for advanced analysis of Ig sequences following 
#' germline segment assignment. Namely, the analysis of SHM. 
#' Which includes:
#'  \itemize{
#'      \item   Statistical analysis of SHM patterns \cr
#'              Computational models and analyses of SHM have separated the process 
#'              into two independent components: 
#'              \enumerate{
#'                  \item  A mutability model that defines where mutations occur.
#'                  \item  A nucleotide substitution model that defines the resulting mutation.
#'              }
#'              Collectively these are what form the targeting model of SHM. \code{shazam} 
#'              provides tools to build these mutability and substitution (i.e. targeting) 
#'              models.
#'                  
#'      \item   BASELINe \cr
#'              Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
#'              novel method for quantifying antigen-driven selection in high-throughput
#'              Ig sequence data. The targeting model created using \code{shazam} is used 
#'              to estimate the null distribution of expected mutation frequencies in 
#'              BASELINe.
#'              
#'      \item   Distance calculations \cr
#'              Based on the underlying SHM targeting (calculated using \code{shazam}) one 
#'              can compute evolutionary distances between sequences or groups of 
#'              sequences. This information is particularly useful in understanding and 
#'              defining clonal relationships.
#'  }
#' 
#' Below are the functions in \code{shazam} broken down by the three main tasks described
#' above:
#' 
#' @section  Targeting models:
#' \itemize{
#'   \item  \link{createTargetingModel}:     Build a 5-mer targeting model.
#'   \item  \link{plotMutability}:           Plot 5-mer mutability rates.
#' }
#' 
#' @section  Mutational profiling:
#' \itemize{
#'   \item  \link{collapseByClone}:    Build clonal consensus sequence.
#'   \item  \link{calcDBObservedMutations}:  Compute observed mutation counts.
#'   \item  \link{calcDBExpectedMutations}:  Compute expected mutation frequencies.
#' }
#'
#' @section  Selection analysis:
#' \itemize{
#'   \item  \link{calcBaseline}:             Calculate the BASELINe probability
#'                                           density functions (PDFs).
#'   \item  \link{groupBaseline}:            Combine PDFs from sequences grouped
#'                                           by biological or experimental relevance.
#'   \item  \link{summarizeBaseline}:        Compute summary statistics from BASELINe PDFs.
#'   \item  \link{plotBaselineDensity}:      Plot the probability density functions
#'                                           resulting from selection analysis.
#'   \item  \link{plotBaselineSummary}:      Plot summary stastistics resulting from 
#'                                           selection analysis.
#' }
#'
#' @section  Distance profiling:
#' \itemize{
#'   \item  \link{distToNearest}:            Tune clonal assignment thresholds by calculating 
#'                                           distances to nearest-neighbors.
#'   \item  \link{calcTargetingDistance}:    Construct a nucleotide distance matrix from a 
#'                                           5-mer targeting model.
#' }
#'
#' @name     shazam
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'  }
#' 
#' @import   alakazam
#' @import   ggplot2
#' @import   graphics
#' @import   methods
#' @import   utils
#' @importFrom  data.table  data.table setkey setkeyv
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       do n desc %>%
#'                          as_data_frame data_frame data_frame_
#'                          bind_cols bind_rows combine
#'                          filter filter_ select select_ arrange arrange_
#'                          group_by group_by_ ungroup
#'                          mutate mutate_ transmute transmute_
#'                          rename rename_ summarize summarize_
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  lazyeval    interp
#' @importFrom  scales      log2_trans log10_trans trans_breaks trans_format
#'                          math_format percent scientific
#' @importFrom  tidyr       gather gather_ spread spread_
#' @importFrom  iterators   icount
#' @importFrom  SDMTools    wt.sd
#' @importFrom  seqinr      c2s s2c words translate
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          approx convolve weighted.mean
#'                          dbeta pbeta qbeta rbeta
#' @importFrom  stringi     stri_dup stri_flatten stri_join stri_length
#'                          stri_count_boundaries stri_count_regex 
#'                          stri_extract_all_regex stri_extract_first_regex  
#'                          stri_replace_all_regex stri_replace_first_regex
NULL


#### Sysdata ####

# Ordered nucleotide character set
# NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

# IMGT V segment length
# VLENGTH <- 312

# 5x312 logical matrix of CDR positions
# CDR_Nuc_Mat

# 5x312 logical matrix of FWR positions
# FWR_Nuc_Mat

# 12x216 matrix of replacement and silent mutation permutations
# CODON_TABLE

# 1x24 vector of amino acid charge classes
# AMINO_ACIDS_CHARGE

# 1x24 vector of amino acid hydropathy classes
# AMINO_ACIDS_HYDROPATHY

# 1x24 vector of amino acid polarity classes
# AMINO_ACIDS_POLARITY

# TODO: What is this?
# CONST_I

# TODO: And what is this?
# BAYESIAN_FITTED

# Add built-in variables to global variables environment
utils::globalVariables(c("M1NDistance", "HS1FDistance", 
                         "U5NModel", "HS5FModel"), package="shazam")