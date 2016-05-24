#' Neuron similarity, search and clustering tools
#'
#' @section Similarity and search:
#'
#'   The main entry point for similarity and search functions is
#'   \code{\link{nblast}}. Traced neurons will normally be converted to the
#'   \code{\link[nat]{dotprops}} format for search. When multiple neurons are
#'   compared they should be in a \code{\link[nat]{neuronlist}} object.
#'
#'   The current nblast version (2) depends on a scoring matrix. Default
#'   matrices trained using \emph{Drosophila} neurons in the FCWB template brain
#'   space are distributed with this package (see \code{\link{smat.fcwb}}); see
#'   \bold{Scoring Matrices} section below for creating new scoring matrices.
#'
#'   \code{nblast} makes use of a more flexible but more complicated function
#'   \code{NeuriteBlast} which includes several additional options. The function
#'   \code{WeightedNNBasedLinesetMatching} provides the primitive functionality
#'   of finding the nearest neighbour distances and absolute dot products for
#'   two sets of segments. Neither of these functions are intended for end use.
#'
#'   Calculating all by all similarity scores is facilitated by the
#'   \code{\link{nblast_allbyall}} function which can take either a neuronlist
#'   as input or a character vector naming (a subset) of neurons in a (large)
#'   neuronlist. The neuronlist containing the input neurons should be resident
#'   in memory i.e. not the \code{neuronlistfh}
#'
#' @section Clustering:
#'
#'   Once an all by all similarity score matrix is available it can be used as
#'   the input to a variety of clustering algorithms. \code{\link{nhclust}}
#'   provides a convenient wrapper for R's hierarchical clustering function
#'   \code{\link{hclust}}. If you wish to use another clustering function, then
#'   you can use the \code{\link{sub_dist_mat}} to convert a raw similarity
#'   score matrix into a normalised distance matrix (or R \code{\link{dist}}
#'   object) suitable for clustering. If you need a similarity matrix or want to
#'   modify the normalisation then you can use \code{\link{sub_score_mat}}.
#'
#'   Note tha raw nblast scores are not symmetric (i.e. S(A,B) is not equal to
#'   S(B,A)) so before clustering we construct a symmetric similarity/distance
#'   matrix \code{1/2 * ( S(A,B)/S(A,A) + S(B,A)/S(B,B) )}. See
#'   \code{\link{sub_score_mat}}'s documentation for details.
#'
#' @section Cached scores:
#'
#'   Although nblast is fast and can be parallelised, it makes sense to cache to
#'   disk all by all similarity scores for a group of neurons that will be
#'   subject to repeated clustering or other analysis. The matrix can simply be
#'   saved to disk and then reloaded using base R functions like
#'   \code{\link{save}} and \code{\link{load}}. \code{\link{sub_score_mat}} and
#'   \code{\link{sub_dist_mat}} can be used to extract a subset of scores from
#'   this raw score matrix. For large matrices, the \code{bigmemory} or
#'   \code{ff} packages allow matrices to be stored on disk and portions loaded
#'   into memory on demand. \code{\link{sub_score_mat}} and
#'   \code{\link{sub_dist_mat}} work equally well for regular in-memory matrices
#'   and these disk-backed matrices.
#'
#'   To give an example, for 16,129 neurons from the flycircuit.tw dataset, the
#'   260,144,641 comparisons took about 250 hours of compute time (half a day on
#'   ~20 cores). When saved to disk as single precision (i.e. 4 bytes per score)
#'   \code{ff} matrix they occupy just over 1Gb.
#'
#' @section Calculating scoring matrices:
#'
#'   The nblast algorithm depends on appropriately calibrated scoring matrices.
#'   These encapsulate the log odds ratio that a pair of segments come from two
#'   structurally related neurons rather than two unrelated neurons, given the
#'   observed distance and absolute dot product of the two segments. Scoring
#'   matrices can be constructed using the \code{\link{create_scoringmatrix}}
#'   function, supplying a set of matching neurons and a set of non-matching
#'   neurons. See the \code{create_scoringmatrix} documentation for links to
#'   lower-level functions that provide finer control over construction of the
#'   scoring matrix.
#'
#' @section Package Options:
#'
#'   There is one package option \code{nat.nblast.defaultsmat} which is
#'   \code{NULL} by default, but could for example be set to one of the scoring
#'   matrices included with the package such as code{"smat.fcwb"} or to a new
#'   user-constructed matrix.
#'
#' @references Costa, M., Ostrovsky, A.D., Manton, J.D., Prohaska, S., and
#'   Jefferis, G.S.X.E. (2014). NBLAST: Rapid, sensitive comparison of neuronal
#'   structure and construction of neuron family databases. Biorxiv preprint.
#'   \href{http://dx.doi.org/10.1101/006346}{doi: 10.1101/006346}.
#'
#' @name nat.nblast-package
#' @aliases nat.nblast
#' @docType package
#' @import methods
#' @keywords package
#' @seealso \code{\link{nblast}}, \code{\link{smat.fcwb}},
#'   \code{\link{nhclust}}, \code{\link{sub_dist_mat}},
#'   \code{\link{sub_score_mat}}, \code{\link{create_scoringmatrix}}
NULL
