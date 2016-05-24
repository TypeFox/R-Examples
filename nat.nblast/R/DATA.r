#' Scoring matrices for neuron similarities in FCWB template brain
#'
#' Scoring matrices quantify the log2 odds ratio that a segment pair with a
#' given distance and absolute dot product come from a pair of neurons of the
#' same type, rather than unrelated neurons.
#'
#' These scoring matrices were generated using all by all pairs from 150 DL2
#' antennal lobe projection neurons from the \url{http://flycircuit.tw} dataset
#' and 5000 random pairs from the same dataset.
#'
#' \itemize{
#'
#' \item \code{smat.fcwb} was trained using nearest-neighbour distance and the
#' tangent vector defined by the first eigen vector of the k=5 nearest
#' neighbours.
#'
#' \item \code{smat_alpha.fcwb} was defined as for \code{smat.fcwb} but weighted
#' by the factor \code{alpha} defined as (l1-l2)/(l1+l2+l3) where l1,l2,l3 are
#' the three eigen values.
#'
#' }
#'
#' Most work on the flyircuit dataset has been carried out using the
#' \code{smat.fcwb} scoring matrix although the \code{smat_alpha.fcwb} matrix
#' which emphasises the significance of matches between linear regions of the
#' neuron (such as axons) may have some advantages.
#'
#' @name smat.fcwb
#' @aliases smat_alpha.fcwb
#' @docType data
NULL

#' 20 traced Drosophila neurons from Chiang et al 2011
#'
#' This R list (which has additional class \code{neuronlist}) contains 15
#' skeletonized \emph{Drosophila} neurons as \code{dotprops} objects. Original
#' data is due to Chiang et al. [1], who have generously shared their raw data
#' at \url{http://flycircuit.tw}. Automated tracing of neuron skeletons was
#' carried out by Lee et al [2]. Image registration and further processing was
#' carried out by Greg Jefferis, Marta Costa and James Manton[3].
#' @name fctraces20
#' @docType data
#' @references [1] Chiang A.S., Lin C.Y., Chuang C.C., Chang H.M., Hsieh C.H.,
#'   Yeh C.W., Shih C.T., Wu J.J., Wang G.T., Chen Y.C., Wu C.C., Chen G.Y.,
#'   Ching Y.T., Lee P.C., Lin C.Y., Lin H.H., Wu C.C., Hsu H.W., Huang Y.A.,
#'   Chen J.Y., et al. (2011). Three-dimensional reconstruction of brain-wide
#'   wiring networks in Drosophila at single-cell resolution. Curr Biol 21 (1),
#'   1--11.
#'
#'   [2] P.-C. Lee, C.-C. Chuang, A.-S. Chiang, and Y.-T. Ching. (2012).
#'   High-throughput computer method for 3d neuronal structure reconstruction
#'   from the image stack of the Drosophila brain and its applications. PLoS
#'   Comput Biol, 8(9):e1002658, Sep 2012. doi: 10.1371/journal.pcbi.1002658.
#'
#'   [3] NBLAST: Rapid, sensitive comparison of neuronal structure and
#'   construction of neuron family databases. Marta Costa, Aaron D. Ostrovsky,
#'   James D. Manton, Steffen Prohaska, Gregory S.X.E. Jefferis. bioRxiv doi:
#'   http://dx.doi.org/10.1101/006346.
NULL
