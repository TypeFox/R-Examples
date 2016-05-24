
#' @name DynClust-package
#' @title Denoising and clustering for dynamical image sequence (2D or 3D)+T
#'
#' @description DynClust is a two-stage procedure for the denoising and clustering
#' of stack of noisy images acquired over time. Clustering only assumes that the data
#' contain an unknown but small number of dynamic features. The method first denoises
#' the signals using local spatial and full temporal information. The clustering step
#' uses the previous output to aggregate voxels based on the knowledge of their spatial
#' neighborhood. Both steps use a single keytool based on the statistical comparison
#' of the difference of two signals with the null signal. No assumption is therefore
#' required on the shape of the signals. The data are assumed to be normally distributed
#' (or at least follow a symmetric distribution) with a known constant variance. Working
#' pixelwise, the method can be time-consuming depending on the size of the data-array
#' but harnesses the power of multicore cpus.
#' 
#' @docType package
#' @author Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
#' @import parallel
NULL
