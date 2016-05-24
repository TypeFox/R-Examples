#' Simulated Data from Buckley (2002)
#' 
#' In Buckley (2002) tissue residue curves for a Meningioma and a Breast Cancer
#' were simulated using the MMID4 model. Note, the model is described in detail
#' by Bassingthwaighte, J.B. \emph{et al.} (1984) and Kroll, K \emph{et al.}
#' (1996).  This model accounts for flow dispersion and heterogeneity, and
#' includes capillaries modeled as axially distributed blood-tissue exchange
#' units.  A plasma concentration-time curve, AKA arterial input function, was
#' simulated as an input to the model using measurements made by Fritz-Hansen
#' \emph{et al.} (1996).
#' 
#' 
#' @aliases buckley breast meningioma
#' @usage data("buckley")
#' @format Two lists are created (breast and meningioma) that contain the
#' simulated time curves and all associated kinetic parameter values.
#' @references Buckley, D.L. (2002) Uncertainty in the Analysis of Tracer
#' Kinetics Using Dynamic Contrast-Enhanced \eqn{T_1}{T1}-weighted MRI,
#' \emph{Magnetic Resonance in Medicine} \bold{47}, 601-606.
#' 
#' Bassingthwaighte, J.B. and Goresky, C.A. (1984) Modelling in the analysis of
#' solute and water exchange in the microvasculature. In: Renkin, E.M., Michel,
#' C.C. and Geiger, S.R., editors. Handbook of physiology. Section 2. The
#' cardiovascular system. Bethesda: American Physiological Society. p549-626.
#' 
#' Kroll, K., Wilke, N., Jerosch-Herold, M., Wang, Y., Zhang Y., Basche, R.J.
#' and Bassingthwaighte, J.B. (1996) Modelling regional myocardial flows from
#' residue functions of an intravascular indicator. \emph{Am J Physiol}
#' \bold{271}, H1643-H1655.
#' 
#' Fritz-Hansen, T., Rostrup, E., Larsson, H.B., Sondergaard, L., Ring, P. and
#' Hendriksen, O. (1996) Measurement of the arterial concentration Gd-DTPA
#' using MRI; a step toward quantitative perfusion imaging. \emph{Magn Reson
#' Med} \bold{36}, 347-357.
#' @source See below.
#' @keywords datasets
"buckley"

#' dcemri: A Package for Medical Image Analysis (S4 implementation)
#' 
#' A collection of routines and documentation that allows one to perform a 
#' quantitative analysis of dynamic contrast-enhanced or diffusion-weighted MRI 
#' data.  Medical imaging data should be organized using either the Analyze or 
#' NIfTI data formats.
#' 
#' Further information is available in the following vignettes: 
#' \tabular{ll}{ \code{dcemriS4} \tab dcemriS4(source, pdf)\cr }
#' 
#' @name dcemriS4-package
#' @aliases dcemriS4-package dcemriS4
#' @docType package
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}\cr Volker Schmid
#'   \email{volkerschmid@@users.sourceforge.net}\cr Andrew Thornton 
#'   \email{zeripath@@users.sourceforge.net}
#'   
#' @references Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang, 
#'   G.-Z.  (2006) Bayesian methods for pharmacokinetic models in dynamic 
#'   contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions on 
#'   Medical Imaging}, \bold{25} (12), 1627-1636.
#'   
#'   Schmid, V., Whitcher, B., Padhani, A.R. and G.-Z. Yang (2009) A 
#'   semi-parametric technique for the quantitative analysis of dynamic 
#'   contrast-enhanced MR images based on Bayesian P-splines, \emph{IEEE 
#'   Transactions on Medical Imaging}, \bold{28} (6), 789-798.
#' @keywords package
#' @examples
#' 
#' \dontrun{
#' demo(avg152T1)
#' demo(avg152T1LR)
#' demo(avg152T1RL)
#' demo(buckley)
#' demo(filtered_func_data)
#' demo(zstat1)
#' }
#' 
#' @import methods
#' @import oro.nifti
NULL
