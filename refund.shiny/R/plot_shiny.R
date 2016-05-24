#' Interactive Plotting for Functional Data
#' 
#' Function for interactive plotting of functional data analysis results. 
#' 
#' This package builds on the \code{refund} package: tools in \code{refund} are used to 
#' conduct analyses and functions in this package create interactive visualizations of the results
#' of those analyses. There are four major categories of analyses that can be viewed:
#' \enumerate{
#' \item{Functional principal components analyses implemented by \code{\link{fpca.sc}}, \code{\link{fpca.face}}, 
#' \code{\link{fpca.ssvd}}, and \code{\link{fpca2s}}. Plots show the mean +/- 2SD times each FPC; scree plots;
#' linear combinations of score values and FPCs; reconstructions for each subject; and score scatterplots.}
#' \item{Function-on-scalar regression analyses implemented by \code{\link{bayes_fosr}}. Plots show the raw data
#' colored by covariate values; fitted values depending on covariates; coefficient functions; and residuals.}
#' \item{Multilevel functional principal components analyses implemented by \code{\link{mfpca.sc}}.  Plots show the 
#' mean +/- 2SD times each FPC; scree plots; linear combinations of score values and FPCs; 
#' reconstructions for each subject; and score scatterplots for levels 1 and 2.}
#' #' \item{Longitudinal functional principal components analyses}
#' }
#' 
#' 
#' 
#' @title plot_shiny The generic function for interactive plots of functional data analyses
#' @param obj object to be plotted. Currently, allowed data types are \code{fpca} \code{mfpca} \code{lfpca} and \code{fosr}.
#' @param ... additional arguments passed to plotting functions
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}, 
#' Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @seealso \code{\link{plot_shiny.fpca}}, \code{\link{plot_shiny.mfpca}},  \code{\link{plot_shiny.fosr}}
#' @import refund
#' @export plot_shiny
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(refund)
#' library(dplyr)
#' 
#' ##### FPCA Example on real data #####
#' 
#' data(cd4)
#' SC = fpca.sc(cd4)
#' plot_shiny(SC)
#' 
#' ##### FPCA Examples on simulated data #####
#' 
#' set.seed(2678695)
#' n = 101
#' m = 101
#' s1 = 20
#' s2 = 10
#' s = 4
#' t = seq(-1, 1, l=m)
#' v1 = t + sin(pi*t)
#' v2 = cos(3*pi*t)
#' V = cbind(v1/sqrt(sum(v1^2)), v2/sqrt(sum(v2^2)))
#' U = matrix(rnorm(n*2), n, 2)
#' D = diag(c(s1^2, s2^2))
#' eps = matrix(rnorm(m*n, sd=s), n, m)
#' Y = U%*%D%*%t(V) + eps
#'
#' SC = fpca.sc(Y)
#' FACE = fpca.face(Y)
#' SSVD = fpca.ssvd(Y, verbose=FALSE)
#' S = fpca2s(Y)
#'  
#' plot_shiny(SC)
#' plot_shiny(FACE)
#' plot_shiny(SSVD)
#' plot_shiny(S)
#' 
#' #' ##### MFPCA Example #####
#' 
#' data(DTI)
#' Y = DTI$cca
#' id = DTI$ID
#' 
#' mfpca.dti = mfpca.sc(Y=Y, id = id, twoway = FALSE)
#' plot_shiny(mfpca.dti)
#' 
#' ##### FoSR Example #####
#' 
#' data(DTI)
#' DTI = DTI[complete.cases(DTI),]
#' fit.fosr = bayes_fosr(cca ~ pasat + sex, data = DTI)
#' plot_shiny(fit.fosr)
#' 
#' ##### FoSR Example with outliers #####
#' 
#' DTI$cca[1,] = DTI$cca[1,] + .4
#' DTI$cca[2,] = DTI$cca[2,] + .4
#' 
#' fosr.dti2 = bayes_fosr(cca ~ pasat + sex, data = DTI)
#' plot_shiny(fosr.dti2)
#' 
#' ##### Longitudinal FoSR Examples #####
#' 
#' data(DTI2)
#' class(DTI2$cca) = class(DTI2$cca)[-1]
#' DTI2 = subset(DTI2, select = c(cca, id, pasat))
#' DTI2 = DTI2[complete.cases(DTI2),]
#' 
#' fosr.dti3 = bayes_fosr(cca ~ pasat + re(id), data = DTI2, Kt = 10, Kp = 4, cov.method = "FPCA")
#' plot_shiny(fosr.dti3)
#' plot_shiny(fosr.dti3$fpca.obj)
#' 
#' ##### LFPCA Example on real data #####
#' 
#' data(DTI)
#' MS <- subset(DTI, case ==1)  # subset data with multiple sclerosis (MS) case
#'
#' index.na <- which(is.na(MS$cca))  
#' Y <- MS$cca; Y[index.na] <- fpca.sc(Y)$Yhat[index.na]; sum(is.na(Y))
#' id <- MS$ID 
#' visit.index <- MS$visit 
#' visit.time <- MS$visit.time/max(MS$visit.time)
#'
#' lfpca.dti1 <- fpca.lfda(Y = Y, subject.index = id,  
#'                        visit.index = visit.index, obsT = visit.time, 
#'                        LongiModel.method = 'lme',
#'                        mFPCA.pve = 0.95)
#' plot_shiny(lfpca.dti1)
#' 
#' lfpca.dti2 <- fpca.lfda(Y = Y, subject.index = id,  
#'                        visit.index = visit.index, obsT = visit.time, 
#'                        LongiModel.method = 'fpca.sc',
#'                        mFPCA.pve = 0.80, sFPCA.pve = 0.80)
#' plot_shiny(lfpca.dti2)
#' 
#' 
#' }
#' 
plot_shiny <- function(obj, ...){
  UseMethod("plot_shiny")
}