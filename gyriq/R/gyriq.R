#' Kinship-Adjusted Survival SNP-Set Analysis
#'
#' SNP-set association testing for censored phenotypes in the presence of 
#' intrafamilial correlation
#'
#' \tabular{ll}{
#' Package: \tab gyriq\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.2\cr
#' Date: \tab 2016-01-06\cr
#' License: \tab GPL (>= 2) \cr
#' }
#'
#' This variance-components test between a set of SNPs and a survival trait
#' is valid for both common and rare variants. A proportional hazards Cox model 
#' (written as a transformation model with censored data; Cheng et al., 1995)
#' is specified for the marginal distribution of the survival trait. The 
#' familial dependence is modelled via a Gaussian copula with a correlation 
#' matrix expressed in terms of the kinship matrix. The statistical procedure
#' has been described in full detail by Leclerc et al. (2015).
#' 
#' Censored values are treated as partially missing data and a multiple 
#' imputation procedure is employed to estimate vectors of residuals. These 
#' residuals and the SNPs in the genomic region under study are used to 
#' compute measures of phenotypic and genotypic similarity between pairs of 
#' subjects. The contribution to the score statistic is maximal when these 
#' measures are both high which corresponds to departure from the null 
#' hypothesis of no association between the set of SNPs and the survival 
#' outcome. The selection of the SNPs forming the SNP set can be based on 
#' biological information such as linkage disequilibrium (LD) blocks or rely on 
#' a sliding window method.
#'
#' The procedure is convenient for GWAS as the multiple imputation procedure
#' for the estimation of a completed vector of residuals has to be performed 
#' only once using the function \link{genComplResid}. A sliding window approach
#' can then be used to examine the evidence of association across the SNP set. 
#' In each run, the p-value is computed with the function \link{testGyriq}.
#' 
#' @name gyriq-package
#' @aliases gyriq-package gyriq
#' @docType package
#' @author Martin Leclerc <martin.leclerc.5@@ulaval.ca> and Lajmi Lakhal Chaieb 
#' <lakhal@@mat.ulaval.ca>
#' @references Cheng SC, Wei LJ, Ying Z. 1995. Analysis of transformation models
#' with censored data. Biometrika 82:835-845.
#' 
#' Leclerc M, The Consortium of Investigators of Modifiers of BRCA1/2, Simard J, 
#' Lakhal-Chaieb L. 2015. SNP set association testing for survival outcomes in 
#' the presence of intrafamilial correlation. Genetic Epidemiology 39:406-414.
#' 
#' Lin X, Zhou Q. 2015. coxKM: Cox Kernel Machine SNP-Set Association Test. R 
#' package version 0.3, URL http://www.hsph.harvard.edu/xlin/software.html#coxkm.
#'
#' Lin X, Cai T, Wu M, Zhou Q, Liu G, Christiani D and Lin X. 2011. Survival 
#' kernel machine SNP-set analysis for genome-wide association studies. Genetic 
#' Epidemiology 35:620-631.
#' 
#' Cai T, Tonini G and Lin X. 2011. Kernel machine approach to testing the 
#' significance of multiple genetic markers for risk prediction. Biometrics 
#' 67:975-986.
#' @keywords package
#' @examples
#' data(simGyriq)
#' for (i in seq_along(simGyriq)) assign(names(simGyriq)[i], simGyriq[[i]])
#' 
#' cr <- genComplResid(U, Delta, Phi, blkID, m=50, X)
#' testGyriq(cr$compResid, G, w, ker="LIN", asv=NULL, method="davies", 
#' starResid=NULL, bsw, tsw, pos)
NULL

#' Simulated SNP-set
#'
#' Simulated dataset of phenotypic, genotypic and kinship data.
#'
#' This dataset was generated under conditions described in Leclerc et al. 
#' (2015).
#' 
#' Samples of n = 600 individuals from 120 families were generated: 40 families 
#' of two parents and one child, 40 families of two parents and two children, 
#' and 40 families of three generations (two grand-parents, four parents, and 
#' two grandchildren). The coefficients of the block diagonal kinship matrix 
#' were fixed at their expected theoretical values. The number of biallelic SNPs
#' was set to s = 50. The minor allele frequencies were randomly 
#' sampled from Unif(0.001, 0.1). The genotypes of the 50 SNPs were 
#' simulated assuming a linkage disequilibrium corresponding to a squared 
#' correlation coefficient of r^2 = 0.5 between consecutive SNPs.
#' 
#' The two covariates follow Bernoulli(0.5) and Uniform(-0.2, 0.2) distributions 
#' respectively. The polygenic heritability parameter was fixed at 0.5. Each 
#' covariate parameter was set equal to 1 and the monotone increasing function 
#' of the transformation model with censored data (Cheng et al., 1995) was fixed
#' at H(t) = log(t) in order to generate the survival traits. The censoring rate 
#' was equal to 50\%. The weight of each SNP was defined as the density function
#' of the Beta (1, 25) evaluated at the corresponding minor allele frequency.
#'
#' The dataset includes simulated positions for the 50 SNPs, and the lower
#' and upper bounds of 4 sliding windows. Each window includes 10 SNPs, 
#' overlapping with the previous and subsequent windows. A vector of size B*n 
#' of permuted row indices is also included, where B=10,000. This is to 
#' be used to compute the p-value of the test following the standard or matching
#' moments permutation approach.
#' 
#' @name simGyriq
#' @docType data
#' @format A list containing the following elements:
#' \describe{
#' \item{U}{ 600x1 vector containing the survival times. \code{U = min(C, T)}
#' where \code{C} is the censoring time, and \code{T} the failure time }
#' \item{Delta}{ 600x1 vector containing the censoring indicator }
#' \item{Phi}{ 600x600 kinship matrix }
#' \item{blkID}{ 600x1 vector with entries identifying correlated groups of 
#' observations }
#' \item{X}{ 600x2 matrix of 2 covariates }
#' \item{G}{ 600x50 matrix containing the set of 50 SNPs }
#' \item{w}{ 50x1 vector of weights for the 50 SNPs }
#' \item{bsw}{ 4x1 vector containing the lower bounds of the 4 sliding windows 
#' considered for the SNP-set }
#' \item{tsw}{ 4x1 vector containing the upper bounds of the 4 sliding windows 
#' considered for the SNP-set }
#' \item{pos}{ 50x1 vector of SNP positions (used for the output only) }
#' \item{indResid}{ 10,000*600x1 vector of permuted row indices }
#' }
#' @keywords dataset
#' @references Cheng SC, Wei LJ, Ying Z. 1995. Analysis of transformation models
#' with censored data. Biometrika 82:835-845.
#' 
#' Leclerc M, The Consortium of Investigators of Modifiers of BRCA1/2, Simard J, 
#' Lakhal-Chaieb L. 2015. SNP set association testing for survival outcomes in 
#' the presence of intrafamilial correlation. Genetic Epidemiology 39:406-414.
#' @examples
#' data(simGyriq)
#' for (i in seq_along(simGyriq)) assign(names(simGyriq)[i], simGyriq[[i]])
#'
#' cr <- genComplResid(U, Delta, Phi, blkID, m=50, X)
#' testGyriq(cr$compResid, G, w, ker="LIN", asv=NULL, method="davies", 
#' starResid=NULL, bsw, tsw, pos)
NULL