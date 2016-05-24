#' Bayesian Meta-Analysis of Diagnostic Test Data
#'
#' Bayesian meta-analysis of diagnostic test data based on a scale mixtures
#' bivariate random-effects model.
#' This package was developed with the aim of simplifying the use of meta-analysis
#' models that up to now have demanded great statistical expertise in Bayesian meta-analysis.
#' The package implements a series of innovative statistical techniques including:
#' the BSROC (Bayesian Summary ROC) curve, the BAUC (Bayesian AUC), predictive surfaces,
#' the use of prior distributions  that avoid boundary estimation problems of component
#' of variance and correlation parameters, analysis of conflict of evidence and robust
#' estimation of model parameters. In addition, the package comes with several published
#' examples of meta-analysis that can be used for illustration or further research in
#' this area.
#'
#' \tabular{ll}{ Package: \tab bamdit\cr Type: \tab Package\cr Version: \tab
#' 2.0\cr Date: \tab 2015-06-15 \cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#' yes\cr }
#'
#' @name bamdit-package
#' @aliases bamdit-package bamdit
#' @docType package
#' @author Dr. Pablo Emilio Verde \email{pabloemilio.verde@@hhu.de}
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. Statistics in Medicine. 29, 3088-3102.
#' @keywords Bayesian modeling
#' JAGS
#' Meta-Analysis
#' Bivariate random-effects
#' Scale mixtures
NULL

#' Diagnosis of appendicities with computer tomography scans
#'
#' This data frame corresponds to 51 clinical studies reporting the accuracy of
#' computer tomography (CT) scans for the diagnosis of appendicities.
#'
#'
#' @name ct
#' @docType data
#' @format A matrix with 51 rows and 16 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{   tp   }{true positives}
#' \item{   n1   }{number of patients with disease}
#' \item{   fp   }{false positives}
#' \item{   n2   }{number of patients without disease}
#' \item{   country   }{Country: EU = 1, others/USA = 2}
#' \item{   hosp   }{Type of hospital: 1 = university, 2 = others}
#' \item{   inclus   }{Inclusion criteria: 1 = Suspected, 2 = appendectomy}
#' \item{    indfind     }{Other CT findings included: 1 = no, 2 = yes}
#' \item{    design     }{Study design: 1 = prospective, 2 = retrospective}
#' \item{    contr     }{Contrast medium: 1 = no, 2 = yes}
#' \item{    localis     }{Localisation: 1 = one area, 2 = more than one area}
#' \item{    child     }{Children included: 1 = no, 2 = yes}
#' \item{    fup.na     }{Followup: 0 = no, 1 = yes }
#' \item{    refer.na     }{Valid reference: 0 = no, 1= yes }
#' \item{    sample.na     }{Sample: 0 = selected, 1= consecutive/random }
#' \item{    gender.na     }{Gender, female: 0 = less than 50\%; 1 = more than 50\% } }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source Ohmann C, Verde PE, Gilbers T, Franke C, Fuerst G, Sauerland S,
#' Boehner H. (2006) Systematic review of CT investigation in suspected acute
#' appendicitis. \emph{Final Report; Coordination Centre for Clinical Trials,
#' Heinrich-Heine University}. Moorenstr. 5, D-40225 Duesseldorf Germany.
#' @keywords datasets
NULL


#' Ectopic pregnancy vs. all other pregnancies data
#'
#' Ectopic pregnancy vs. all other pregnancies data
#' Table III Mol et al. 1998
#'
#' @name ep
#' @docType data
#' @format A matrix with 21 rows. and 8 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{    tp   }{true positives}
#' \item{    n1   }{number of patients with disease}
#' \item{    fp   }{false positives}
#' \item{    n2   }{number of patients without disease}
#' \item{    d1   }{Prospective vs retrospective}
#' \item{    d2   }{Cohort vs case-control}
#' \item{    d3   }{Consecutive sampling patients series vs non-consecutive} }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source Table III Mol et al. 1998
#' @keywords datasets
NULL


#' Radiological evaluation of lymph node metastases in patients with cervical cancer: a
#' meta-analysis.
#'
#' This data frame summarizes the tables 1-3 of Scheidler et al. 1997.
#'
#'
#' @name scheidler
#' @docType data
#' @format A matrix with 46 rows and 7 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{    tp   }{true positives}
#' \item{    n1   }{number of patients with disease}
#' \item{    fp   }{false positives}
#' \item{    n2   }{number of patients without disease}
#' \item{    author   }{first author of the study}
#' \item{    year   }{publication date}
#' \item{    test   )}{test method used in the study} }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source The data were obtainded from
#'
#' Scheidler J, Hricak H, Yu KK, Subak L, Segal MR. (1997) Radiological
#' evaluation of lymph node metastases in patients with cervical cancer: a
#' meta-analysis. \emph{The Journal of the American Medical Association};
#' \bold{278}:1096-1101.
#' @keywords datasets
NULL


#' Accuracy of Positron Emission Tomography for Diagnosis of Pulmonary
#' Nodules and Mass Lesions
#'
#' Data from a Meta-Analysis of Studies on
#' Quality of Studies of FDG-PET for Diagnosis of SPNs and Mass Lesions
#'
#' @name gould
#' @docType data
#' @format A matrix with 31 rows and 6 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{    tp   }{true positives}
#' \item{    n1   }{number of patients with disease}
#' \item{    fp   }{false positives}
#' \item{    n2   }{number of patients without disease}
#' \item{    author   }{first author of the study}
#' \item{    year   }{publication
#'  date} }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source The data were obtainded from
#'
#' Gould MK, Maclean CC, Kuschner WG, Rydzak CE, Owens Dk. (2001) Accuracy of
#' positron emission tomography for diagnosis of pulmonary nodules and mass
#' lesions: a meta-analysis. \emph{The Journal of the American Medical
#' Association};\bold{285}:914-24.
#' @keywords datasets
NULL

#' Tumor markers in the diagnosis of primary bladder cancer.
#'
#' Outcome of individual studies evaluating urine markers
#'
#'
#' @name glas
#' @docType data
#' @format A matrix with 46 rows and 7 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{    tp   }{true positives}
#' \item{    n1   }{number of patients with disease}
#' \item{    fp   }{false positives}
#' \item{    n2   }{number of patients without disease}
#' \item{    author   }{first author of the study}
#' \item{    cutoff   }{cutoff in U/ml}
#'  \item{    marker   }{test method used in the study} }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source The data were obtainded from
#'
#' Glas AS, Roos D, Deutekom M, Zwindermann AH, Bossuyt PM, Kurth KH. (2003)
#' Tumor markers in the diagnosis of primary bladder cancer. A systematic
#' review. \emph{Journal of Urology}; \bold{169}:1975-82.
#' @keywords datasets
NULL

#' Diagnosis of Intravascular Device-Related Bloodstream Infection
#'
#' Outcome of individual studies evaluating intravascular device-related
#' bloodstream infection
#'
#'
#' @name safdar05
#' @docType data
#' @format A matrix with 78 rows and 8 columns. Each row represents study
#' results, the columns are:
#' \describe{
#' \item{    tp   }{true positives}
#' \item{    n1   }{number of patients with disease}
#' \item{    fp   }{false positives}
#' \item{    n2   }{number of patients without disease}
#' \item{    author   }{first author of the study}
#' \item{    year   }{publication date}
#' \item{    technique   }{diagnostic technique used in the study}
#' \item{    duration   }{duration of catheterization: short term or long term}
#' }
#' @references Verde P. E. (2010). Meta-analysis of diagnostic test data: A
#' bivariate Bayesian modeling approach. \emph{Statistics in Medicine}.
#' \bold{29}, 3088-3102.
#' @source The data were obtainded from
#'
#' Safdar N, Fine JP, Maki DG. (2005) Meta-analysis: methods for diagnosing
#' intravascular device-related bloodstream infection.
#'  \emph{Ann Intern Med.}; \bold{142}:451-66.
#' @keywords datasets
NULL

#' @docType data
#' @name mri
#' @title Diagnosis of lymph node metastasis with magnetic resonance imaging
#' @format A matrix with 10 rows and 4 columns. Each row represents study results, the columns are:
#' \describe{
#' \item{tp}{true positives}
#' \item{n1}{number of patients with disease}
#' \item{fp}{false positives}
#' \item{n2}{number of patients without disease} }
#' @source
#'  The data were obtainded from
#'
#'  Scheidler J, Hricak H, Yu KK, Subak L, Segal MR. (1997)
#'  Radiological evaluation of lymph node metastases in
#'  patients with cervical cancer: a meta-analysis. \emph{The
#'    Journal of the American Medical Association};
#'  \bold{278}:1096-1101.
#
#'@references
#'  Verde P. E. (2010). Meta-analysis of diagnostic test
#'  data: A bivariate Bayesian modeling approach.
#'  \emph{Statistics in Medicine}. \bold{29}, 3088-3102.
#'
#'@keywords datasets
NULL

