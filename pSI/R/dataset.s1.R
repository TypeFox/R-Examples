#'@title
#'Supplementary Data Set 1
#' 
#'@description
#'Differentially expressed genes (up-regulated & down-regulated) in cortices of autsitic vs control subjects from human transcriptomic data.
#'NOTE:Supplementary data (human & mouse expression sets, calculated pSI datasets, etc.) can be found in \code{pSI.data} package located at the following URL:
#'\url{http://genetics.wustl.edu/jdlab/psi_package/}

#'@details
#'\code{dataset.s1} is a list which contains two character vectors of differentially expressed genes. Differential expression 
#'was assessed using the SAM package (Significance Analysis of Microarrays), for FDR <.05 and fold change >1.3. 
#'\itemize{
#'  \item{\code{dataset.s1$up.regulated}}\cr{Up-regulated genes (N=234)}
#'  \item{\code{dataset.s1$down.regulated}}\cr{ - Down-regulated genes (N=229)}
#'}
#'
#'@examples
#'data(dataset.s1)
#'
#'@docType data
#'@keywords datasets
#'@format List containing two character vectors
#'@name dataset.s1
#'@source
#'Voineagu I, Wang X, Johnston P, Lowe JK, Tian Y, Horvath S, et al. (2011): Transcriptomic 
#'analysis of autistic brain reveals convergent molecular pathology. Nature. 474:380-384.
NULL