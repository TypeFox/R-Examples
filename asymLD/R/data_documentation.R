#' Pairwise haplotype frequencies for 7 HLA loci
#'
#' HLA haplotype frequencies for 21 pairs of HLA loci in a set of 300 controls
#' from a study of myopericarditis incidence following smallpox vaccination. 
#'
#'
#'
#' @format A data frame with 3063 observations on the following 5 variables.
#'  \describe{
#' \item{\code{haplo.freq}}{a numeric vector}
#' \item{\code{locus1}}{a character vector}
#' \item{\code{locus2}}{a character vector}
#' \item{\code{allele1}}{a character vector}
#' \item{\code{allele2}}{a character vector}
#' }
#' 
#' @name hla.freqs
#' @docType data
#' @source Wilson, C., 2010 Identifying polymorphisms associated with risk for the development of myopericarditis following smallpox vaccine. The Immunology Database and Analysis Portal (ImmPort), Study #26.
#' 
#' 
#' @references https://immport.niaid.nih.gov/immportWeb/clinical/study/displayStudyDetails.do?itemList=SDY26
#' 
#' 
#' @keywords data
#' @usage data(hla.freqs)
#' 
NULL


#' Pairwise haplotype frequencies for 4 SNP loci
#'
#' HLA haplotype frequencies for 6 pairs of SNP loci from the de Bakker et al. 2006 data for 90 unrelated individuals with European ancestry (CEU) from the Centre d'Etude du Polymorphisme Humain (CEPH)
#' collection obtained from the Tagger/MHC webpage.
#'
#' @details For bi-allelic SNP data the ALD measures are symmetric and equivalent to the r correlation measure of LD.
#'
#' @format A data frame with 20 observations on the following 5 variables. 
#'   \describe{
#'     \item{\code{locus1}}{a character vector}
#'     \item{\code{locus2}}{a character vector}
#'     \item{\code{allele1}}{a character vector}
#'     \item{\code{allele2}}{a character vector}
#'     \item{\code{haplo.freq}}{a numeric vector}
#'   }
#' 
#' @name snp.freqs
#' @docType data
#' @source de Bakker, P. I., G. McVean, P. C. Sabeti, M. M. Miretti, T. Green et al., 2006 A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC. Nat.Genet. 38: 1166-1172.
#' 
#' 
#' @references http://www.broadinstitute.org/mpg/tagger/mhc.html
#' 
#' 
#' @keywords data
#' @usage data(snp.freqs)
#' 
#' 
NULL
