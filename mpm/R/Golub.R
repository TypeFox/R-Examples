#' Golub (1999) Data
#' Golub et al. (1999) data on gene expression profiles of 38 patients
#' suffering from acute leukemia and a validation sample of 34 patients.
#' 
#' The original data of Golub et al. (1999) were preprocessed as follows: genes
#' that were called 'absent' in all samples were removed from the data sets,
#' since these measurements are considered unreliable by the manufacturer of
#' the technology.  Negative measurements in the data were set to 1.
#' 
#' The resulting data frame contains 5327 genes of the 6817 originally reported
#' by Golub et al. (1999).
#' 
#' @name Golub
#' @aliases Golub Golub.grp
#' @docType data
#' @format The expression data are available in data frame \code{Golub} with
#'   5327 observations on the following 73 variables.  
#' \describe{
#'   \item{list("Gene")}{a character vector with gene identifiers}  
#'   \item{list("1")}{gene expression data for sample 1}
#'   \item{list("2")}{gene expression data for sample 2} 
#'   \item{list("3")}{gene expression data for sample 3} 
#'   \item{list("4")}{gene expression data for sample 4} 
#'   \item{list("5")}{gene expression data for sample 5}
#'   \item{list("6")}{gene expression data for sample 6} 
#'   \item{list("7")}{gene expression data for sample 7} 
#'   \item{list("8")}{gene expression data for sample 8} 
#'   \item{list("9")}{gene expression data for sample 9}
#'   \item{list("10")}{gene expression data for sample 10}
#'   \item{list("11")}{gene expression data for sample 11}
#'   \item{list("12")}{gene expression data for sample 12}
#'   \item{list("13")}{gene expression data for sample 13}
#'   \item{list("14")}{gene expression data for sample 14}
#'   \item{list("15")}{gene expression data for sample 15}
#'   \item{list("16")}{gene expression data for sample 16}
#'   \item{list("17")}{gene expression data for sample 17}
#'   \item{list("18")}{gene expression data for sample 18}
#'   \item{list("19")}{gene expression data for sample 19}
#'   \item{list("20")}{gene expression data for sample 20}
#'   \item{list("21")}{gene expression data for sample 21}
#'   \item{list("22")}{gene expression data for sample 22}
#'   \item{list("23")}{gene expression data for sample 23}
#'   \item{list("24")}{gene expression data for sample 24}
#'   \item{list("25")}{gene expression data for sample 25}
#'   \item{list("26")}{gene expression data for sample 26}
#'   \item{list("27")}{gene expression data for sample 27}
#'   \item{list("34")}{gene expression data for sample 34}
#'   \item{list("35")}{gene expression data for sample 35}
#'   \item{list("36")}{gene expression data for sample 36}
#'   \item{list("37")}{gene expression data for sample 37}
#'   \item{list("38")}{gene expression data for sample 38}
#'   \item{list("28")}{gene expression data for sample 28}
#'   \item{list("29")}{gene expression data for sample 29}
#'   \item{list("30")}{gene expression data for sample 30}
#'   \item{list("31")}{gene expression data for sample 31}
#'   \item{list("32")}{gene expression data for sample 32}
#'   \item{list("33")}{gene expression data for sample 33}
#'   \item{list("39")}{gene expression data for sample 39}
#'   \item{list("40")}{gene expression data for sample 40}
#'   \item{list("42")}{gene expression data for sample 42}
#'   \item{list("47")}{gene expression data for sample 47}
#'   \item{list("48")}{gene expression data for sample 48}
#'   \item{list("49")}{gene expression data for sample 49}
#'   \item{list("41")}{gene expression data for sample 41}
#'   \item{list("43")}{gene expression data for sample 43}
#'   \item{list("44")}{gene expression data for sample 44}
#'   \item{list("45")}{gene expression data for sample 45}
#'   \item{list("46")}{gene expression data for sample 46}
#'   \item{list("70")}{gene expression data for sample 70}
#'   \item{list("71")}{gene expression data for sample 71}
#'   \item{list("72")}{gene expression data for sample 72}
#'   \item{list("68")}{gene expression data for sample 68}
#'   \item{list("69")}{gene expression data for sample 69}
#'   \item{list("67")}{gene expression data for sample 67}
#'   \item{list("55")}{gene expression data for sample 55}
#'   \item{list("56")}{gene expression data for sample 56}
#'   \item{list("59")}{gene expression data for sample 59}
#'   \item{list("52")}{gene expression data for sample 52}
#'   \item{list("53")}{gene expression data for sample 53}
#'   \item{list("51")}{gene expression data for sample 51}
#'   \item{list("50")}{gene expression data for sample 50}
#'   \item{list("54")}{gene expression data for sample 54}
#'   \item{list("57")}{gene expression data for sample 57}
#'   \item{list("58")}{gene expression data for sample 58}
#'   \item{list("60")}{gene expression data for sample 60}
#'   \item{list("61")}{gene expression data for sample 61}
#'   \item{list("65")}{gene expression data for sample 65}
#'   \item{list("66")}{gene expression data for sample 66}
#'   \item{list("63")}{gene expression data for sample 63}
#'   \item{list("64")}{gene expression data for sample 64}
#'   \item{list("62")}{gene expression data for sample 62} }
#' 
#' The classes are in a separate numeric vector \code{Golub.grp} with values
#'   \code{1} for the 38 ALL B-Cell samples, \code{2} for the 9 ALL T-Cell
#'   samples and \code{3} for the 25 AML samples.
#' @note Luc Wouters et al. (2003), p. 1134 contains a typo concerning the
#'   sample sizes of AML- and ALL-type and erroneously reported
#' @references Luc Wouters et al. (2003). Graphical Exploration of Gene
#'   Expression Data: A Comparative Study of Three Multivariate Methods,
#'   Biometrics, 59, 1131-1139.
#' @source Golub, T. R., Slonim, D. K., Tamayo, P., et al. (1999). Molecular
#'   classification of cancer: Class discovery and class prediction by gene
#'   expression monitoring. Science 286, 531 -- 537.
#' @keywords datasets
NULL



