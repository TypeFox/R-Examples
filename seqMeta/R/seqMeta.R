#' seqMeta: Meta-Analysis of Region-Based Tests of Rare DNA Variants
#'
#' Computes necessary information to meta analyze region-based tests for rare genetic variants (e.g. SKAT, T1) in individual studies, and performs meta analysis.
#'
#'
#' To learn more about seqMeta, start with the vignettes:
#' \code{browseVignettes(package = "seqMeta")}
#'
#' @name seqMeta
#' @docType package
#' @useDynLib seqMeta
#' @import methods Matrix CompQuadForm survival
#' @importFrom coxme lmekin
#' @importFrom stats na.action na.omit dbeta dchisq qchisq pchisq pnorm uniroot 
#' integrate gaussian model.matrix lm glm residuals model.frame var cov coef
#' @importFrom utils capture.output setTxtProgressBar txtProgressBar
#' @export prepCondScores skatMeta burdenMeta singlesnpMeta skatOMeta
#' @exportMethod c.seqMeta c.seqMeta
NULL


# 
# #' Illumina HumanExome BeadChip SNP Information file
# #' 
# #' Contains standard Names and associated genes for the Illumina HumanExome
# #' BeadChip
# #' 
# #' @format   A data frame with 247504 observations on the following 2 variables.
# #'   \describe{ 
# #'   \item{\code{Chr}}{Chromosome} 
# #'   \item{\code{Name}}{factor: illumina variant names} 
# #'   \item{\code{SKATgene}}{factor: gene names}
# #'   }
# #'   
# #' @details There are several non-exonic SNPs included. For these SNPs the
# #'   `gene` name is the same as the illumina variant name.
# #'   
# #' @references Grove ML, Cochran BJ, Haritunians T, Bis JC, Taylor KD, Hansen M,
# #'   O'Donnell CJ, Rotter JI, Boerwinkle E, CHARGE Exome Chip Genotyping
# #'   Committee; Best practices and joint calling of the Illumina HumanExome
# #'   BeadChip: the CHARGE consortium; (Abstract/Program #1445W). Presented at
# #'   the 62nd Annual Meeting of The American Society of Human Genetics (ASHG),
# #'   November 7, 2012, San Francisco, CA.
# #'   
# #'   Grove ML, Yu B, Cochran BJ, Haritunians T, Bis JC, Taylor KD, Hansen M,
# #'   Borecki I, Cupples LA, Fornage M, Gudnason V, Harris T, Katherisan S,
# #'   Kraaij R, Launer LJ, Levy D, Liu Y, Mosley T, Peloso GM, Psaty BM, Rich SS,
# #'   Rivadeneira F, Siscovick DS, Smith AV, Uitterlinden A, van Duijn CM, Wilson
# #'   JG, O'Donnell CJ, Rotter JI, Boerwinkle E. Best practices and joint calling
# #'   of the Illumina HumanExome BeadChip: the CHARGE consortium. PLoS One
# #'   [submitted]
# #'
# #' @examples 
# #' data(SNPInfo)
# #' ##summary of the data set:
# #' summary(as.numeric(table(SNPInfo$SKATgene)))
# #' hist(table(SNPInfo$SKATgene),nclass = 300,xlim=c(0,50), main = "SNPs per gene", xlab ="#SNPs", ylab = "#Genes")
# #' 
# #' @keywords datasets
# "SNPInfo"