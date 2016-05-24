#' An R package for computing gene and pathway p-values using the Adaptive Rank Truncated test (ARTP). 
#' This package can be used to analyze pathways/genes based on a genetic association study, 
#' with a binary case-control outcome or a survival outcome. This package is an extension of the ARTP method developped
#' by Kai Yu (Genet Epidemiol. 2009) for gene- and pathway-environment interaction analysis.
#' 
#' The statistical significance of the pathway-level test statistics is evaluated using a highly 
#' efficient permutation algorithm that remains computationally feasible irrespective of the size 
#' of the pathway and complexity of the underlying test statistics for summarizing SNP- and gene-level
#'  associations. The function  \code{\link{ARTP.GE}} is used to compute gene and pathway p-values provided 
#'  that the observed and permutation p-values for each SNP already exist in files. The input files required for  \code{\link{ARTP.GE}}
#'  could be obtained by calling the function  \code{\link{permutation.snp}} and the function  \code{\link{compute.p.snp.obs}}.
#' @name PIGE-package
#' @aliases PIGE
#' @docType package
#' @title Gene and pathway p-values using the Adaptive Rank Truncated Product test
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}
#' @keywords package
#' @references  Yu K, Li Q, Berger AW, Pfeiffer R, Rosenberg P, Caporaso N, Kraft P, Chatterjee N (2009). Pathway analysis by adaptive combination of P-values. Genet Epidemiol 33:700-709.
#' @seealso  \code{\link{ARTP.GE}},  \code{\link{permutation.snp}},  \code{\link{compute.p.snp.obs}}
NULL
#' Sample data of a case-control study  
#' 
#' This fictive data set contains one environmental variable (var_int), 4 fixed covariates (cov) and SNP variables (coded 0,1,2).   
#' @name data.pige
#' @docType data
#' @format A dataframe with 1000 rows and 134 columns
NULL
#' Fictive list for the case-control study example containing the names of the snp for each gene included in the studied pathways
#' @name list.gene.snp
#' @docType data
#' @format A list containing the names of the SNPs belonging to each gene analysed. 
NULL
#' Fictive data frame containing the names of the gene (first column) and the names of the SNP belonging to the SNPS (second column)
#' @name snp.gene.select
#' @docType data
#' @format A data frame of 130 entries corresponding of the gene's name of each SNP.
NULL
#' Fictive data frame corresponding to the gene x pathway analysis for the example of case-control study.
#' @name data.pathway
#' @docType data
#' @format A data frame of 0 and 1 indicating in which pathway (column) 
#' each gene (rows) is included.
NULL
#' Fictive data set for sample data of a survival analysis.  
#' 
#' This fictive data set contains a time to event (TIME), an indicator variable (named EVENT) coded 0/1 (censored/event), one environmental variable (var_int), and SNP variables.   
#' @name data.surv
#' @docType data
#' @format A dataframe with 66 rows and 887 columns
NULL
#' Fictive list (for survival example) containing the names of the snp for each gene included in the studied pathways
#' @name list.gene.snp.surv
#' @docType data
#' @format A list containing the names of the SNPs belonging to each gene analysed. 
NULL
#' Fictive data frame corresponding to the gene x pathway analysis for the example of survival analysis.
#' @name data.pathway.surv
#' @docType data
#' @format A data frame of 0 and 1 indicating in which pathway (column) 
#' each gene (rows) is included.
NULL
