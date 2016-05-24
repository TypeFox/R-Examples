##' A package to perform robust quantitative traits association testing of copy number variants. It provides two methods for association study: first, the observed probe intensity measurement can be directly used to detect the association of CNV with phenotype of interest. Second, the most probable copy number is estimated with the proposed likelihood and the association of the most probable copy number with phenotype is tested. Also, it can be used to determine the optimal clustering number and clustering assignment for each individuals. This method can be applied to both the independent and correlated population.
##'
##' \tabular{ll}{
##' Package: \tab PedCNV\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1\cr
##' Date: \tab 2013-09-03\cr
##' License: \tab MIT \cr
##' Main functions:
##' \tab AssoTestProc \cr
##' \tab ClusProc \cr
##' \tab STE \cr
##' \tab STIM \cr
##' \tab print.asso\cr
##' \tab print.clus\cr
##' \tab plot.clus\cr
##' }
##'
##' @name PedCNV-package
##' @docType package
##' @title CNV association implementation
##' @author Meiling Liu, Sungho Won and Weicheng Zhu
##' @useDynLib PedCNV
##' @references On the association analysis of CNV data: fast and efficient method with family-based samples
NULL

##' The simulated intensity measurements. The order of the row in this file must consistent with the second column in FAM file.
##' 
##' @title CNV simulated intensity measurements
##' @name signal
##' @docType data
##' @author Meiling Liu
NULL

##' The simulated environmental file which contains the possible environmental variables. The order of the row in this file must consistent with the second column in FAM file.
##' 
##' @title CNV simulated environmental variables
##' @name envirX
##' @docType data
##' @author Meiling Liu
NULL

##' The simulated FAM file. The first six columns of FAM file are mandatory: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown) and Phenotype.
##' 
##' @title CNV simulated data
##' @name fam
##' @docType data
##' @author Meiling Liu
NULL


##' Empirical/kinship correlation matrix between individuals. This correlation matrix can be calculated based on the familial relationship between individuals or large-scale SNP data by omic data analysis toolkit FQLS. The free software FQLS can be downloaded from \url{http://biostat.cau.ac.kr/fqls/}. If correlation matrix is estimated with the large-scale SNP data, the proposed method becomes robust under the presence of population substructure.
##' 
##' @name phi
##' @title Empirical correlation matrix
##' @docType data
##' @author Meiling Liu
##' @references FQLS \url{http://biostat.cau.ac.kr/fqls/}
##' @examples
##' data(phi)
NULL


