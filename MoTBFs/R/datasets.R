#' Data set Thyroid Disease (thyroid0387)
#' 
#' This data set if one of the several databases about Thyroid avalaible at the UCI repository. 
#' The task is to detect is a given patient is normal (1) or suffers from hyperthyroidism (2) 
#' or hypothyroidism (3) .
#' 
#' \describe{
#'   \item{Age}{Age of the patient (0.01--0.97). Continuous variable.}
#'   \item{Sex}{Sex of the patient, 0 (Male) 1 (Female). Binary variable. }
#'   \item{On_thyroxine}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Query_on_thyroxine}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{On_antithyroid_medication}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Sick}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Pregnant}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Thyroid_surgery}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{I131_treatment}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Query_hypothyroid}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Query_hyperthyroid}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Lithium}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Goitre}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Tumor}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Hypopituitary}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{Psych}{0 (FALSE) 1 (TRUE). Binary variable.}
#'   \item{TSH}{amount of TSH (0.0--0.53). Continuous variable.}
#'   \item{T3}{amount of T3 (0.0005--0.18). Continuous variable.}
#'   \item{TT4}{amount of TT4 (0.002--0.6). Continuous variable.}
#'   \item{T4U}{amount of T4U (0.017--0.233). Continuous variable.}
#'   \item{FTI}{amount of FTI (0.002--0.642). Continuous variable.}
#'   \item{Class}{1 (normal) 2 (hyperthyroidism) 3 (hypothyroidism). Class variable.}
#' }
#' @format A data frame with 7200 rows, 21 variables and the class.
#' @source \url{http://archive.ics.uci.edu/ml/datasets/Thyroid+Disease}
#' @name thyroid
NULL

#' Data set Ecoli: Protein Localization Sites
#' 
#' This data set contains information of Escherichia coli. It is a 
#' bacterium of the genus Escherichia that is commonly found in the 
#' lower intestine of warm-blooded organism.
#' 
#' \describe{
#'   \item{Sequence Name}{Accession number for the SWISS-PROT database.}
#'   \item{mcg}{McGeoch's method for signal sequence recognition.} 
#'   \item{gvh}{Von Heijne's method for signal sequence recognition.} 
#'   \item{lip}{Von Heijne's Signal Peptidase II consensus sequence score. Binary attribute.}
#'   \item{chg}{Presence of charge on N-terminus of predicted lipoproteins. Binary attribute.}
#'   \item{aac}{Score of discriminant analysis of the amino acid content of outer membrane and periplasmic proteins.}
#'   \item{alm1}{Score of the ALOM membrane spanning region prediction program.}
#'   \item{alm2}{Score of ALOM program after excluding putative cleavable signal regions from the sequence.}
#'   \item{Class}{Class variable. 8 possibles states.}
#' }
#' @format A data frame with 336 rows, 8 variables and the class.
#' @source \url{http://archive.ics.uci.edu/ml/datasets/Ecoli}
#' @name ecoli
NULL

