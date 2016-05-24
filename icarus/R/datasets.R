#' @name dataPop
#' @title Test population for Icarus.
#' @description This data set features a generated population of 50000 units.
#' 11 characteristics of interest for all units in population are featured. 
#' These characteristics of interest are variously correlated to one another.
#' A stratified random sampling (with a proportional allocation on variable Y3) 
#' of fixed size 1000 is selected. Among the 1000 units in the selected sample, only
#' 718 are respondant to the survey. These responding units are selected using a dummy
#' logit model.
#' @docType data
#' @usage dataPop
#' @format 1 column "ident" with unique id for all units.
#' 11 columns with various characteristics of interest for units in the population.
#' 1 column "weight", with sampling weights . Weights equal to zero means that the unit
#' is not selected in the sample.
#' 1 column "simul_nr" indicates the probability that each unit will respond to the survey.
#' 1 column "responding". For sampled units, indicates whether unit is respondant to survey (1)
#' or not (0). Variable is also equal to 0 for units not selected in sample
#' 50000 rows, 1 row per unit in the population.
#' @author Antoine Rebecq
#' @references Rebecq, A., & Merly-Alpa, T. Pourquoi minimiser la dispersion des 
#' poids en sondage. preprint.
NULL

#' @name poptest_calmar
#' @title Calibration on population test - made on Calmar2
#' @description This data set features calibration weights for the sample test of \code{\link{dataPop}}
#' (using margins tables \code{\link{table_margins_1}} and \code{\link{table_margins_2}}). Calibration is
#' is computed using the SAS Macro Calmar2, for test purposes.
#' @docType data
#' @usage poptest_calmar
#' @format 1000 rows, one per unit in the sample.
#' 1 column "ident", with a unique id for every unit in the sample
#' 3 methods of calibration are used (linear, raking, and logit with bounds LO=0.2 and UP=1.3) for
#' two different margins tables \code{\link{table_margins_1}} and \code{\link{table_margins_2}}, which results in 
#' 6 columns of weights.
#' @author Antoine Rebecq
#' @references Le Guennec, J., and Sautory, O. (2002). Calmar 2: Une nouvelle version 
#' de la macro calmar de redressement d'echantillon par calage. 
#' Journees de Methodologie Statistique, Paris. INSEE.
NULL

#' @name poptest_calmar_nr
#' @title Calibration with nonresponse on population test - made on Calmar2
#' @description This data set features calibration weights for the sample test of \code{\link{dataPop}}
#' (using margins tables \code{\link{table_margins_1}} and \code{\link{table_margins_2}}). Calibration is
#' is computed using the SAS Macro Calmar2, for test purposes. Only the 718 responding units are
#' taken into account.
#' @docType data
#' @usage poptest_calmar_nr
#' @format 718 rows, one per unit in the sample.
#' 1 column "ident", with a unique id for every unit in the sample
#' 3 methods of calibration are used (linear, raking, and logit with bounds LO=0.1 and UP=2.0 and parameter ECHELLE=0) for
#' two different margins tables \code{\link{table_margins_1}} and \code{\link{table_margins_2}}, which results in 
#' 6 columns of weights.
#' @author Antoine Rebecq
#' @references Le Guennec, J., and Sautory, O. (2002). Calmar 2: Une nouvelle version 
#' de la macro calmar de redressement d'echantillon par calage. 
#' Journees de Methodologie Statistique, Paris. INSEE.
NULL

#' @name table_margins_1
#' @title Margins for calibration of test population
#' @description This table features calibration margins for the sample of the test population
#' of \code{\link{dataPop}}
#' @docType data
#' @usage table_margins_1
#' @format A margins table written in the Icarus format.
#' @author Antoine Rebecq
NULL


#' @name table_margins_2
#' @title Margins for calibration of test population
#' @description This table features calibration margins for the sample of the test population
#' of \code{\link{dataPop}}. Margins for categorical variables are entered in percentages.
#' @docType data
#' @usage table_margins_2
#' @format A margins table written in the Icarus format.
#' @author Antoine Rebecq
NULL

#' @name data_ex2
#' @title A small example sample for calibration with Icarus
#' @description This table features a samples of 15 units (drawn from a population of size 300), 
#' used in a small survey to determine how frequently the employees of a firm go the movies (column "cinema").
#' Some auxiliary variables are given, which allows the use of calibration to improve
#' estimates. Margins for these auxiliary variables are known:
#' categ: 80 (modality 1) ; 90 (modality 2) ; 60 (modality 3)
#' sexe: 140 (modality 1) ; 90 (modality 2)
#' service: 100 (modality 1) ; 130 (modality 2)
#' salaire : 470000
#' @docType data
#' @usage data_ex2
#' @format 15 rows, one per unit in sample.
#' 1 column "id", unique id for each unit.
#' 4 columns of auxiliary variables ("service", "categ", "sexe", "salaire").
#' 1 column "cinema" - the variable of interest
#' 1 column "weight" - the Horvitz-Thompson weights
#' @author Antoine Rebecq
NULL

#' @name calWeights_ex2
#' @title Calibration weights for \code{\link{data_ex2}}
#' @description Calibration weights computed with Calmar2 for the small example \code{\link{data_ex2}}.
#' @docType data
#' @usage calWeights_ex2
#' @format 1 column "id", unique id for each of the 15 units in sample.
#' 3 columns with calibration weights using 3 different methods (linear, raking, and logit
#' with bounds LO=0.4, UP=2.2)
#' @author Antoine Rebecq
NULL