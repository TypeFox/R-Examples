#' Copenhagen Stroke Study
#' 
#' This data set contains a subset of the data from the Copenhagen stroke
#' study.
#'
#' @name cost
#' @docType data
#' @format This data frame contains the observations of 518 stroke patients :
#' \describe{ \item{age}{Age of the patients in years.} \item{sex}{A factor
#' with two levels \code{female} and \code{male}.} \item{hypTen}{Hypertension,
#' a factor with two levels \code{no} and \code{yes}.} \item{ihd}{History of
#' ischemic heart disease at admission, a factor with two levels \code{no} and
#' \code{yes}.} \item{prevStroke}{History of previous strokes before admission,
#' a factor with two levels \code{no} and \code{yes}.}
#' \item{othDisease}{History of other disabling diseases (e.g. severe
#' dementia), a factor with two levels \code{no} and \code{yes}.}
#' \item{alcohol}{Daily alcohol consumption, a factor with two levels \code{no}
#' and \code{yes}.} \item{diabetes}{Diabetes mellitus status indicating if the
#' glucose level was higher than 11 mmol/L, a factor with two levels \code{no}
#' and \code{yes}.} \item{smoke}{Daily smoking status, a factor with two levels
#' \code{no} and \code{yes}.} \item{atrialFib}{Atrial fibrillation, a factor
#' with two levels \code{no} and \code{yes}.} \item{hemor}{Hemorrhage (stroke
#' subtype), a factor with two levels \code{no} (infarction) and \code{yes}
#' (hemorrhage).} \item{strokeScore}{Scandinavian stroke score at admission to
#' the hospital. Ranges from 0 (worst) to 58 (best).}
#' \item{cholest}{Cholesterol level} \item{time}{Survival time (in days).}
#' \item{status}{Status (\code{0}: censored, \code{1}: event).} }
#' @references Joergensen HS, Nakayama H, Reith J, Raaschou HO, and Olsen TS.
#' Acute stroke with atrial fibrillation. The Copenhagen Stroke Study. Stroke,
#' 27(10):1765-9, 1996.
#' 
#' Mogensen UB, Ishwaran H, and Gerds TA. Evaluating random forests for
#' survival analysis using prediction error curves. Technical Report 8,
#' University of Copenhagen, Department of Biostatistics, 2010.
#' @useDynLib pec
#' @importFrom foreach %dopar%
#' @importFrom survival Surv
#' @importFrom prodlim Hist
#' @importFrom grDevices col2rgb gray
#' @importFrom graphics abline axis box legend lines mtext par plot points segments text title
#' @importFrom stats model.frame model.response as.formula coef family formula median model.matrix na.fail na.omit pnorm predict quantile rbinom rexp runif sd smooth terms time update update.formula var wilcox.test
#' @importFrom utils capture.output head select.list
#' @keywords datasets
NULL


#' German Breast Cancer Study Group 2
#' 
#' A data frame containing the observations from the GBSG2 study.
#' 
#' 
#' @name GBSG2
#' @docType data
#' @format This data frame contains the observations of 686 women: \describe{
#' \item{horTh}{hormonal therapy, a factor at two levels \code{no} and
#' \code{yes}.} \item{age}{of the patients in years.}
#' \item{menostat}{menopausal status, a factor at two levels \code{pre}
#' (premenopausal) and \code{post} (postmenopausal).} \item{tsize}{tumor size
#' (in mm).} \item{tgrade}{tumor grade, a ordered factor at levels \code{I < II
#' < III}.} \item{pnodes}{number of positive nodes.}
#' \item{progrec}{progesterone receptor (in fmol).} \item{estrec}{estrogen
#' receptor (in fmol).} \item{time}{recurrence free survival time (in days).}
#' \item{cens}{censoring indicator (0- censored, 1- event).} }
#' @references M. Schumacher, G. Basert, H. Bojar, K. Huebner, M. Olschewski,
#' W. Sauerbrei, C. Schmoor, C. Beyerle, R.L.A. Neumann and H.F. Rauschecker
#' for the German Breast Cancer Study Group (1994), Randomized \eqn{2\times2}
#' trial evaluating hormonal treatment and the duration of chemotherapy in
#' node-positive breast cancer patients.  \emph{Journal of Clinical Oncology},
#' \bold{12}, 2086--2093.
#' @keywords datasets
NULL


#' Pbc3 data
#' 
#' PBC3 was a multi-centre randomized clinical trial conducted in six European
#' hospitals. Between 1 Jan. 1983 and 1 Jan. 1987, 349 patients with the liver
#' disease primary biliary cirrhosis (PBC) were randomized to either treatment
#' with Cyclosporin A (CyA, 176 patients) or placebo (173 patients). The
#' purpose of the trial was to study the effect of treatment on the survival
#' time. However, during the course of the trial an increased use of liver
#' transplantation for patients with this disease made the investigators
#' redefine the main response variable to be time to ``failure of medical
#' treatment'' defined as either death or liver transplantation. Patients were
#' then followed from randomization until treatment failure, drop-out or 1 Jan,
#' 1989; 61 patients died (CyA: 30, placebo: 31), another 29 were transplanted
#' (CyA: 14, placebo: 15) and 4 patients were lost to follow-up before 1 Jan.
#' 1989. At entry a number of clinical, biochemical and histological variables,
#' including serum bilirubin, serum albumin, sex, age were recorded.
#' 
#' 
#' @name Pbc3
#' @docType data
#' @format A data frame with 349 observations on the following 15 variables.
#' \describe{ \item{ptno}{patient identification}
#' \item{unit}{hospital (1: Hvidovre, 2: London, 3: Copenhagen, 4:
#' Barcelona, 5: Munich, 6: Lyon)} \item{tment}{treatment (0: placebo,
#' 1: CyA)} \item{sex}{(1: males, 0: females)} \item{age}{age
#' in years} \item{stage}{histological stage (1, 2, 3, 4)}
#' \item{gibleed}{previous gastrointestinal bleeding (1: yes, 0: no)}
#' \item{crea}{creatinine (micromoles/L)} \item{alb}{albumin
#' (g/L)} \item{bili}{bilirubin (micromoles/L)}
#' \item{alkph}{alkaline phosphatase (IU/L)}
#' \item{asptr}{aspartate transaminase (IU/L)}
#' \item{weight}{body weight (kg)} \item{days}{observation time
#' (days)} \item{status}{status at observation time (0: censored, 1:
#' liver transplantation, 2 : dead)} }
#' @references Andersen and Skovgaard. Regression with linear predictors.
#' Springer, 2010.
#' @source
#' 
#' \url{http://192.38.117.59/~linearpredictors/?page=datasets&dataset=Pbc3}
#' @keywords datasets
#' @examples
#' 
#' data(Pbc3)
#' 
NULL



