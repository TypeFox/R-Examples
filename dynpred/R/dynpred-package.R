

#' Data from the European Society for Blood and Marrow Transplantation (EBMT)
#' 
#' @format A data frame of 2279 patients transplanted at the EBMT between 1985 and
#' 1998. These data were used in Fiocco, Putter & van Houwelingen (2008) and
#' van Houwelingen & Putter (2008). The included variables are \describe{
#' \item{id}{Patient identification number} \item{rec}{Time in days from
#' transplantation to recovery or last follow-up} \item{rec.s}{Recovery status;
#' 1 = recovery, 0 = censored} \item{ae}{Time in days from transplantation to
#' adverse event (AE) or last follow-up} \item{ae.s}{Adverse event status; 1 =
#' adverse event, 0 = censored} \item{recae}{Time in days from transplantation
#' to both recovery and AE or last follow-up} \item{plag.s}{Recovery and AE
#' status; 1 = both recovery and AE, 0 = no recovery or no AE or censored}
#' \item{rel}{Time in days from transplantation to relapse or last follow-up}
#' \item{rel.s}{Relapse status; 1 = relapse, 0 = censored} \item{srv}{Time in
#' days from transplantation to death or last follow-up} \item{srv.s}{Relapse
#' status; 1 = dead, 0 = censored} \item{year}{Year of transplantation; factor
#' with levels "1985-1989", "1990-1994", "1995-1998"} \item{agecl}{Patient age
#' at transplant; factor with levels "<=20", "20-40", ">40"}
#' \item{proph}{Prophylaxis; factor with levels "no", "yes"}
#' \item{match}{Donor-recipient gender match; factor with levels "no gender
#' mismatch", "gender mismatch"} }
#' 
#' 
#' @name EBMT data
#' @aliases ALL
#' @docType data
#' @references Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
#' proportional hazards regression and simulation-based prediction for
#' multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' 
#' van Houwelingen HC, Putter H (2008). Dynamic predicting by landmarking as an
#' alternative for multi-state modeling: an application to acute lymphoid
#' leukemia data. \emph{Lifetime Data Anal} \bold{14}, 447--463.
#' @source We gratefully acknowledge the European Society for Blood and Marrow
#' Transplantation (EBMT) for making available these data. Disclaimer: these
#' data were simplified for the purpose of illustration of the analysis of
#' competing risks and multi-state models and do not reflect any real life
#' situation. No clinical conclusions should be drawn from these data.
#' @keywords datasets
NULL





#' The companion package of the book "Dynamic Prediction in Survival Analysis"
#' 
#' The companion package of the book "Dynamic Prediction in Survival Analysis".
#' 
#' \tabular{ll}{ Package: \tab dynpred\cr Type: \tab Package\cr Version: \tab
#' 0.1.2\cr Date: \tab 2014-11-10\cr License: \tab GPL (>= 2)\cr } An overview
#' of how to use the package, including the most important functions.
#' 
#' @name dynpred-package
#' @aliases dynpred-package dynpred
#' @docType package
#' @author Hein Putter Maintainer: Hein Putter <H.Putter@@lumc.nl>
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#' Clinical Survival Analysis. Chapman & Hall.
#' @keywords package
NULL





#' Clinical and follow-up data of breast cancer patients as collected in the
#' Dutch Cancer Institute (NKI) in Amsterdam
#' 
#' A data frame of 295 patients with breast cancer. The included variables are
#' \describe{ \item{patnr}{Patient identification number} \item{d}{Survival
#' status; 1 = death; 0 = censored} \item{tyears}{Time in years until death or
#' last follow-up} \item{diameter}{Diameter of the primary tumor}
#' \item{posnod}{Number of positive lymph nodes} \item{age}{Age of the patient}
#' \item{mlratio}{Estrogen level?} \item{chemotherapy}{Chemotherapy used
#' (yes/no)} \item{hormonaltherapy}{Hormonal therapy used (yes/no)}
#' \item{typesurgery}{Type of surgery (excision or mastectomy)}
#' \item{histolgrade}{Histological grade (Intermediate, poorly, or well
#' differentiated)} \item{vasc.invasion}{Vascular invasion (-, +, or +/-)}
#' \item{crossval.clin.class}{??} \item{PICV}{Estrogen level?} }
#' 
#' 
#' @name NKI breast cancer clinical data
#' @aliases nki
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references van't Veer LJ, Dai HY, van de Vijver MJ, He YDD, Hart AAM, Mao
#' M, Peterse HL, van der Kooy K, Marton MJ, Witteveen AT, Schreiber GJ,
#' Kerkhoven RM, Roberts C, Linsley PS, Bernards R \& Friend SH (2002). Gene
#' expression profiling predicts clinical outcome of breast cancer.
#' \emph{Nature} \bold{415}, 530--536.
#' 
#' van de Vijver MJ, He YD, van `t Veer LJ, Dai H, Hart AAM, Voskuil DW,
#' Schreiber GJ, Peterse JL, Roberts C, Marton MJ, Parrish M, Atsma D,
#' Witteveen A, Glas A, Delahaye L, van der Velde T, Bartelink H, Rodenhuis S,
#' Rutgers ET, Friend SH \& Bernards R (2002). A gene-expression signature as a
#' predictor of survival in breast cancer. \emph{New England Journal of
#' Medicine} \bold{347}, 1999--2009.
#' 
#' van Houwelingen HC, Bruinsma T, Hart AAM, van't Veer LJ \& Wessels LFA
#' (2006).  Cross-validated Cox regression on microarray gene expression data.
#' \emph{Statistics in Medicine} \bold{25}, 3201--3216.
#' @keywords datasets
NULL





#' Data originate from two clinical trials on the use of different combination
#' chemotherapies, carried out in The Netherlands around 1980
#' 
#' A data frame of 358 patients with ovarian cancer. The included variables are
#' \describe{ \item{tyears}{Time in years until death or last follow-up}
#' \item{d}{Survival status; 1 = death; 0 = censored} \item{Karn}{Karnofsky
#' score} \item{Broders}{Broders score: factor with levels "unknown", "1", "2",
#' "3", "4"} \item{FIGO}{FIGO stage; factor with levels "III", "IV"}
#' \item{Ascites}{Presence of ascires; factor with levels "unknown", "absent",
#' "present"} \item{Diam}{Diameter of the tumor; factor with levels "micr.",
#' "<1cm", "1-2cm", "2-5cm", ">5cm"} }
#' 
#' 
#' @name Ovarian cancer data
#' @aliases ova
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Neijt, J. P., ten Bokkel Huinink, W. W., van der Burg, M. E.,
#' van Oosterom, A. T., Vriesendorp, R., Kooyman, C. D., van Lindert, A. C.,
#' Hamerlynck, J. V., van Lent, M. & van Houwelingen, J. C. (1984), `Randomised
#' trial comparing two combination chemotherapy regimens (Hexa-CAF vs CHAP- 5)
#' in advanced ovarian carcinoma', Lancet 2, 594--600.
#' 
#' Neijt, J. P., ten Bokkel Huinink, W. W., van der Burg, M. E., van Oosterom,
#' A. T., Willemse, P. H., Heintz, A. P., van Lent, M., Trimbos, J. B., Bouma,
#' J. & Vermorken, J. B. (1987), `Randomized trial comparing two combination
#' chemotherapy regimens (CHAP-5 vs CP) in advanced ovarian carcinoma', Journal
#' of Clinical Oncology 5, 1157--1168.
#' 
#' van Houwelingen, J. C., ten Bokkel Huinink, W. W., van der Burg, M. E., van
#' Oosterom, A. T. & Neijt, J. P. (1989), `Predictability of the survival of
#' patients with advanced ovarian cancer.', Journal of Clinical Oncology 7,
#' 769--773.
#' @keywords datasets
NULL





#' Data from the Benelux CML study
#' 
#' A data frame of 210 patients with Chronic Myeloid Leukemia from the Benelux
#' CML study (Kluin-Nelemans et al. 1998). Data have been used in two
#' methodological papers, de Bruijne et al. (2001) and van Houwelingen (2007),
#' and in the book van Houwelingen \& Putter (2011), especially Chapter 8. More
#' background is given in Appendix A.2 of van Houwelingen \& Putter (2011).
#' Interest is in the time-dependent covariate White Blood Cell count (WBC).
#' Data set wbc1 contains the follow-up data and time-fixed covariates, while
#' \code{\link{wbc2}} contains the WBC measurements. The included variables in
#' wbc1 are \describe{ \item{patnr}{Patient identification number}
#' \item{tyears}{Time in years from randomization to death or last follow-up}
#' \item{d}{Survival status; 1 = dead, 0 = censored} \item{sokal}{Clinical
#' index based on spleen size, percentage of circulating blasts, platelet and
#' age at diagnosis} \item{age}{Age at diagnosis} }
#' 
#' 
#' @name WBC follow-up data
#' @aliases wbc1
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Kluin-Nelemans JC, Delannoy A, Louwagie A, le Cessie S, Hermans
#' J, van der Burgh JF, Hagemeijer AM, van den Berghe H \& Benelux CML Study
#' Group (1998). Randomized study on hydroxyurea alone versus hydroxyurea
#' combined with low-dose interferon-alpha 2b for chronic myeloid leukemia.
#' \emph{Blood} \bold{91}, 2713--2721.
#' 
#' de Bruijne MHJ, le Cessie S, Kluin-Nelemans HC \& van Houwelingen HC (2001).
#' On the use of Cox regression in the presence of an irregularly observed
#' time-dependent covariate. \emph{Statistics in Medicine} \bold{20},
#' 3817--3829.
#' 
#' van Houwelingen HC (2007). Dynamic prediction by landmarking in event
#' history analysis. \emph{Scandinavian Journal of Statistics} \bold{34},
#' 70--85.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Predicting in Clinical Survival
#' Analysis. Chapman \& Hall.
#' @keywords datasets
NULL





#' Data from the Benelux CML study
#' 
#' A data frame of 210 patients with Chronic Myeloid Leukemia from the Benelux
#' CML study (Kluin-Nelemans et al. 1998). Data have been used in two
#' methodological papers, de Bruijne et al. (2001) and van Houwelingen (2007),
#' and in the book van Houwelingen \& Putter (2011), especially Chapter 8. More
#' background is given in Appendix A.2 of van Houwelingen \& Putter (2011).
#' Interest is in the time-dependent covariate White Blood Cell count (WBC).
#' Data set \code{\link{wbc1}} contains the follow-up data and time-fixed
#' covariates, while wbc2 contains the WBC measurements. The included variables
#' in wbc2 are \describe{ \item{patnr}{Patient identification number}
#' \item{tyears}{Time of WBC measurement in years from randomization}
#' \item{lwbc}{Log-transformed and standardized WBC measurement, more
#' precisely, defined as lwbc=log10(wbc)-0.95} }
#' 
#' 
#' @name WBC measurements data
#' @aliases wbc2
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Kluin-Nelemans JC, Delannoy A, Louwagie A, le Cessie S, Hermans
#' J, van der Burgh JF, Hagemeijer AM, van den Berghe H \& Benelux CML Study
#' Group (1998). Randomized study on hydroxyurea alone versus hydroxyurea
#' combined with low-dose interferon-alpha 2b for chronic myeloid leukemia.
#' \emph{Blood} \bold{91}, 2713--2721.
#' 
#' de Bruijne MHJ, le Cessie S, Kluin-Nelemans HC \& van Houwelingen HC (2001).
#' On the use of Cox regression in the presence of an irregularly observed
#' time-dependent covariate. \emph{Statistics in Medicine} \bold{20},
#' 3817--3829.
#' 
#' van Houwelingen HC (2007). Dynamic prediction by landmarking in event
#' history analysis. \emph{Scandinavian Journal of Statistics} \bold{34},
#' 70--85.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Predicting in Clinical Survival
#' Analysis. Chapman \& Hall.
#' @keywords datasets
NULL



