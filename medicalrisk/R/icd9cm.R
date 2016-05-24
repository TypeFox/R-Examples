#' Create Deyo map of ICD-9-CM to Charlson comorbidities
#' 
#' Function that generates a data frame linking ICD-9-CM codes to the Charlson
#' comorbidity categories using the Deyo mapping.
#' 
#' NOTE: The input vector of ICD-9-CM codes must be unique, because the output dataframe
#' uses the ICD-9-CM code as row.name.
#' 
#' Uses regular expressions created from the paper by Deyo in 1992.
#' 
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and one logical column for each 
#' comorbidity in \code{\link{charlson_list}}
#' @references 1. Deyo RA, Cherkin DC, Ciol MA: Adapting a clinical comorbidity index for
#' use with ICD-9-CM administrative databases. Journal of clinical epidemiology
#' 1992; 45:613-9
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/1607900}
#
#' @seealso \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_charlson_romano}},
#'    \code{\link{icd9cm_elixhauser_quan}}, \code{\link{icd9cm_elixhauser_ahrq37}},
#'    \code{\link{charlson_weights}},
#' @examples
#' # Identify Charlson categories in ICD-9-CM listing
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' cases_with_cm <- merge(cases, icd9cm_charlson_deyo(levels(cases$icd9cm)), 
#'   by.x="icd9cm", by.y="row.names", all.x=TRUE)
#' 
#' # generate crude comorbidity summary for each patient
#' library(plyr)
#' ddply(cases_with_cm, .(id), 
#'   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) })
#' @export
icd9cm_charlson_deyo <- function(icd9) {
  data.frame(
    row.names = icd9,
    mi = grepl('^D(410|412)', icd9),
    chf = grepl('^D428[0-9]', icd9),
    perivasc = grepl('^(D(4439|441[0-9]|7854|V434)|P3848)', icd9),
    cvd = grepl('^D43[0-8]', icd9),
    dementia = grepl('^D290[0-9]', icd9),
    chrnlung = grepl('^D(49[0-6]|50[0-5]|5064)', icd9),
    rheum = grepl('^D(710[0,1,4]|714[0-2]|71481|725)', icd9),
    ulcer = grepl('^D53[1-4]', icd9),
    liver = grepl('^D571[2456]',  icd9),
    dm = grepl('^D250[0-3,7]', icd9),
    dmcx = grepl('^D250[4-6]', icd9),
    para = grepl('^D(3441|342)', icd9),
    renal = grepl('^D(58[2568]|583[0-7])', icd9),
    tumor = grepl('^D(1[4-6]|17[0-24-9]|18|19[0-5]|20[0-8])', icd9),
    modliver = grepl('^D(456[0-2]|572[2-8])', icd9),
    mets = grepl('^D19[6-9]', icd9),
    aids = grepl('^D04[2-4]', icd9))
}

#' Create Romano map of ICD-9-CM to Charlson comorbidities
#' 
#' Function that creates a dataframe which links ICD-9-CM codes to the Charlson 
#' comorbidity categories using the Romano mapping.
#' 
#' NOTE: The input vector of ICD-9-CM codes must be unique, because the output dataframe
#' uses the ICD-9-CM code as row.name.
#' 
#' Uses regular expressions created from the paper by Romano in 1993.
#' 
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and one logical column for each 
#' comorbidity in \code{\link{charlson_list}}
#' @references 1. Romano PS, Roos LL, Jollis JG: Adapting a clinical comorbidity index for
#' use with ICD-9-CM administrative data: differing perspectives. Journal of
#' clinical epidemiology 1993; 46:1075-9; discussion 1081-90
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/8410092}
#
#' @seealso \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_charlson_deyo}},
#'    \code{\link{icd9cm_elixhauser_quan}}, \code{\link{icd9cm_elixhauser_ahrq37}},
#'    \code{\link{charlson_weights}},
#' @examples
#' # Identify Charlson categories in ICD-9-CM listing
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' cases_with_cm <- merge(cases, icd9cm_charlson_romano(levels(cases$icd9cm)), 
#'   by.x="icd9cm", by.y="row.names", all.x=TRUE)
#' 
#' # generate crude comorbidity summary for each patient
#' library(plyr)
#' ddply(cases_with_cm, .(id), 
#'   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) })
#' @export
icd9cm_charlson_romano <- function(icd9) {
  data.frame(
    row.names = icd9,
    mi = grepl('^D(410|412)', icd9),
    chf = grepl('^D(402[01]1|40291|425|428|4293)', icd9),
    perivasc = grepl('^(D(44[012]|443[1-9]|4471|7854)|P(38(1[3468]|3[3468]|4[3468])|392[2-69]))', icd9),
    cvd = grepl('^(D(36234|43[0-6]|437[019]|438|7814|7843|9970)|P38[14]2)', icd9),
    dementia = grepl('^D(290[0-9]|331[0-2])', icd9),
    chrnlung = grepl('^D(4150|416[89]|49[1-46])', icd9),
    rheum = F,
    ulcer = grepl('^D53[1-4]', icd9),
    liver = grepl('^D571[25689]',  icd9),
    dm = grepl('^D250[0-3]', icd9),
    dmcx = grepl('^D250[4-9]', icd9),
    para = grepl('^D(34[2,4])', icd9),
    renal = grepl('^(D(58[56]|V420|V451|V56)|P(39(27|42|9[3-5])|5498))', icd9),
    tumor = grepl('^(D(1[4-6]|17[0-14-9]|18|19[0-5]|20[0-8]|273[03]|V1046)|P(605|6241))', icd9),
    modliver = grepl('^(D(456[0-2]|572[2-4])|P(391|4291))', icd9),
    mets = grepl('^D19[6-9]', icd9),
    aids = grepl('^D04[2-4]', icd9))
}

#' Create Quan map of ICD-9-CM to Charlson comorbidities
#' 
#' Function that creates a dataframe that links ICD-9-CM codes to the Charlson comorbidity 
#' categories using Quan's method.
#' 
#' NOTE: The input vector of ICD-9-CM codes must be unique, because the output dataframe
#' uses the ICD-9-CM code as row.name.
#' 
#' Uses regular expressions created from the paper by Quan in 2005.
#' 
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and one logical column for each 
#' comorbidity in \code{\link{charlson_list}}
#' @references 1. Quan H, Sundararajan V, Halfon P, Fong A, Burnand B, Luthi 
#'   J-C, Saunders LD, Beck CA, Feasby TE, Ghali WA: Coding algorithms for 
#'   defining comorbidities in ICD-9-CM and ICD-10 administrative data. Medical 
#'   care 2005; 43:1130-9 
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/16224307}
#
#' @seealso \code{\link{icd9cm_charlson_deyo}}, \code{\link{icd9cm_charlson_romano}},
#'    \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_elixhauser_quan}}
#' @examples
#' # Identify Charlson categories in ICD-9-CM listing
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' cases_with_cm <- merge(cases, icd9cm_charlson_quan(levels(cases$icd9cm)), 
#'   by.x="icd9cm", by.y="row.names", all.x=TRUE)
#' 
#' # generate crude comorbidity summary for each patient
#' library(plyr)
#' ddply(cases_with_cm, .(id), 
#'   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) })
#' @export
icd9cm_charlson_quan <- function(icd9) { 
  data.frame(
    row.names = icd9,
    mi = grepl('^D(410|412)', icd9),
    chf = grepl('^D(39891|402[019]1|404[019][13]|425[4-9]|428)', icd9),
    perivasc = grepl('^D(0930|4373|44[01]|443[1-9]|4471|557[19]|V434)', icd9),
    cvd = grepl('^D(36234|43[0-8])', icd9),
    dementia = grepl('^D(290|2941|3312)', icd9),
    chrnlung = grepl('^D(416[89]|49|50[0-5]|5064|508[18])', icd9),
    rheum = grepl('^D(4465|710[0-4]|714[0-2]|7148|725)', icd9),
    ulcer = grepl('^D53[1-4]', icd9),
    liver = grepl('^D(070[23][23]|070[45]4|070[69]|57[01]|573[3489]|V427)',  icd9),
    dm = grepl('^D250[0-389]', icd9),
    dmcx = grepl('^D250[4-7]', icd9),
    para = grepl('^D(3341|34[23]|344[0-69])', icd9),
    renal = grepl('^D(403[019]1|404[019][23]|582|583[0-7]|58[56]|5880|V420|V451|V56)', icd9),
    tumor = grepl('^D(1[4-6]|17[0-24-9]|18|19[0-5]|20[0-8]|2386)', icd9),
    modliver = grepl('^D(456[0-2]|572[2-8])', icd9),
    mets = grepl('^D19[6-9]', icd9),
    aids = grepl('^D04[2-4]', icd9))
}

#' Create Quan map of ICD-9-CM to Elixhauser comorbidities
#' 
#' Function to make a dataframe that links ICD-9-CM codes to the Elixhauser comorbidity 
#' categories using the Quan mapping.
#' 
#' Uses regular expressions created from the Quan paper from 2005.
#' 
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' Some ICD-9-CM codes will correspond to more than one category.  For example, 
#' \code{404.03} (Hypertensive heart and chronic kidney disease ... stage V) is 
#' in both \code{chf} and \code{renlfail} categories.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and one logical column for each 
#' comorbidity in \code{\link{elixhauser_list}}
#' @references 1. Quan H, Sundararajan V, Halfon P, Fong A, Burnand B, Luthi 
#'   J-C, Saunders LD, Beck CA, Feasby TE, Ghali WA: Coding algorithms for 
#'   defining comorbidities in ICD-9-CM and ICD-10 administrative data. Medical 
#'   care 2005; 43:1130-9 
#'   \url{http://www.ncbi.nlm.nih.gov/pubmed/16224307}
#
#' @seealso \code{\link{icd9cm_charlson_deyo}}, \code{\link{icd9cm_charlson_romano}},
#'    \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_elixhauser_ahrq37}}
#' @examples
#' # Identify Elixhauser categories
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' cases_with_cm <- merge(cases, icd9cm_elixhauser_quan(levels(cases$icd9cm)), 
#'   by.x="icd9cm", by.y="row.names", all.x=TRUE)
#' 
#' # generate crude comorbidity summary for each patient
#' library(plyr)
#' ddply(cases_with_cm, .(id), 
#'   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) })
#' @export
icd9cm_elixhauser_quan <- function(icd9) {
  data.frame(
    row.names = icd9,
    chf = grepl('^D(39891|402[019]1|404[019][13]|425[4-9]|428)', icd9),
    arrhythmia = grepl('^D(426([079]|13)|426(10|12)|427[0-46-9]|7850|9960[14]|V450|V533)', icd9),
    valve = grepl('^D(0932|39[4567]|424|746[3-6]|V422|V433)', icd9),
    pulmcirc = grepl('^D(415[01]|416|417[089])', icd9),
    perivasc = grepl('^D(0930|4373|44[01]|443[1-9]|4471|557[19]|V434)', icd9),
    htn = grepl('^D401', icd9),
    htncx = grepl('^D40[2-5]', icd9),
    para = grepl('^D(3341|34[23]|344[0-69])', icd9),
    neuro = grepl('^D(3319|332[01]|333[45]|33392|33[45]|3362|34[015]|348[13]|7803|7843)', icd9),
    chrnlung = grepl('^D(416[89]|49|50[0-5]|5064|508[18])', icd9),
    dm = grepl('^D250[0-3]', icd9),
    dmcx = grepl('^D250[4-9]', icd9),
    hypothy = grepl('^D(2409|24[34]|246[18])', icd9),
    renlfail = grepl('^D(403[019]1|404[019][23]|58[56]|5880|V420|V451|V56)', icd9),
    liver = grepl('^D(070[23][23]|070[45]4|070[69]|456[0-2]|57[01]|572[2-8]|573[3489]|V427)',  icd9),
    ulcer = grepl('^D53[1-4][79]', icd9),
    aids = grepl('^D04[2-4]', icd9),
    lymph = grepl('^D(20[0-2]|2030|2386)', icd9),
    mets = grepl('^D19[6-9]', icd9),
    tumor = grepl('^D(1[4-6]|17[0-24-9]|18|19[0-5])', icd9),
    rheum = grepl('^D(446|7010|710[0-489]|7112|714|7193|72[05]|7285|72889|72930)', icd9),
    coag = grepl('^D(286|287[13-5])', icd9),
    obese = grepl('^D2780', icd9),
    wghtloss = grepl('^D(26[0-3]|7832|7994)', icd9),
    lytes = grepl('^D(2536|276)', icd9),
    bldloss = grepl('^D2800', icd9),
    anemdef = grepl('^D(280[1-9]|281)', icd9),
    alcohol = grepl('^D(2652|291[1-35-9]|303[09]|3050|3575|4255|5353|571[0-3]|980|V113)', icd9),
    drug = grepl('^D(292|304|305[2-9]|V6542)', icd9),
    psych = grepl('^D(2938|295|296[0145]4|29[78])', icd9),
    depress = grepl('^D(296[235]|3004|309|311)', icd9))  
}

#' Create AHRQ v3.7 map of ICD-9-CM to Elixhauser comorbidities
#' 
#' Function makes a dataframe that links ICD-9-CM codes to the Elixhauser comorbidity 
#' categories using the AHRQ v3.7 mapping.
#' 
#' Uses regular expressions based on the file "comformat2012-2013.txt" from AHRQ.
#' 
#' The Agency for Healthcare Research and Quality (AHRQ) has developed 
#' Comorbidity Software as part of the Healthcare Cost and Utilization Project 
#' (HCUP).  The software was developed to report on the comorbidity measures 
#' reported by Elixhauser (1998).
#' 
#' The AHRQ software has two parts, one that classifies ICD-9-CM codes by 
#' comorbidity, and another that performs heuristics to eliminate duplicate 
#' comorbidities and ignore comorbidities which are the primary reason for the 
#' hospital visit, as per the DRG.
#' 
#' This table is a translation of the first part of the software, the 
#' classifier, as implemented in the SAS file \code{Comformat2012-2013.txt}.
#' 
#' ICD-9-CM codes must have periods removed.  Diagnostic codes are prefixed with 
#' 'D' while procedure codes are prefixed with 'P'. So, diagnostic code
#' \code{404.03} should be \code{"D40403"}.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and one logical column for each 
#' comorbidity in \code{\link{elixhauser_list}}
#' @references 1. \url{http://www.hcup-us.ahrq.gov/toolssoftware/comorbidity/comorbidity.jsp}
#' 
#' @seealso \code{\link{icd9cm_charlson_deyo}}, \code{\link{icd9cm_charlson_romano}},
#'    \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_elixhauser_quan}}
#' @examples
#' # Identify Elixhauser categories
#' cases <- data.frame(id=c(1,1,1,2,2,2),
#'   icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403"))
#' cases_with_cm <- merge(cases, icd9cm_elixhauser_ahrq37(levels(cases$icd9cm)), 
#'   by.x="icd9cm", by.y="row.names", all.x=TRUE)
#' 
#' # generate crude comorbidity summary for each patient
#' library(plyr)
#' ddply(cases_with_cm, .(id), 
#'   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) })
#' @export
icd9cm_elixhauser_ahrq37 <- function(icd9) {
  # the below is taken from comformat2012-2013.txt
  # htncx: htncx, htnpreg, htnwochf, hrenworf, hhrwohrf
  # htncx & chf: htnwchf, hhrwchf
  # htncx & renlfail: hrenwrf, hhrwrf
  # htncx & chf & renlfail: hhrwhrf
  htnpreg <- grepl('^D6422', icd9)
  htnwochf <- grepl('^D(402[019]0|405[019]9)', icd9)
  htnwchf <- grepl('^D402[019]1', icd9)
  hrenworf <- grepl('^D(403[019]0|405[019]1|6421)', icd9)
  hrenwrf <- grepl('^D403[019]1', icd9)
  hhrwohrf <- grepl('^D404[019]0', icd9)
  hhrwchf <- grepl('^D404[019]1', icd9)
  hhrwrf <- grepl('^D404[019]2', icd9)
  hhrwhrf <- grepl('^D404[019]3', icd9)
  ohtnpreg <- grepl('^D642[79]', icd9)
  
  data.frame(
    row.names = icd9,
    chf = grepl('^D(39891|428)', icd9) | htnwchf | hhrwchf | hhrwhrf,
    arrhythmia = F,
    valve = grepl('^D(0932|39[456]|397[019]|424|746[3-6]|V422|V433)', icd9),
    pulmcirc = grepl('^D(4151|416|4179)', icd9),
    perivasc = grepl('^D(44[012]|443[1-9]|4442|4471|449|557[19]|V434)', icd9),
    htn = grepl('^D(401[19]|6420)', icd9),
    htncx = grepl('^D(4010|4372)', icd9) | 
      htnpreg | htnwochf | hrenworf | hhrwohrf |
      htnwchf | hhrwchf | hrenwrf | hhrwrf | hhrwhrf,
    para = grepl('^D(34[234]|438[2-5]|78072)', icd9),    
    neuro = grepl('^D(33[01]|3320|333[457]|33385|33394|33[45]|3380|34[0157]|6494|7687|7803|78097|7843)', icd9),
    chrnlung = grepl('^D(49|50[0-5]|5064)', icd9),
    dm = grepl('^D(250[0-3]|6480|249[0-3])', icd9),
    dmcx = grepl('^D(250[4-9]|7751|249[4-9])', icd9),
    hypothy = grepl('^D(243|244[01289])', icd9),
    renlfail = grepl('^D(585[3-9]|586|V420|V451|V56)', icd9) | hrenwrf | hhrwrf,
    liver = grepl('^D(070[23][23]|070[45]4|456[0-2]|571[02-9]|572[38]|5735|V427)', icd9),
    ulcer = grepl('^D(53[1234][456]1|53[1234]7[01]|53[1234]91)', icd9),
    aids = grepl('^D04[2-4]', icd9),
    lymph = grepl('^D(20[01]|202[0-35-9]|203|2386|2733)', icd9),
    mets = grepl('^D(19[6-8]|199[01]|2097|78951)', icd9),
    tumor = grepl('^D(1[4-6]|17[012459]|18|19[0-5]|209[0-3]|2580)', icd9),
    rheum = grepl('^D(7010|71[04]|720|725)', icd9),
    coag = grepl('^D(286|287[13-5]|6493|28984)', icd9),
    obese = grepl('^D(2780|6491|V85[34]|V8554|79391)', icd9),
    wghtloss = grepl('^D(26[0-3]|7832)', icd9),
    lytes = grepl('^D276', icd9),
    bldloss = grepl('^D(2800|6482)', icd9),
    anemdef = grepl('^D(280[1-9]|281|285[29])', icd9),
    alcohol = grepl('^D(291[0-35-9]|303[09]|3050)', icd9),
    drug = grepl('^D(292|304|305[2-9]|6483)', icd9),
    psych = grepl('^D(29[5678]|2991)', icd9),
    depress = grepl('^D(3004|30112|309[01]|311)', icd9))  
}

#' Create Map of ICD-9-CM to Revised Cardiac Risk Index classes
#' 
#' Function to generate data frame that links ICD-9-CM codes to the RCRI comorbidity categories.
#' 
#' Lee et al in 1999 published a "Revised Cardiac Risk Index" based on the work 
#' on Goldman in 1997. The RCRI is used to determine the major cardiac 
#' complication risk for a patient about to undergo major noncardiac surgery. 
#' The six predictors that make up the RCRI are: 1. high-risk surgery 2. history
#' of ischemic heart disease 3. history of congestive heart failure 4. history 
#' of cerebrovascular disease 5. preoperative treatment with insulin 6. 
#' preoperative serum creatinine with Cr > 2 mg/dL.
#' 
#' In 2005 Boersma et al demonstrated that the Lee indexed can be adapted to use
#' administrative data to predict cardiovascular mortality.  They used the 
#' following for each point above: 1. retroperitoneal, intrathoracic, or 
#' suprainguinal vascular procedure; 2. Ischemia: ICD-9 codes 410.*, 411.*, 
#' 412.*, 413.*, 414.*; 3. CVA: ICD-9 428.*; 4. CHF: ICD-9 943.0; 5. DM: ICD-9 
#' 425.0; 6. Renal: ICD-9 958.0.
#' 
#' This function merges the ICD-9 guidelines 
#' used by Boersma with some of the other ICD-9 classifiers in this package. 
#' This data set uses the following for each aspect of the RCRI:
#' 1. procedure is left to you
#' 2. 'ischemia' as defined in Boersma
#' 3. 'cvd' as defined by Quan in \code{\link{icd9cm_charlson_quan}}
#' 4. 'chf' as defined by AHRQ in \code{\link{icd9cm_elixhauser_ahrq37}}
#' 5. 'dm' as  defined by AHRQ (both 'dm' and 'dmcx')
#' 6. renlfail' as defined by AHRQ.
#' 
#' @param icd9 a unique character vector of ICD-9-CM codes
#' @return A data frame, with ICD9 codes as row names and logical columns for
#' \code{chf}, \code{cvd}, \code{dm}, \code{ischemia}, and \code{renlfail}.
#' @references 1. Lee TH, Marcantonio ER, Mangione CM, Thomas EJ, Polanczyk CA, 
#'   Cook EF, Sugarbaker DJ, Donaldson MC, Poss R, Ho KK, Ludwig LE, Pedan A, 
#'   Goldman L: Derivation and prospective validation of a simple index for 
#'   prediction of cardiac risk of major noncardiac surgery. Circulation 1999; 
#'   100:1043-9 \url{http://www.ncbi.nlm.nih.gov/pubmed/10477528}
#'   
#'   2. Boersma E, Kertai MD, Schouten O, Bax JJ, Noordzij P, Steyerberg EW, 
#'   Schinkel AFL, Santen M van, Simoons ML, Thomson IR, Klein J, Urk H van, 
#'   Poldermans D: Perioperative cardiovascular mortality in noncardiac surgery:
#'   validation of the Lee cardiac risk index. The American journal of medicine 
#'   2005; 118:1134-41 \url{http://www.ncbi.nlm.nih.gov/pubmed/16194645}
#'   
#' @seealso \code{\link{icd9cm_charlson_quan}}, \code{\link{icd9cm_elixhauser_quan}}, 
#'    \code{\link{icd9cm_elixhauser_ahrq37}}
#' @export
icd9cm_rcri <- function(icd9) { 
  ahrq <- icd9cm_elixhauser_ahrq37(icd9)                                
  quan <- icd9cm_charlson_quan(icd9)
  
  data.frame(
    row.names = icd9,
    chf = ahrq$chf,
    cvd = quan$cvd,
    dm = ahrq$dm | ahrq$dmcx,
    ischemia = grepl('^D41[01234]', icd9),
    renlfail = ahrq$renlfail)
}

#' Convert ICD-9-CM code list to dataframe
#' 
#' \code{melt_icd9list} uses \code{\link{ddply}} to melt a column of comma-separated ICD-9-CM
#' codes into a series of rows, one for each code.
#' 
#' @param df a data frame with at least two columns, specified as \code{idvar}
#'   and \code{icd9var}.
#' @param idvar string with name of ID variable within \code{df} (defaults to "id")
#' @param icd9var string with name of ICD code variable within \code{df} (defaults to "icd9cm")
#' @param .progress passed to \code{\link{ddply}}
#' @param .parallel passed to \code{\link{ddply}}
#' @param .paropts passed to \code{\link{ddply}}
#' @return a dataframe with two columns, \code{idvar} and \code{"icd9cm"}
#' @importFrom plyr ddply
#' @examples
#' cases <- data.frame(id=c(1,2), icd9list=c('162.4,070.30,155.0,401.9','996.52,E878.8,V45.86'))
#' melt_icd9list(cases, "id", "icd9list")
#' @export
melt_icd9list <- function(df, idvar="id", icd9var="icd9cm",
    .progress="none", .parallel=FALSE, .paropts=NULL) {
    ddply(df, idvar, function(x) { 
        data.frame(
            icd9cm=unlist(strsplit(
                paste(gsub('.','',x[,icd9var],fixed=TRUE),collapse=','),',',fixed=TRUE))) 
    }, .progress=.progress, .parallel=.parallel, .paropts=.paropts)
}

#' Merge ICD-9-CM diagnostic and procedure codes
#' 
#' Merges a dataframe containing ICD-9-CM diagostic codes with a dataframe containing ICD-9 procedure codes
#' Diagnostic codes are prefixed with 'D', while procedure codes are prefixed with 'P'
#' 
#' @param dx_df a data frame with at least two columns, specified as \code{idvar}
#'   and \code{icd9dxvar}, where the values are ICD-9 diagnostic codes
#' @param proc_df a data frame with at least two columns, specified as \code{idvar}
#'   and \code{icd9pvar}, where the values are ICD-9 procedure codes
#' @param icd9dxvar name of icd9 diagnostic code column, default "icd9cm"
#' @param icd9pvar name of icd9 procedure code column, default "icd9cm"
#' @return a merged dataframe with common columns and \code{"icd9cm"}
#' @examples
#' cases <- data.frame(id=c(1,2),
#'                     icd9dxlist=c('162.4,070.30,155.0,401.9','996.52,E878.8,V45.86'), 
#'                     icd9plist=c('38.16','38.42'))
#' dx_df <- melt_icd9list(cases, "id", "icd9dxlist")
#' proc_df <- melt_icd9list(cases, "id", "icd9plist")
#' merge_icd9_dx_and_procs(dx_df, proc_df)
#' @export
merge_icd9_dx_and_procs <- function(dx_df, proc_df, icd9dxvar="icd9cm", icd9pvar="icd9cm") {
  if (!is.null(dx_df)) {
    dx_df$icd9cm <- paste('D',dx_df[,icd9dxvar],sep='')
  }
  if (!is.null(proc_df)) {
    proc_df$icd9cm <- paste('P',proc_df[,icd9pvar],sep='')
  }
  if (!is.null(proc_df) & !is.null(dx_df)) {
    commonCols <- intersect(names(dx_df),names(proc_df))
    df <- rbind(dx_df[,commonCols],proc_df[,commonCols])
  } else if (is.null(proc_df)) {
    df <- dx_df
  } else if (is.null(dx_df)) {
    df <- proc_df
  }
  df$icd9cm <- factor(df$icd9cm)
  df
}

#' Generate a comorbidity dataframe
#' 
#' Merges a given DF of IDs and ICD-9-CM codes to one of the ICD9CM maps, 
#' removes redundant comorbidities, and returns a dataframe.
#' 
#' Redundancy rules:
#' * If "tumor" and "mets", only "mets" will be returned.
#' * If "htn" and "htncx", only "htncx" will be returned.
#' * If "dm" and "dmcx", only "dmcx" will be returned.
#' * If "liver" and "modliver", only "modliver" will be returned.
#' 
#' Van Walraven has a modification adopted here where the following "dmcx" codes
#' are downgraded to "dm" if the specific DM complication is separately coded:
#' * D2(49|50)4x is DM w renal
#' * D2(49|50)6x is DM w neuro
#' * D2(49|50)7x is DM w PVD
#' 
#' Cases without any comorbidities will not appear in the returned data
#' frame.
#' 
#' @param df a data frame with at least two columns, specified as \code{idvar} 
#'   and \code{icd9var}.
#' @param idvar string with name of ID variable within \code{df} (defaults to 
#'   "id")
#' @param icd9var string with name of ICD code variable within \code{df} 
#'   (defaults to \code{icd9cm})
#' @param icd9mapfn Function to generate comorbidity data frame from ICD-9 codes
#'   (defaults to \code{\link{icd9cm_charlson_quan}})
#' @param .progress passed to \code{\link{ddply}}
#' @param .parallel passed to \code{\link{ddply}}
#' @param .paropts passed to \code{\link{ddply}}
#' @return a dataframe with column \code{idvar} and a logical column for each comorbidity
#' @importFrom plyr ddply
#' @examples
#' cases <- data.frame(id=c(1,1,1,2,2,2,2,2), 
#'          icd9cm=c("D20206","D24220","D4439","D5064","DE8788","D40403","D1960","D1958"))
#' generate_comorbidity_df(cases)
#' # generate categories for patients in the \code{\link{vt_inp_sample}}
#' generate_comorbidity_df(vt_inp_sample)
#' # in this example, D25071 is reduced to "dm" from "dmcx" because D4439 already codes perivasc
#' # also, D20206 "tumor" and D1970 "mets" lead to just "mets"
#' # D25001 and D25040 are just "dmcx"
#' # D45621 and D570 are just "modliver"
#' cases <- data.frame(id=c(1,1,1,1,2,2,2,2),
#'   icd9cm=c("D1970","D20206","D25071","D4439","D25001","D25040","D45621","D570"))
#' generate_comorbidity_df(cases)
#' @export
generate_comorbidity_df <- function(df, idvar="id", icd9var="icd9cm", 
                                 icd9mapfn=icd9cm_charlson_quan,
                                 .progress="none", .parallel=FALSE, .paropts=NULL) {

  # create a unique vector for icd9 codes
  uniq_icd9map <-
    if (is.factor(df[,icd9var])) levels(df[,icd9var]) else unique(as.character(df[,icd9var]))

  icd9map <- icd9mapfn(uniq_icd9map)
  cases_merged <- merge(df[,c(idvar,icd9var)], icd9map, by.x=icd9var, by.y="row.names", all.x=T)
  
  ddply(cases_merged, c(idvar), 
        function(df) {
          rows <- c()
          # preprocess df with van Walraven step
          if ((all(c("renal","dm","dmcx") %in% names(df)) & any(df$renal)) | 
                (all(c("renlfail","dm","dmcx") %in% names(df)) & any(df$renlfail))) {
            rows <- grep('^D2(49|50)4', df[,icd9var])
            if (length(rows)) 
              df[rows, c("dm","dmcx")] <- c(T,F)
          }
          if ((all(c("neuro","dm","dmcx") %in% names(df)) & any(df$neuro)) | 
                (all(c("cvd","dm","dmcx") %in% names(df)) & any(df$cvd))) {
            rows <- grep('^D2(49|50)6', df[,icd9var])
            if (length(rows)) 
              df[rows, c("dm","dmcx")] <- c(T,F)
          }
          if (all(c("perivasc","dm","dmcx") %in% names(df)) & any(df$perivasc)) {
            rows <- grep('^D2(49|50)7', df[,icd9var])
            if (length(rows)) 
              df[rows, c("dm","dmcx")] <- c(T,F)
          }

          # now merge rows together with any
          out <- data.frame(lapply(df[,names(icd9map)], any))
          
          # postprocess: eliminate redundancy
          if (all(c("mets","tumor") %in% names(out))) {
            rows <- which(out$mets & out$tumor)
            if (length(rows))
              out$tumor[rows] <- F
          }
          if (all(c("htn","htncx") %in% names(out))) {
            rows <- which(out$htn & out$htncx)
            if (length(rows))
              out$htn[rows] <- F
          }
          if (all(c("dm","dmcx") %in% names(out))) {
            rows <- which(out$dm & out$dmcx)
            if (length(rows))
              out$dm[rows] <- F
          }
          if (all(c("liver","modliver") %in% names(out))) {
            rows <- which(out$liver & out$modliver)
            if (length(rows))
              out$liver[rows] <- F
          }
          out
        },
        .progress=.progress, .parallel=.parallel, .paropts=.paropts)
}

#' Calculate the Charlson Comorbidity Index
#' 
#' \code{generate_charlson_index_df} merges a data frame of Charlson 
#' comorbidities with \code{\link{charlson_weights}} and sums the results per 
#' patient.
#' 
#' @param df a data frame with ID column \code{idvar} and logical columns for each 
#' comorbidity, such as that generated by \code{\link{generate_comorbidity_df}}
#' @param idvar string with name of ID variable within \code{df}
#' @param weights defaults to \code{\link{charlson_weights}}
#' @return a dataframe with two columns, \code{idvar} and \code{"index"}
#' @seealso \code{\link{generate_comorbidity_df}},
#'   \code{\link{charlson_weights}}, \code{\link{charlson_weights_orig}}
#' @importFrom plyr empty
#' @examples
#' # calculate Charlson Comorbidity Index for all patients in the \code{\link{vt_inp_sample}}
#' data(vt_inp_sample)
#' generate_charlson_index_df(generate_comorbidity_df(vt_inp_sample))
#' @export
generate_charlson_index_df <- function(df, idvar="id", weights=medicalrisk::charlson_weights) {
    if (empty(df)) {
      return(df)
    }
        
    if (!exists(idvar,where=df)) {
      stop("Data frame must have variable specified in idvar")
    }
    
    missing_cols <- setdiff(names(weights), names(df))
    if (length(missing_cols)) {
      stop(sprintf("Missing Charlson columns: %s", paste(missing_cols, collapse=",")))
    }

    # to generate index score, multiply weights by logical columns, then take row sum of that
    rv <- data.frame(id=df[,idvar], 
                     index=rowSums(df[,names(weights)] * weights))
    names(rv)[1] <- idvar
    rv
}
