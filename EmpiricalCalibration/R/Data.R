# @file Data.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Incidence rate ratios from Self-Controlled Case Series
#'
#' @details
#' A dataset containing the incidence rate ratios (and standard errors) produced using a
#' Self-Controlled Case Series (SCCS) design. The outcome is upper GI bleeding, the drug of interest
#' (groundTruth = 1) is sertraline. Also included are 45 negative control drugs, for which we believe
#' there to be no causal relation with upper GI bleeding. We used a database of medical records from
#' general practices in the USA, the General Electric (GE) Centricity database, which contains data on
#' 11.2 million subjects. We restricted on study period (start of 1990 through November 2003), age
#' requirements (18 years or older), available time prior to event (180 days), and risk definition
#' window (30 days following the prescription). Time 30 days prior to the first prescription was
#' removed to account for possible contra-indications. Cases of upper GI bleeding were identified on
#' the basis of the occurrence of ICD-9 diagnosis codes in the problem list. These codes pertain to
#' esophageal, gastric, duodenal, peptic, and gastrojejunal ulceration, perforation, and hemorrhage,
#' as well as gastritis and non-specific gastrointestinal hemorrhage. For more information on this set
#' see Schuemie et al (2014).
#'
#' @docType data
#' @keywords datasets
#' @name sccs
#' @usage
#' data(sccs)
#' @format
#' A data frame with 46 rows and 4 variables: \describe{ \item{drugName}{Name of the drug}
#' \item{groundTruth}{Whether the drug is a positive (1) or negative (0) control} \item{logRr}{The log
#' of the incidence rate ratio} \item{seLogRr}{The standard error of the log of the incidence rate
#' ratio} }
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
NULL

#' Odds ratios from a case-control design
#'
#' @details
#' A dataset containing the odds ratios (and standard errors) produced using a case-control design.
#' The outcome is upper GI bleeding, the drug of interest (groundTruth = 1) is sertraline. Also
#' included are 46 negative control drugs, for which we believe there to be no causal relation with
#' upper GI bleeding. We used a database of medical records from general practices in the USA, the
#' General Electric (GE) Centricity database, which contains data on 11.2 million subjects. We
#' restricted on study period (start of 1990 through November 2003), age requirements (18 years or
#' older), available time prior to event (180 days), number of controls per case (6), and risk
#' definition window (30 days following the prescription). Controls were matched on age and sex. Cases
#' of upper GI bleeding were identified on the basis of the occurrence of ICD-9 diagnosis codes in the
#' problem list. These codes pertain to esophageal, gastric, duodenal, peptic, and gastrojejunal
#' ulceration, perforation, and hemorrhage, as well as gastritis and non-specific gastrointestinal
#' hemorrhage. For more information on this set see Schuemie et al (2014).
#'
#' @docType data
#' @keywords datasets
#' @name caseControl
#' @usage
#' data(caseControl)
#' @format
#' A data frame with 47 rows and 4 variables: \describe{ \item{drugName}{Name of the drug}
#' \item{groundTruth}{Whether the drug is a positive (1) or negative (0) control} \item{logRr}{The log
#' of the incidence rate ratio} \item{seLogRr}{The standard error of the log of the incidence rate
#' ratio} }
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
NULL

#' Relative risks from a new-user cohort design
#'
#' @details
#' A dataset containing the relative risks (and standard errors) produced using a new-user cohort
#' design. The outcome is acute liver injury, the drug of interest (groundTruth = 1) is Isoniazid Also
#' included are 30 negative control drugs, for which we believe there to be no causal relation with
#' acute liver injury. We used the Thomson MarketScan Medicare Supplemental Beneficiaries database,
#' which contains data on 4.6 million subjects. We selected two groups (cohorts): (1) all subjects
#' exposed to isoniazid and (2) all subjects having the ailment for which isoniazid is indicated, in
#' this case tuberculosis, and having received at least one drug that is not known to cause acute
#' liver injury. We removed all subjects who belonged to both groups and subjects for which less than
#' 180 days of observation time was available prior to their first exposure to the drug in question.
#' Acute liver injury was identified on the basis of the occurrence of ICD-9-based diagnosis codes
#' from inpatient and outpatient medical claims and was defined broadly on the basis of codes
#' associated with hepatic dysfunction, as have been used in prior observational database studies. The
#' time at risk was defined as the length of exposure + 30 days, and we determined whether subjects
#' experienced an acute liver injury during their time at risk. Using propensity score stratification,
#' the cohorts were divided over 20 strata, and an odds ratio over all strata was computed using a
#' Mantel-Haenszel test. The propensity score was estimated using Bayesian logistic regression using
#' all available drug, condition, and procedure covariates occurring in the 180 days prior to first
#' exposure, in addition to age, sex, calendar year of first exposure, Charlson index, number of
#' drugs, number of visit days, and number of procedures. For more information on this set see
#' Schuemie et al (2014).
#'
#' @docType data
#' @keywords datasets
#' @name cohortMethod
#' @usage
#' data(cohortMethod)
#' @format
#' A data frame with 31 rows and 4 variables: \describe{ \item{drugName}{Name of the drug}
#' \item{groundTruth}{Whether the drug is a positive (1) or negative (0) control} \item{logRr}{The log
#' of the incidence rate ratio} \item{seLogRr}{The standard error of the log of the incidence rate
#' ratio} }
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
NULL
