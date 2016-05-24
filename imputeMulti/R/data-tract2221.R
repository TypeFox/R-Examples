

#' @title Observational data on individuals living in census tract 2221
#' @description A dataset containing attributes of 3974 individuals living in census tract 2221
#' in Los Angeles County, CA. Data comes from the 5-year American Community Survey with end year
#' 2014. Missing values have been inserted.
#' @format A \code{data.frame} with 3974 rows and 10 variables. All variables are of class \code{factor}:
#' \describe{
#'    \item{age}{The individual's age coded in roughly 5 year age buckets.}
#'    \item{gender}{The indiviudals gender -- Male, Female}
#'    \item{marital_status}{The individuals marital status. Takes one of 5 levels: 
#'      \code{never_mar} never married; \code{married} married; \code{mar_apart} married but living apart;
#'      \code{divorced} divorced; and \code{widowed} widowed}
#'    \item{edu_attain}{The individual's educational attainment. Takes one of 7 levels: 
#'      \code{lt_hs} less than high school; \code{some_hs} completed some high school but did not graduate;
#'      \code{hs_grad} high school graduate; \code{some_col} completed some college but did not graduate;
#'      \code{assoc_dec} completed an associates degree; \code{ba_deg} obtained a bachelors degree;
#'      \code{grad_deg} obtained a graduate or professional degree}
#'    \item{nativity}{The individual's nativity status. Takes one of 4 values: \code{born_state_residence}
#'      born in the state of residence; \code{born_other_state} born in another US state; \code{born_out_us}
#'      a US citizen born outside the US; \code{foreigner} foreign born}
#'    \item{pov_status}{The individual's poverty status in the past year. Takes one of 2 levels: 
#'      \code{below_pov_level} below the poverty level; \code{at_above_pov_level} at or above the poverty level}
#'    \item{geog_mobility}{The individual's geographic mobility in the last year. Takes one of 5 values:
#'      \code{same house} lived in the same house; \code{same county} moved within the same county;
#'      \code{same state} moved within the same state; \code{same state} moved from a different county
#'      within the same state; \code{diff state} moved from a different state; \code{moved from abroad}
#'      moved from another country}
#'    \item{ind_income}{The individual's annual income. Takes one of 9 levels: \code{no_income} no income;
#'      \code{1_lt10k} income <$10,000; \code{10k_lt15k} $10000-$14999; \code{15k_lt25k} $15000-$24999;
#'      \code{25k_lt35k} $25000-$34999; \code{35k_lt50k} $35000-$49999; \code{50k_lt65k} $50000-$64999;
#'      \code{65k_lt75k} $65000-$74999; \code{gt75k} $75000+}
#'    \item{race}{The individual's ethnicity.}
#' }
"tract2221"