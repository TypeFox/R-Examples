utils::globalVariables(c("adv", "grade", "count", "id", "gradeP", "vals", "prof"))
#' Evaluation of educational policy tools
#' @description Make common tasks for educational evaluation easier to do!
#' @name eeptools
#' @details \tabular{ll}{
#' Package: \tab eeptools\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1 \cr
#' Date: \tab 2012-07-21\cr
#' License: \tab GPL-3 \cr
#' }
#' his package has a number of useful shortcuts for common tasks. It includes 
#' some themes for ggplot2 plots, processing arbitrary text files of data, 
#' calculating student characteristics, and finding thresholds within vectors. 
#' Future development work will include methods for tuning and evaluating early 
#' warning system models.
#' @author Jared E. Knowles
#' @note This package is still in beta and function names may change in the next release.
#' @docType package
#' @examples 
#' gender<-c("M","M","M","F","F","F")
#' statamode(gender)
#' statamode(gender[1:5])
#' 
#' missing_data<-c(NA,NA,NA)
#' max_mis(missing_data)
#' 
#' makenum(gender)
#' gender <- factor(gender)
#' defac(gender)
NULL


#' A dataframe of aggregate test scores for schools in a Midwest state.
#' @description This data comes from publicly available aggregated test scores 
#' of a large midwestern state. Each row represents scores for school A in grade X 
#' and then scores in school A and grade X+1. Additionally, some regression 
#' diagnostics and results from a predictive model of test scores in grade 
#' X+1 are included.
#' @format A data frame with 19985 observations on the following 16 variables.
#' \describe{
#'  \item{\code{district_id}}{a numeric vector}
#'  \item{\code{school_id}}{a numeric vector}
#'  \item{\code{subject}}{a factor with levels \code{math} \code{read} representing the subject of the test scores in the row}
#'  \item{\code{grade}}{a numeric vector}
#'  \item{\code{n1}}{a numeric vector for the count of students in the school and grade in t}
#'  \item{\code{ss1}}{a numeric vector for the scale score in t}
#'  \item{\code{n2}}{a numeric vector for the count of students in the school and grade in t+1}
#'  \item{\code{ss2}}{a numeric vector for the mean scale score in t+1}
#'  \item{\code{predicted}}{a numeric vector of the predicted ss2 for this observation}
#'  \item{\code{residuals}}{a numeric vector of residuals from the predicted ss2}
#'  \item{\code{resid_z}}{a numeric vector of standardized residuals}
#'  \item{\code{resid_t}}{a numeric vector of studentized residuals}
#'  \item{\code{cooks}}{a numeric vector of cooks D for the residuals}
#'  \item{\code{test_year}}{a numeric vector representing the year the test was taken}
#'  \item{\code{tprob}}{a numeric vector representing the probability of a residual appearing}
#'  \item{\code{flagged_t95}}{a numeric vector}
#' }
#' @details These data were fit with a statistical model by a large newspaper to 
#' investigate unusual gains in test scores. Fifty separate models were fit 
#' representing all unique combinations of grade,year, and subject
#' @examples 
#' data(midsch)
#' head(midsch)
"midsch"

#' Student Attributes from the Strategic Data Project Toolkit
#' @description A synthetic dataset of student attributes from the Strategic Data 
#' Project which includes records with errors to practice data cleaning and 
#' impelmenting business rules for consistency in data. 
#' @format A data frame with 87534 observations on the following 9 variables.
#' \describe{
#'   \item{\code{sid}}{a numeric vector of the unique student ID}
#'   \item{\code{school_year}}{a numeric vector of the school year}
#'   \item{\code{male}}{a numeric vector indicating 1 = male}
#'   \item{\code{race_ethnicity}}{a factor with levels \code{A} \code{B} \code{H} \code{M/O} \code{W}}
#'   \item{\code{birth_date}}{a numeric vector of the student birthdate}
#'   \item{\code{first_9th_school_year_reported}}{a numeric vector of the first year a student is reported in 9th grade}
#'   \item{\code{hs_diploma}}{a numeric vector}
#'   \item{\code{hs_diploma_type}}{a factor with levels \code{} \code{Alternative Diploma} \code{College Prep Diploma} \code{Standard Diploma}}
#'   \item{\code{hs_diploma_date}}{a factor with levels \code{} \code{12/2/2008} \code{12/21/2008} \code{4/14/2008} \code{4/18/2008} ...}
#'   }
#' @details This is the non-clean version of the data to allow for implementing 
#' business rules to clean data.
#' @source Available from the Strategic Data Project online at 
#' \url{http://sdp.cepr.harvard.edu/toolkit-effective-data-use}
#' @references Visit the Strategic Data Project online at:  \url{http://sdp.cepr.harvard.edu/}
#' @examples 
#' data(stuatt)
#' head(stuatt)
"stuatt"

#' A synthetic data set of K-12 student attributes.
#' @description A small dataset of synthetic data on K-12 students with 2700 
#' observations. 1200 individual students are represented, nested within 
#' 4 districts and 2 schools.
#' @format A data frame with 2700 observations on the following 32 variables.
#' \describe{
#'  \item{\code{X}}{a numeric vector}
#'  \item{\code{school}}{a numeric vector}
#'  \item{\code{stuid}}{a numeric vector}
#'  \item{\code{grade}}{a numeric vector}
#'  \item{\code{schid}}{a numeric vector}
#'  \item{\code{dist}}{a numeric vector}
#'  \item{\code{white}}{a numeric vector}
#'  \item{\code{black}}{a numeric vector}
#'  \item{\code{hisp}}{a numeric vector}
#'  \item{\code{indian}}{a numeric vector}
#'  \item{\code{asian}}{a numeric vector}
#'  \item{\code{econ}}{a numeric vector}
#'  \item{\code{female}}{a numeric vector}
#'  \item{\code{ell}}{a numeric vector}
#'  \item{\code{disab}}{a numeric vector}
#'  \item{\code{sch_fay}}{a numeric vector}
#'  \item{\code{dist_fay}}{a numeric vector}
#'  \item{\code{luck}}{a numeric vector}
#'  \item{\code{ability}}{a numeric vector}
#'  \item{\code{measerr}}{a numeric vector}
#'  \item{\code{teachq}}{a numeric vector}
#'  \item{\code{year}}{a numeric vector}
#'  \item{\code{attday}}{a numeric vector}
#'  \item{\code{schoolscore}}{a numeric vector}
#'  \item{\code{district}}{a numeric vector}
#'  \item{\code{schoolhigh}}{a numeric vector}
#'  \item{\code{schoolavg}}{a numeric vector}
#'  \item{\code{schoollow}}{a numeric vector}
#'  \item{\code{readSS}}{a numeric vector}
#'  \item{\code{mathSS}}{a numeric vector}
#'  \item{\code{proflvl}}{a factor with levels \code{advanced} \code{basic} \code{below basic} \code{proficient}}
#'  \item{\code{race}}{a factor with levels \code{A} \code{B} \code{H} \code{I} \code{W}}
#'  }
#' @details This data is synthetically generated to reflect student test scores 
#' and demographic attributes. 
#' @source The script to generate this synthetic dataset can be found and modified 
#' at \url{https://github.com/jknowles/r_tutorial_ed}
#' @examples 
#' data(stulevel)
#' head(stulevel)
"stulevel"
