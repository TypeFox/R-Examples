#' @name icu
#' @docType data
#' @title Intensive Care Unit study data
#'
#' @format A \code{data.frame} with
#' \eqn{200} observations (rows)
#' and \eqn{14} variables (columns).
#'
#' @details A sample of 200 subjects who were part of a
#' study on survival of patients admitted to an adult
#' intensive care unit (ICU).
#' \cr
#' The observed variable values were modified to
#' protect patient confidentiality.
#' \cr \cr
#' Columns are:
#' \describe{
#'  \item{ID}{Identification code.}
#'  \item{STA}{Vital status (\code{factor}):
#'   \describe{
#'    \item{0}{lived}
#'    \item{1}{died}}}
#'  \item{AGE}{Age (years).}
#'  \item{SEX}{Gender (\code{factor}):
#'   \describe{
#'    \item{0}{male}
#'    \item{1}{female}}}
#'  \item{RACE}{Race (\code{factor}):
#'   \describe{
#'    \item{1}{white}
#'    \item{2}{black}
#'    \item{3}{other}}}
#'  \item{SER}{Service, when admitted to ICU (\code{factor}):
#'   \describe{
#'    \item{0}{Medical}
#'    \item{1}{Surgical}}}
#' \item{CAN}{Cancer part of present problem? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{CRN}{Chronic renal failure? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{INF}{Infection probable
#'            when admitted to ICU? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{CPR}{Cardiopulmonary resuscitataion prior
#'            to ICU admission? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{SYS}{Systolic blood pressure (mmHG)
#'            when admitted to ICU.}
#' \item{HRA}{Heart rate when admitted to ICU.}
#' \item{PRE}{Previous admission to ICU
#'            within 6 months? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{TYP}{Type of admission (\code{factor}):
#'   \describe{
#'    \item{0}{elective}
#'    \item{1}{emergency}}}
#' \item{FRA}{Fracture present (long bone, multiple,
#'             neck, single area or hip)? (\code{factor}):
#'   \describe{
#'    \item{0}{no}
#'    \item{1}{yes}}}
#' \item{PO2}{pO2 from initial blood gases (\code{factor}):
#'   \describe{
#'    \item{0}{>60}
#'    \item{1}{<=60}}}
#' \item{PH}{pH from initial blood gases (\code{factor}):
#'   \describe{
#'    \item{0}{>=7.25}
#'    \item{1}{<7.25}}}
#' \item{PCO}{pCO2 from initial blood gases (\code{factor}):
#'   \describe{
#'    \item{0}{>=18}
#'    \item{1}{<18}}}
#' \item{CRE}{Creatinine from
#'            initial blood gases (\code{factor}):
#'   \describe{
#'    \item{0}{<=2}
#'    \item{1}{>2}}}
#' \item{LOC}{Level of consciousness
#'            when admitted to ICU (\code{factor}):
#'   \describe{
#'    \item{0}{no_coma}
#'    \item{1}{deep_stupor}
#'    \item{2}{coma}}}
#' }
#'
#' @keywords datasets
#'
#' @source
#' \href{ftp://ftp.wiley.com/public/sci_tech_med/logistic}{
#'       Wiley FTP}
#' @references
#' \bold{H&L 2nd ed.} Page 22, Section 1.6.1.
#'
#' Lemeshow S,  Teres D, Avrunin JS, Pastides H 1988.
#' Predicting the outcome of intensive care unit patients.
#' \emph{Journal of the American Statistical Association}.
#' \bold{83}(402):348--356.
#' \href{http://www.jstor.org.cuhsl.creighton.edu/stable/2288849}{
#'       JSTOR (free)}
#'
#' Lemeshow S, Teres D, Klar J, Avrunin JS,
#' Gehlbach SH , Rapoport John 1993.
#' Mortality Probability Models (MPM II) based on
#' an international cohort of intensive care unit patients.
#' \emph{Journal of the American Medical Association}.
#' \bold{270}(20):2478--2486.
#' \href{http://dx.doi.org/10.1001/jama.1993.03510200084037}{
#'       JAMA (paywall)}
#'
#' Lemeshow S, Le Gall J 1994.
#' Modeling the severity of illness of ICU patients:
#' a systems update.
#' \emph{Journal of the American Medical Association}.
#' \bold{272}(13):1049--1055.
#' \href{http://dx.doi.org/10.1001/jama.1994.03520130087038}{
#'       JAMA (paywall)}
NULL



