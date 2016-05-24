#' M36NT2PM.dat - data file in the NAEP Primer data set
#'
#' The NAEP Primer is designed to simplify access to the NAEP database and make
#' its technologies more user-friendly in NAEP data analysis. The Primer data
#' contain a random mini-sample of real data from the 2005 grade 8 mathematics
#' assessment that have been approved for public use. Only public schools are
#' included in this subsample that contains selected variables for about 10 percent
#' of the schools and students in this assessment. All students who participated
#' in NAEP in the selected public schools are included. The mini-sample contains
#' assessment data for students (assessed and excluded), as well as questionnaire
#' data for teachers and schools. It has one combined data file with student data,
#' teacher data, and school data at the student level. This mini-sample is not
#' sufficient to make state comparisons. In addition, to ensure confidentiality,
#' no state, school, or student identifiers are included (see The NAEP Primer
#' (NCES 2011463) for further information). 
#'
#' Together M36NT2PM.dat and M36NC2PM.dat are the data for the NAEP Primer--for
#' students and schools, respectively. The files M36NT2PM.fr2 M36NC2PM.fr2
#' are the layout files for the same data.
#'
#' Traditionally the NAEP Primer would be a single file but for exposition purposes we have split it up into a 
#' M36NT2PM.dat contains 300 variables and 17606 rows, and has the majority of student, teacher, and school data.
#' M36NC2PM.dat contains 2 variables and 4092 rows (one per school), and has two school variables. 
#' The two school variables on M36NC2PM are:
#'
#' \itemize{
#'  \item{"SSCRPSU"}{Scrambled PSU and school code (used for linking)}
#'  \item{"C052601"}{Percent enrolled in math for remediation}
#' }
#'
#' @aliases M36NC2PM.dat M36NT2PM.fr2 M36NC2PM.fr2 
#' @format The .dat files are in a fixed width dat file. The format is described by the .fr2 files.
#' @docType data
#' @keywords datasets
#' @name M36NT2PM.dat
#' @references \url{https://nces.ed.gov/pubsearch/pubsinfo.asp?pubid=2011463}
NULL

