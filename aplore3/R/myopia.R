#' MYOPIA data
#' 
#' myopia dataset.
#' 
#' @format A data.frame with 618 rows and 18 variables:
#' \describe{
#' \item{id}{Subject identifier (1-1503)}
#' \item{studyyear}{Year subject entered the study (Year)}
#' \item{myopic}{Myopia within the first five years of follow up (1: No, 2: Yes)}
#' \item{age}{Age at first visit (Years)}
#' \item{gender}{Gender (1: Male, 2: Female)}
#' \item{spheq}{Spherical Equivalent Refraction (diopter)}
#' \item{al}{Axial Length (mm)}
#' \item{acd}{Anterior Chamber Depth (mm)}
#' \item{lt}{Lens Thickness (mm)}
#' \item{vcd}{Vitreous Chamber Depth (mm)}
#' \item{sporthr}{How many hours per week outside of school the child spent
#' engaging in sports/outdoor activities (Hours per week)}
#' \item{readhr}{How many hours per week outside of school the child spent reading for pleasure (Hours per week)}
#' \item{comphr}{How many hours per week outside of school the child spent playing video/computer games or working on the computer (Hours per week)}
#' \item{studyhr}{How many hours per week outside of school the child spent reading or studying for school assignments (Hours per week)}
#' \item{tvhr}{How many hours per week outside of school the child spent watching television (Hours per week)}
#' \item{diopterhr}{Composite of near-work activities (Hours per week)}
#' \item{mommy}{Was the subject's mother myopic? (1: No, 2: Yes)}
#' \item{dadmy}{Was the subject's father myopic? (1: No, 2: Yes)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(myopia, n = 10)
#' summary(myopia)
"myopia"
