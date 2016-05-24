

#'Growth of Dutch boys
#'
#'Height, weight, head circumference and puberty of 7482 Dutch boys.
#'
#'The complete sample of cross-sectional data from boys 0-21 years used to
#'construct the Dutch growth references 1997. Variables \code{gen} and
#'\code{phb} are ordered factors. \code{reg} is a factor. Note: A 10\% sample
#'from this data is available in data set \code{boys} in the \code{mice}
#'package.
#'
#'@name boys7482
#'@docType data
#'@format A data frame with 7482 rows on the following 9 variables: 
#'\describe{
#'\item{age}{Decimal age (0-21 years)} 
#'\item{hgt}{Height (cm)}
#'\item{wgt}{Weight (kg)} 
#'\item{bmi}{Body mass index}
#'\item{hc}{Head circumference (cm)} 
#'\item{gen}{Genital Tanner stage (G1-G5)} 
#'\item{phb}{Pubic hair (Tanner P1-P6)}
#'\item{tv}{Testicular volume (ml)} 
#'\item{reg}{Region (north, east, west, south, city)} 
#'}
#'@author Stef van Buuren, 2012
#'@source Fredriks, A.M,, van Buuren, S., Burgmeijer, R.J., Meulmeester JF,
#'Beuker, R.J., Brugman, E., Roede, M.J., Verloove-Vanhorick, S.P., Wit, J.M.
#'(2000) Continuing positive secular growth change in The Netherlands
#'1955-1997.  \emph{Pediatric Research}, \bold{47}, 316-323. 
#'
#'Fredriks, A.M., van Buuren, S., Wit, J.M., Verloove-Vanhorick, S.P. (2000).
#'Body index measurements in 1996-7 compared with 1980.  \emph{Archives of
#'Disease in Childhood}, \bold{82}, 107-112. 
#'@keywords datasets
NULL





#'Reference tables from CDC 2000
#'
#'Reference tables from CDC 2000
#'
#'The models were fitted by the LMS model. Parameters are stored as type
#'\code{LMS}. Tabulated values are point ages.
#'
#'The naming conventions are as follows: \describe{
#'\item{list("cdc.hgt")}{Combined length/height (cm) for Age, 0-20 years.
#'Measures <2 years apply to length (lying), while ages >= 2 years apply to
#'height, or stature (standing).} 
#'\item{list("cdc.wgt")}{Weight (kg) for Age, 0-20 years.} 
#'\item{list("cdc.bmi")}{Body Mass Index (kg/m2) for Age, 2-20 years.}}
#'
#'@name References CDC
#'@aliases cdc.hgt cdc.wgt cdc.bmi
#'@docType data
#'@format A data frame with seven variables: \describe{
#'\item{list("pop")}{Study Population} \item{list("sub")}{Subpopulation}
#'\item{list("sex")}{Sex (M,F)} \item{list("x")}{Decimal age (0-5 years)}
#'\item{list("L")}{Lambda (skewness) curve} \item{list("M")}{Median curve}
#'\item{list("S")}{Coefficient of Variation curve} }
#'@seealso \code{\link{nl4.wgt}}, \code{\link{nl4.hgt}}, \code{\link{nl4.bmi}},
#'\code{\link{who.wgt}}
#'@source Kuczmarski RJ, Ogden CL, Guo SS, Grummer-Strawn LM, Flegal KM, Mei Z,
#'Wei R, Curtin LR, Roche AF, Johnson CL.  2000 CDC growth charts for the
#'United States: methods and development.  \emph{Vital Health Stat}, 2002,
#'\bold{11}, 246, 1-190.
#'@keywords datasets
NULL





#'Reference tables from Third Dutch Growth Study 1980
#'
#'Reference table from the Third Dutch Growth Study 1980
#'
#'The model was fitted by the LMS model. Parameters are stored as type
#'\code{LMS}. Tabulated values are point ages.
#'
#'Height follows a normal distribution, with all lambda parameters set equal to
#'1. The standard deviation (in cm) is obtained as \code{S*M}.
#'
#'The naming conventions are as follows: \describe{
#'\item{list("nl4.hgt")}{Length/Height (cm) for Age}
#'\item{list("nl4.wgt")}{Weight (kg) for Age} \item{list("nl4.wfh")}{Weight
#'(kg) for Height (cm)} \item{list("nl4.bmi")}{Head circumference (cm) for Age}
#'\item{list("nl4.lgl")}{Leg Length (cm) for Age} \item{list("nl4.hip")}{Hip
#'circumference (cm) for Age} \item{list("nl4.wst")}{Waist circumference (cm)
#'for Age} \item{list("nl4.whr")}{Waist/Hip ratio for Age}
#'\item{list("nl4.sit")}{Sitting Height for Age} \item{list("nl4.shh")}{Sitting
#'Height/Height ratio for Age} }
#'
#'@name References NL3
#'@aliases nl3.bmi
#'@docType data
#'@format A data frame with seven variables: \describe{
#'\item{list("pop")}{Study Population} \item{list("sub")}{Subpopulation, e.g.
#'ethnicity or age group (for \code{nl4.wfh})} \item{list("sex")}{Sex (M,F)}
#'\item{list("x")}{Decimal age (0-21 years) or Height (for \code{nl4.wfh})}
#'\item{list("L")}{Lambda (skewness) curve} \item{list("M")}{Median curve}
#'\item{list("S")}{Coefficient of Variation curve} }
#'@seealso \code{\link{cdc.wgt}}, \code{\link{who.wgt}}
#'@source Fredriks, A.M,, van Buuren, S., Burgmeijer, R.J., Meulmeester JF,
#'Beuker, R.J., Brugman, E., Roede, M.J., Verloove-Vanhorick, S.P., Wit, J.M.
#'(2000) Continuing positive secular growth change in The Netherlands
#'1955-1997.  \emph{Pediatric Research}, \bold{47}, 316-323. 
#'
#'Fredriks, A.M., van Buuren, S., Wit, J.M., Verloove-Vanhorick, S.P. (2000).
#'Body index measurements in 1996-7 compared with 1980.  \emph{Archives of
#'Disease in Childhood}, \bold{82}, 107-112. 
#'
#'@keywords datasets
NULL





#'Reference tables from Fourth Dutch Growth Study 1997
#'
#'Reference table from the Fourth Dutch Growth Study 1997
#'
#'The model was fitted by the LMS model. Parameters are stored as type
#'\code{LMS}. Tabulated values are point ages.
#'
#'Height follows a normal distribution, with all lambda parameters set equal to
#'1. The standard deviation (in cm) is obtained as \code{S*M}.
#'
#'The naming conventions are as follows: \describe{
#'\item{list("nl4.hgt")}{Length/Height (cm) for Age}
#'\item{list("nl4.wgt")}{Weight (kg) for Age} \item{list("nl4.wfh")}{Weight
#'(kg) for Height (cm)} \item{list("nl4.bmi")}{Head circumference (cm) for Age}
#'\item{list("nl4.lgl")}{Leg Length (cm) for Age} \item{list("nl4.hip")}{Hip
#'circumference (cm) for Age} \item{list("nl4.wst")}{Waist circumference (cm)
#'for Age} \item{list("nl4.whr")}{Waist/Hip ratio for Age}
#'\item{list("nl4.sit")}{Sitting Height for Age} \item{list("nl4.shh")}{Sitting
#'Height/Height ratio for Age} }
#'
#'@name References NL4
#'@aliases nl4.hgt nl4.wgt nl4.wfh nl4.bmi nl4.hdc nl4.lgl nl4.hip nl4.wst
#'nl4.whr nl4.sit nl4.shh
#'@docType data
#'@format A data frame with seven variables: \describe{
#'\item{list("pop")}{Study Population} \item{list("sub")}{Subpopulation, e.g.
#'ethnicity or age group (for \code{nl4.wfh})} \item{list("sex")}{Sex (M,F)}
#'\item{list("x")}{Decimal age (0-21 years) or Height (for \code{nl4.wfh})}
#'\item{list("L")}{Lambda (skewness) curve} \item{list("M")}{Median curve}
#'\item{list("S")}{Coefficient of Variation curve} }
#'@seealso \code{\link{cdc.wgt}}, \code{\link{who.wgt}}
#'@source Fredriks, A.M,, van Buuren, S., Burgmeijer, R.J., Meulmeester JF,
#'Beuker, R.J., Brugman, E., Roede, M.J., Verloove-Vanhorick, S.P., Wit, J.M.
#'(2000) Continuing positive secular growth change in The Netherlands
#'1955-1997.  \emph{Pediatric Research}, \bold{47}, 316-323. 
#'
#'Fredriks, A.M., van Buuren, S., Wit, J.M., Verloove-Vanhorick, S.P. (2000).
#'Body index measurements in 1996-7 compared with 1980.  \emph{Archives of
#'Disease in Childhood}, \bold{82}, 107-112. 
#'
#'@keywords datasets
NULL





#'References WHO
#'
#'Reference tables, combined from the WHO Multicentre Growth Reference 
#'Study (MGRS) (ages 0-5 years) and the WHO 2007 reference (5-19 years).
#'
#'The data were fitted by the LMS model. Parameters are stored as type
#'\code{LMS}. Tabulated values are point ages.
#'
#'The naming conventions are as follows: 
#'\describe{
#'\item{who.hgt}{Length (cm, 0-2 Yrs) or height (cm, 2-19 years)} 
#'\item{who.wgt}{Weight (kg) for age (0-10 years)}
#'\item{who.bmi}{BMI (kg/m^2) for age (0-19 years)}
#'\item{who.wfh}{Weight (kg) for height (65-120 cm)}
#'\item{who.wfl}{Weight (kg) for length (45-110 cm)}
#'}
#'
#'@name References WHO
#'@aliases who.wgt who.hgt who.bmi who.wfh who.wfl
#'@docType data
#'@format A data frame with seven variables: 
#'\describe{
#'\item{pop}{Study Population (always \code{"who"})} 
#'\item{sub}{Subpopulation (always \code{"N"})}
#'\item{sex}{Sex (M, F)}
#'\item{x}{Decimal age, height (cm) or length(cm)}
#'\item{L}{Lambda (skewness) curve}
#'\item{M}{Median curve}
#'\item{S}{Coefficient of variation}
#'}
#'@seealso \code{\link{nl4.wgt}}, \code{\link{cdc.wgt}}, 
#'\url{http://www.who.int/childgrowth/mgrs/en/},
#'\url{http://www.who.int/growthref/en/} 
#'@source 
#'WHO Multicentre Growth Reference Study Group.  WHO Child Growth
#'Standards based on length/height, weight and age.  \emph{Acta Paediatr},
#'Suppl. 2006, 450, 76-85.
#'
#'de Onis M, Onyango AW, Borghi E, Siyam A, Nishida C, Siekmann J.
#'Development of a WHO growth reference for school-aged children and adolescents
#'\emph{Bulletin of the World Health Organization}, 2007;85:660-7.
#'@keywords datasets
NULL



