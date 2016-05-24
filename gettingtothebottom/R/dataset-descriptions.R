#' Baltimore Youth Indicators - 2010 and 2011
#'
#' Data from the 2011 release of Baltimore's Community Health Profiles and Healthy Baltimore 2015 Safe Homes and Families indicators. Data provided by the Maryland Department of Health and Mental Hygiene, Maryland Department of the Environment, Baltimore Substance Abuse Systems (BSAS), the Baltimore City Health Department, and the United States Bureau of the Census
#'
#' @source Data obtained from the Baltimore Neighborhood Indicators Alliance using the Children and Family Health & Well-Being and Education indicators.
#' \url{http://www.bniajfi.org/data_downloads}
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 55 rows and 79 variables
#' @name baltimoreyouth
#' 
#' @param CSA2010 Name of Community Statistical Areas in 2010
#' @param eattendXX Number of Students Ever Enrolled 1st - 5th Grade
#' @param mattendXX Number of Students Ever Enrolled 6th - 8th Grade
#' @param heattendXX Number of Students Ever Enrolled 9th - 12th Grade
#' @param eenrolXX Number of Students Enrolled in 1st - 5th Grade
#' @param menrolXX Number of Students Enrolled in 6th - 8th Grade
#' @param henrolXX Number of Students Enrolled in 9th - 12th Grade
#' @param aastudXX Percent of Students that are African American
#' @param wstudXX  Percent of Students that are White (non-Hispanic)
#' @param hstudXX	Percent of Students that are Hispanic
#' @param abseXX	Percent of 1st-5th Grade Students that are Chronically Absent (Missing at least 20 days)
#' @param absmdXX	Percent of 6th-8th Grade Students that are Chronically Absent (Missing at least 20 days)
#' @param abshsXX	Percent of 9th-12th Grade Students that are Chronically Absent (Missing at least 20 days)
#' @param suspXX	Percentage of Students Suspended or Expelled During School Year
#' @param farmsXX	Percentage of Students Receiving Free or Reduced Meals
#' @param spedXX	Percentage of Students Enrolled in Special Education Programs
#' @param readyXX	Kindergarten School Readiness
#' @param 3mathXX	Percentage of 3rd Grade Students Passing MSA Math
#' @param 3readXX	Percentage of 3rd Grade Students Passing MSA Reading
#' @param 5mathXX	Percentage of 5th Grade Students Passing MSA Math
#' @param 5readXX	Percentage of 5th Grade Students Passing MSA Reading
#' @param 8mathXX	Percentage of 8th Grade Students Passing MSA Math
#' @param 8readXX	Percentage of 8th Grade Students Passing MSA Reading
#' @param hsaengXX Percentage of Students Passing H.S.A. English 
#' @param hsabioXX Percentage of Students Passing H.S.A. Biology
#' @param hsagovXX Percentage of Students Passing H.S.A. Government
#' @param hsaalgXX Percentage of Students Passing H.S.A. Algebra
#' @param dropXX High School Dropout/Withdrawl Rate
#' @param complXX High School Completion Rate
#' @param sclswXX Percent of Students Switching Schools within School Year
#' @param sclempXX Percentage of Population aged 16-19 in School and/or Employed
#' @param teenbirXX Teen Birth Rate per 1,000 Females (aged 15-19)
#' @param termbirXX Percent of Births Delivered at Term (37-42 Weeks)
#' @param birthwtXX Percent of Babies Born with a Satisfactory Birth Weight
#' @param prenatalXX Percent of Births Where the Mother Received Early Prenatal Care (First Trimester)
#' @param leadtestXX Number of Children (aged 0-6) Tested for Elevated Blood Lead Levels
#' @param ebllXX Percent of Children (aged 0-6) with Elevated Blood Lead Levels
#' @param leadvXX Percent of Lead Violations per 1,000 Residential Units
#' @param tanfXX Percent of Families Receiving TANF
#' @param liquorXX Liquor Outlet density (per 1,000 Residents)
#' @param fastfdXX Fast Food Outlet Density (per 1,000 Residents)
#' @param LifeExpXX Life Expectancy
#' @param mort1_XX Mortality by Age (Less than 1 year old)
#' @param mort14_XX Mortality by Age (1-14 years old)
#' @param mort24_XX Mortality by Age (15-24 years old)
#' @param mort44_XX Mortality by Age (25-44 years old)
#' @param mort64_XX Mortality by Age (45-64 years old)
#' @param mort84_XX Mortality by Age (65-84 years old)
#' @param mort85_XX Mortality by Age (85 and over)
#' 
NULL

#' Movie ratings and budget database derived from data from IMDB.com
#'
#' A dataset containing movie ratings and budget data for 5,183 movies.
#'
#' @source Data obtained from Hadley Wickham.  The data in this package contains only those movies not exceeding 400 minutes in length and those with known total budgets.
#' \url{http://had.co.nz/data/movies/}.
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 5183 rows and 24 variables
#' @name moviebudgets
#' 
#' @param title Title of the movie
#' @param year Year the movie was released
#' @param budget Total budget (if known) in U.S. dollars
#' @param length Length of movie (in minutes)
#' @param rating Average IMDB user rating
#' @param votes Number of IMDB users who rated the movie
#' @param r1-10 Distribution of votes for each rating, to mid point of nearest decile: 0 = no votes, 4.5 = 1-9 percent votes, 14.5 = 11-19 percent of votes, etc. Due to rounding errors these may not sum to 100.
#' @param mpaa MPAA rating
#' @param genre Binary variables indicating whether movie belongs to any of the following genres: action, animation, comedy, drama, documentary, romance, short
#' 
NULL

#' Movie ratings database derived from data from IMDB.com
#'
#' A dataset containing movie ratings for 58,771 movies.
#'
#' @source Data obtained from Hadley Wickham.  The data in this package contains only those movies not exceeding 400 minutes in length.
#' \url{http://had.co.nz/data/movies/}.
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 58771 rows and 23 variables
#' @name movieratings
#' 
#' @param title Title of the movie
#' @param year Year the movie was released
#' @param length Length of movie (in minutes)
#' @param rating Average IMDB user rating
#' @param votes Number of IMDB users who rated the movie
#' @param r1-10 Distribution of votes for each rating, to mid point of nearest decile: 0 = no votes, 4.5 = 1-9 percent votes, 14.5 = 11-19 percent of votes, etc. Due to rounding errors these may not sum to 100.
#' @param mpaa MPAA rating
#' @param genre Binary variables indicating whether movie belongs to any of the following genres: action, animation, comedy, drama, documentary, romance, short
#' 
NULL

#' The Diet Problem: "Nutritive Values of Common Foods per Dollar of Expenditure, August 15, 1944", from George Stigler's 1945 paper on "The Cost of Subsistence"
#'
#' A dataset of 10 rows and 14 columns describing 1944 prices and nutritional data for 14 food commodities.  Obtained from George Stigler's 1945 paper on "The Cost of Subsistence".
#'
#' @source George J. Stigler, "The Cost of Subsistence", Journal of Farm Economics, Vol. 27, No. 2 (May, 1945), pp. 303-314.
#' \url{http://www.jstor.org/stable/1231810}
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 10 rows and 14 columns
#' @name stigler
#' 
#' @param Calories Calories (in kilocalories) per dollar
#' @param Protein  Protein (in grams) per dollar
#' @param Calcium Calcium (in grams) per dollar
#' @param Iron Iron (in milligrams) per dollar
#' @param Vitamin.A Vitamin A (in 1000 Interntional Units) per dollar
#' @param Thiamine Thiamine (in milligrams) per dollar
#' @param Riboflavin Riboflavin (in milligrams) per dollar
#' @param Niacin Niacin (in milligrams) per dollar
#' @param Ascorbic.Acid Ascorbin Acid (in milligrams) per dollar
#' 
NULL

#' The Diet Problem: "Daily Allowances of Nutrients for a Moderately Active Man (weighing 154 pounds)" from George Stigler's 1945 paper on "The Cost of Subsistence"
#'
#' A vector describing 1943 dietary requirements for a "moderately active" man of 154 pounds.  Obtained from Table 1 of George Stigler's 1945 paper on "The Cost of Subsistence".
#'
#' @source George J. Stigler, "The Cost of Subsistence", Journal of Farm Economics, Vol. 27, No. 2 (May, 1945), pp. 303-314.
#' \url{http://www.jstor.org/stable/1231810}
#'
#' @docType data
#' @keywords datasets
#' @format A vector of length 9
#' @name nutrition
#' 
#' @param Calories Calories (in kilocalories)
#' @param Protein  Protein (in grams)
#' @param Calcium Calcium (in grams)
#' @param Iron Iron (in milligrams)
#' @param Vitamin.A Vitamin A (in 1000 Interntional Units)
#' @param Thiamine Thiamine (in milligrams)
#' @param Riboflavin Riboflavin (in milligrams)
#' @param Niacin Niacin (in milligrams)
#' @param Ascorbic.Acid Ascorbin Acid (in milligrams)
#'
NULL

#' Engel's Law - Engel Food Expenditures Data from the quantreg package for R
#'
#' Data on income and food expenditure for 235 working class households in 1857 Belgium.
#'
#' @source This dataset was used in Koenker and Bassett (1982) and obtained from the quantreg package for R.  Citations:  Koenker, R. and Bassett, G (1982) Robust Tests of Heteroscedasticity based on Regression Quantiles; Econometrica 50, 43-61.
#' \url{http://CRAN.R-project.org/package=quantreg}
#'
#' @docType data
#' @keywords datasets
#' @format A dataset containing 235 observations on 2 variables
#' @name engel
#'
#' @param income Annual household income (Belgian francs) 
#' @param foodexp Annual household food expenditure (Belgian francs)
#'
NULL
