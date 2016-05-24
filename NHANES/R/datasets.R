#' NHANES 2009-2012 with adjusted weighting
#' 
#' This is survey data collected by the US National Center for Health Statistics
#' (NCHS) which has conducted a series of health and nutrition surveys since the
#' early 1960's. Since 1999 approximately 5,000 individuals of all ages are
#' interviewed in their homes every year and complete the health examination
#' component of the survey. The health examination is conducted in a mobile
#' examination centre (MEC).
#' 
#'
#' @docType data
#' @name NHANES
#' @aliases NHANESraw
#' @usage data(NHANES)
#' @usage data(NHANESraw)
#' @examples
#' # Due to the sampling design, some races were over/under-sampled.
#' rbind(
#'   NHANES = table(NHANES$Race1) / nrow(NHANES),
#'   NHANESraw = table(NHANESraw$Race1) / nrow(NHANESraw),
#'   diff = (table(NHANES$Race1) - table(NHANESraw$Race1)) / nrow(NHANESraw)
#' )
#' # SmokeNow is only asked of people who answer Yes to Smoke100
#' if (require(mosaic)) {
#'   nhanes <- 
#'     NHANES %>% 
#'     mutate(
#'       SmokingStatus = derivedFactor(
#'         Current = SmokeNow == "Yes",
#'         Former = SmokeNow == "No",
#'         Never  = Smoke100 == "No"
#'         )
#'       )
#'    tally( ~SmokingStatus, data = nhanes )
#' }
#'
#' @format data frames with raw and resampled versions of the NHANES data.  See below for details 
#' and descriptions of the varaibles.
#'
#'
#' @section NHANES warning:
#' The following warning comes directly from the NHANES web site:
#' 
#' For NHANES datasets, the use of sampling weights and sample design variables 
#' is recommended for all analyses because the sample design is a clustered design 
#' and incorporates differential probabilities of selection.  
#' If you fail to account for the sampling parameters, you may obtain biased estimates 
#' and overstate significance levels.
#' 
#' @section Disclamer:
#' Please note that the data sets provided in this package are derived from the 
#' NHANES database and have been adapted for educational purposes.  
#' As such, they are NOT
#' suitable for use as a research database.  For research purposes you should
#' download original data files from the NHANES website and follow the analysis
#' instructions given there. Further details and relevant documentation can be
#' found on the following NHANES websites 
#' \itemize{
#' \item
#' \url{http://www.cdc.gov/nchs/nhanes.htm},
#' \item
#' \url{http://wwwn.cdc.gov/nchs/nhanes/search/nhanes11_12.aspx}, and 
#' \item
#' \url{http://wwwn.cdc.gov/nchs/nhanes/search/nhanes09_10.aspx}.
#' }
#' 
#' @source
#' These data were originally assembled by Michelle Dalrymple of Cashmere High School 
#' and Chris Wild of the University of Auckland, New Zealand for use in teaching 
#' statistics.
#' 
#'
#' @details
#'
#' The NHANES target population is "the non-institutionalized civilian resident
#' population of the United States".  
#' NHANES, (American National Health
#' and Nutrition Examination surveys), use complex survey designs (see
#' http://www.cdc.gov/nchs/data/series/sr_02/sr02_162.pdf) that oversample certain
#' subpopulations like racial minorities. Naive analysis of the original NHANES
#' data can lead to mistaken conclusions. The percentages of people from each
#' racial group in the data, for example, are quite different from the way they
#' are in the population.  
#' 
#' \code{NHANES} and \code{NHANESraw} 
#' each include 75 variables available for the 2009-2010 and 2011-2012 sample years.
#' \code{NHANESraw} has 20,293 observations of these variables plus four additional
#' variables that describe that sample weighting scheme employed. 
#' \code{NHANES} contains 10,000 rows of data resampled from 
#' \code{NHANESraw} to undo these oversampling effects. 
#' \code{NHANES} can be treated, for educational purposes, 
#' as if it were a simple random sample from the American population.  
#' 
#' A list of the variables in
#' the data set follows appears below along with variable descriptions and links to the original
#' NHANES documentation.
#'
#'

#' @section Study Variables:
#'   \describe{
#'     \item{SurveyYr}{Which survey the participant participated in.}
#'     \item{ID}{Participant identifier.}
#'   }
#'   
#' @section Demographic Variables:
#' For more information on these demographic variables, see
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2009-2010/DEMO_F.htm} or
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2011-2012/DEMO_G.htm}.
#' \describe{
#' 	\item{Gender}{Gender (sex) of study participant	coded as \code{male} or \code{female}}
#' 	\item{Age}{Age in years at screening of study participant.  Note:  Subjects 80 years or older were
#' 	  recorded as 80.}
#' 	\item{AgeDecade}{Categorical variable derived from age with levels \code{0-9}, \code{10-19}, \dots \code{70+}}	
#' 	\item{AgeMonths}{Age in months at screening of study participant.  Reported
#' 	  for participants aged 0 to 79 years for 2009 to 2010 data Reported for
#' 	  participants aged 0 to 2 years for 2011 to 2012 data.} 
#' 	\item{Race1}{Reported race of study participant: Mexican, Hispanic, White, Black, or Other.}
#' 	\item{Race3}{Reported race of study participant, including non-Hispanic Asian category:
#' 		Mexican, Hispanic, White, Black, Asian, or Other.  Not availale for 2009-10.}
#' 	\item{Education}{Educational level of study participant Reported for participants aged 20 years or older.
#' 		One of \code{8thGrade}, \code{9-11thGrade}, \code{HighSchool}, \code{SomeCollege}, or \code{CollegeGrad}.}
#'     \item{MaritalStatus}{Marital status of study participant.  Reported for participants aged 20 years or older.
#' 	  One of \code{Married}, \code{Widowed}, \code{Divorced}, \code{Separated}, \code{NeverMarried}, or \code{LivePartner} (living with partner).}
#' 
#' 	\item{HHIncome}{Total annual gross income for the household in US dollars.  One of 
#'       \code{0 - 4999}, \code{5000 - 9,999}, 
#' 	  \code{10000 - 14999}, \code{15000 - 19999}, \code{20000 - 24,999},
#' 	  \code{25000 - 34999}, \code{35000 - 44999}, \code{45000 - 54999}, \code{55000 - 64999}, \code{65000 - 74999},
#' 	  \code{75000 - 99999}, or \code{100000 or More}.}
#' 	\item{HHIncomeMid}{Numerical version of \code{HHIncome} derived from the middle income in each category}
#' 	\item{Poverty}{A ratio of family income to poverty guidelines.  Smaller numbers indicate more poverty}
#' 	\item{HomeRooms}{How many rooms are in home of study participant (counting kitchen but not bathroom).
#' 	  13 rooms = 13 or more rooms.}
#' 	 \item{HomeOwn}{One of \code{Home}, \code{Rent}, or \code{Other} indicating whether 
#' 	 the home of study participant or someone in their family is owned, rented or occupied 
#' 	 by some other arrangement.}	
#' }
#' 
#' @section Physical Measurements:
#' For more information on body measurements, see 
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2009-2010/BMX_F.htm} and 
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2011-2012/BMX_G.htm}.
#' \describe{
#' \item{Weight}{Weight in kg}
#' \item{Length}{Recumbent length in cm. Reported for participants aged 0 - 3 years.
#' }
#' \item{HeadCirc}{Head circumference in cm.  
#' Reported for participants aged 0 years (0 - 6 months).}
#' \item{Height}{Standing height in cm. Reported for participants aged 2 years or older.}
#' \item{BMI}{Body mass index (weight/height2 in kg/m2). 
#' Reported for participants aged 2 years or older.}
#' \item{BMICatUnder20yrs}{Body mass index category. 
#'       Reported for participants aged 2 to 19 years.
#'       One of 
#'       \code{UnderWeight} (BMI < 5th percentile)
#'       \code{NormWeight} (BMI 5th to < 85th percentile),
#'       \code{OverWeight} (BMI 85th to < 95th percentile),
#'       \code{Obese} (BMI >= 95th percentile).}
#' \item{BMI_WHO}{Body mass index category. 
#' Reported for participants aged 2 years or older.
#' One of \code{12.0_18.4}, \code{18.5_24.9}, \code{25.0_29.9}, or \code{30.0_plus}.}
#' \item{Pulse}{60 second pulse rate}
#' \item{BPSysAve}{Combined systolic blood pressure reading, 
#' following the procedure outlined for BPXSAR.}
#' \item{BPDiaAve}{Combined diastolic blood pressure reading, 
#' following the procedure outlined for BPXDAR.}
#' \item{BPSys1}{Systolic blood pressure in mm Hg -- first reading}
#' \item{BPDia1}{Diastolic blood pressure in mm Hg -- second reading (consecutive readings)}
#' \item{BPSys2}{Systolic blood pressure in mm Hg -- second reading (consecutive readings)}
#' \item{BPDia2}{Diastolic blood pressure in mm Hg -- second reading}
#' \item{BPSys3}{Systolic blood pressure in mm Hg  third reading (consecutive readings)}
#' \item{BPDia3}{Diastolic blood pressure in mm Hg -- third reading (consecutive readings)}
#' \item{Testosterone}{Testerone total (ng/dL). Reported for participants aged 6 years or 
#'       older. Not available for 2009-2010.}
#' }
#' @section Health Variables:
#' For more information on these variables, see
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2009-2010/HDL_F.htm} or
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2011-2012/HDL_G.htm}.
#' \describe{
#' \item{DirectChol}{Direct HDL cholesterol in mmol/L.
#'       Reported for participants aged 6 years or older.}
#' \item{TotChol}{Total HDL cholesterol in mmol/L.
#'       Reported for participants aged 6 years or older.}
#' \item{UrineVol1}{Urine volume in mL -- first test.
#'       Reported for participants aged 6 years or older.}
#' \item{UrineFlow1}{Urine flow rate 
#'       (urine volume/time since last urination) 
#'       in mL/min -- first test. 
#'       Reported for participants aged 6 years or older.}
#' \item{UrineVol2}{Urine volume in mL -- second test.
#'       Reported for participants aged 6 years or older.}
#' \item{UrineFlow2}{Urine flow rate 
#'       (urine volume/time since last urination) 
#'       in mL/min -- second test. 
#'       Reported for participants aged 6 years or older.}
#' \item{Diabetes}{Study participant told by a doctor or health professional 
#' that they have diabetes. Reported for participants aged 1 year or older
#' as \code{Yes} or \code{No}.}
#' \item{DiabetesAge}{Age of study participant when first told they had diabetes.
#'       Reported for participants aged 1 year or older.}
#' \item{HealthGen}{Self-reported rating of participant's health in general
#'             Reported for participants aged 12 years or older.  
#'             One of  \code{Excellent}, \code{Vgood}, \code{Good}, \code{Fair}, or \code{Poor}.}
#' \item{DaysPhysHlthBad}{Self-reported number of days participant's physical health was 
#'       not good out of the past 30 days. Reported for participants aged 12 years or older.}
#' \item{DaysMentHlthBad}{Self-reported number of days participant's mental health was not 
#'       good out of the past 30 days. Reported for participants aged 12 years or older.}
#' \item{LittleInterest}{Self-reported number of days where participant had little 
#' interest in doing things. Reported for participants aged 18 years or older.
#' One of  \code{None}, \code{Several}, \code{Majority} (more than half the days), 
#' or \code{AlmostAll}.}
#' \item{Depressed}{Self-reported number of days where participant felt down, 
#'       depressed or hopeless. Reported for participants aged 18 years or older.
#' One of  \code{None}, \code{Several}, \code{Majority} (more than half the days), 
#' or \code{AlmostAll}.}
#' \item{nPregnancies}{How many times participant has been pregnant.
#'       Reported for female participants aged 20 years or older.}
#' \item{nBabies}{How many of participants deliveries resulted in live births.
#'                Reported for female participants aged 20 years or older.}
#' \item{PregnantNow}{Pregnancy status at the time of the health examination 
#' was ascertained for females 8-59 years of age. 
#' Due to disclosure risks pregnancy status was only be released for women 20-44 
#' years of age. The information used included urine pregnancy test results 
#' and self-reported pregnancy status. Urine pregnancy tests were performed prior 
#' to the dual energy x-ray absorptiometry (DXA) exam. 
#' Persons who reported they were pregnant at the time of exam were assumed to 
#' be pregnant. As a result, if the urine test was negative, but the subject 
#' reported they were pregnant, the status was coded as \code{"Yes"}.
#' If the urine pregnancy results were negative and the respondent stated that they 
#' were not pregnant, the respondent was coded as \code{"No"} If the urine pregnancy 
#' results were negative and the respondent did not know her pregnancy status, 
#' the respondent was coded \code{"unknown"}  Persons who were interviewed, 
#' but not examined also have a value of \code{"unknown"}. In addition
#' there are missing values. 
#' }
#' \item{Age1stBaby}{Age of participant at time of first live birth. 
#'                  14 years or under = 14,  45 years or older = 45. 
#'                  Reported for female participants aged 20 years or older.}
#' \item{SleepHrsNight}{Self-reported number of hours study participant usually gets 
#'                     at night on weekdays or workdays. Reported for participants aged 
#'                     16 years and older.}
#' \item{SleepTrouble}{Participant has told a doctor or other health professional that they 
#'                     had trouble sleeping. Reported for participants aged 16 years and older.
#'                     Coded as \code{Yes} or \code{No}.}
#' }
#' @section Lifestyle Variables:
#' More information about these variables is available at
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2009-2010/SMQ_F.htm} or 
#' \url{http://www.cdc.gov/nchs/nhanes/nhanes2011-2012/SMQ_G.htm}.
#' \describe{
#' \item{PhysActive}{Participant does moderate or vigorous-intensity sports, fitness or recreational 
#' activities (Yes or No).  Reported for participants 12 years or older.}
#' \item{PhysActiveDays}{Number of days in a typical week that participant does moderate or 
#' vigorous-intensity activity.  Reported for participants 12 years or older.}
#' \item{TVHrsDay}{Number of hours per day on average participant watched TV over the 
#' past 30 days. Reported for participants 2 years or older.  
#' One of  \code{0_to_1hr}, \code{1_hr}, \code{2_hr}, \code{3_hr}, \code{4_hr}, \code{More_4_hr}.
#' Not available 2009-2010.} 
#' \item{CompHrsDay}{Number of hours per day on average participant used a computer or gaming 
#' device over the past 30 days.  Reported for participants 2 years or older.  One of 
#' \code{0_hrs}, \code{0_to_1hr}, \code{1_hr}, \code{2_hr}, \code{3_hr}, \code{4_hr}, \code{More_4_hr}.
#' Not available 2009-2010.} 
#' \item{TVHrsDayChild}{Number of hours per day on average participant watched TV over the past 30 days.
#'  Reported for participants 2 to 11 years.
#' Not available 2011-2012.} 
#' \item{CompHrsDayChild}{Number of hours per day on average participant used a computer or gaming device 
#' over the past 30 days.  Reported for participants 2 to 11 years old.
#' Not available 2011-2012.}
#' \item{Alcohol12PlusYr}{Participant has consumed at least 12 drinks of any type of alcoholic beverage 
#' in any one year. Reported for participants 18 years or older as Yes or No.}
#' \item{AlcoholDay}{Average number of drinks consumed on days that participant drank alcoholic 
#' beverages.  Reported for participants aged 18 years or older.}
#' \item{AlcoholYear}{Estimated number of days over the past year that participant drank 
#' alcoholic beverages. Reported for participants aged 18 years or older.}
#' \item{SmokeNow}{Study participant currently smokes cigarettes regularly. 
#' Reported for participants aged 20 years or older as \code{Yes} or \code{No}, provieded they
#' answered Yes to having somked 100 or more cigarettes in their life time.  All subjects who 
#' have not smoked 100 or more cigarettes are listed as \code{NA} here.}
#' \item{Smoke100}{Study participant has smoked at least 100 cigarettes in their entire life.
#'  Reported for participants aged 20 years or older as \code{Yes} or \code{No}.}
#' \item{SmokeAge}{Age study participant first started to smoke cigarettes fairly regularly. 
#' Reported for participants aged 20 years or older.}
#' \item{Marijuana}{Participant has tried marijuana. Reported for participants aged 18 to 59 years as
#' \code{Yes} or \code{No}.}
#' \code{AgeFirstMarij}{Age participant first tried marijuana. Reported for participants aged 18 to 59 years.}
#' \item{RegularMarij}{Participant has been/is a regular marijuana user (used at least once a month for a year).
#'  Reported for participants aged 18 to 59 years as \code{Yes} or \code{No}.}
#' \item{AgeRegMarij}{Age of participant when first started regularly using marijuana. 
#' Reported for participants aged 18 to 59 years.}
#' \item{HardDrugs}{Participant has tried cocaine, crack cocaine, heroin or methamphetamine. 
#' Reported for participants aged 18 to 69 years as \code{Yes} or \code{No}.}
#' \item{SexEver}{Participant had had vaginal, anal, or oral sex.  
#' Reported for participants aged 18 to 69 years as \code{Yes} or \code{No}.}
#' \item{SexAge}{Age of participant when had sex for the first time. 
#' Reported for participants aged 18 to 69 years.}
#' \item{SexNumPartnLife}{Number of opposite sex partners participant has had any kind of sex with 
#' over their lifetime.  Reported for participants aged 18 to 69 years.}
#' \item{SexNumPartYear}{Number of opposite sex partners participant has had any kind of sex with over 
#' the past 12 months. Reported for participants aged 18 to 59 years.}
#' \item{SameSex}{Participant has had any kind of sex with a same sex partner. 
#' Reported for participants aged 18 to 69 years ad \code{Yes} or \code{No}.}
#' \item{SexOrientation}{participant's sexual orientation (self-described).    
#' Reported for participants aged 18 to 59 years.
#'  One of \code{Heterosexual}, \code{Homosexual}, \code{Bisexual}.}
#' }
#' 
#' @section Weighting Variables (\code{NHANESraw} only):
#' \describe{
#'   \item{WTINT2YR, WTMEC2YR,  SDMVPSU,  SDMVSTRA}{Sample weighting variables.  
#'   For more details see one of the following.
#'   \itemize{
#'   \item
#'   \url{http://www.cdc.gov/Nchs/tutorials/environmental/orientation/sample_design/index.htm}
#'   \item
#'   \url{http://www.cdc.gov/nchs/nhanes/nhanes2009-2010/DEMO_F.htm#WTINT2YR} and 
#'   \item
#'   \url{http://www.cdc.gov/nchs/nhanes/nhanes2011-2012/DEMO_G.htm#WTINT2YR} 
#'   }
#'   }
#' }
#' 

NA
