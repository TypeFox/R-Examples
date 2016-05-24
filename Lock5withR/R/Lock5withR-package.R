#' Lock5 Datasets
#' 
#' Datasets for Statistics: Unlocking the Power of Data \cr by Lock, Lock,
#' Lock, Lock and Lock
#' 
#' @name Lock5Data-package
#' @aliases Lock5Data-package Lock5Data
#' @docType package
#' @author Robin Lock \email{rlock@@stlawu.edu}
#' 
NA

#' American Community Survey
#' 
#' Data from a sample of individuals in the American Community Survey
#' 
#' The American Community Survey, administered by the US Census Bureau, is
#' given every year to a random sample of about 3.5 million households (about
#' 3\% of all US households).  Data on a random sample of 1\% of all US
#' residents are made public (after ensuring anonymity), and we have selected a
#' random sub-sample of n = 1000 from the 2010 data for this dataset.
#' 
#' @name ACS
#' @docType data
#' @format A dataset with 1000 observations on the following 9 variables.
#' \itemize{
#' \item{\code{Sex}} { a factor with levels \code{Female} and \code{Male} }
#' \item{\code{Age}} {age in years }
#' \item{\code{Married}} {a factor with levels \code{Married} and \code{Not Married} }
#' \item{\code{Income}}
#' {Wages and salary for the past 12 months (in \\$1,000's) }
#' \item{\code{HoursWk}}
#' {Hours of work per week }
#' \item{\code{Race}} { a factor with levels \code{asian}, \code{black},
#' \code{white}, or \code{other} } 
#' \item{\code{USCitizen}} {a factor wit two levels (\code{Citizen} and 
#' \code{Noncitizen}).}
#' \item{\code{HealthInsurance}} {a factor with two levels (\code{Insured} and 
#' \code{Uninsured}}
#' \item{\code{Language}} { a factor with two levels (\code{English} and \code{Other})}
#' }
#' @source The full public dataset can be downloaded at 
#' \url{http://www.census.gov/programs-surveys/acs/data/pums.html}, and the
#' full list of variables are at 
#' \url{http://www.census.gov/programs-surveys/acs/guidance.html}.
# @source The full public dataset can be downloaded at 
# \url{http://www.census.gov/acs/www/data_documentation/pums_data/}
# and the full list of variables is at 
# \url{http://www.census.gov/acs/www/Downloads/data/documentation/pums/DataDict/PUMSDataDict10.pdf}.
#' @keywords datasets
NA





#' AllCountries
#' 
#' Data on the countries of the world
#' 
#' Most data from 2008 to avoid many missing values in more recent years.
#' 
#' @name AllCountries
#' @docType data
#' @format A dataset with 213 observations on the following 18 variables.
#' \itemize{
#' \item{\code{Country}} {Name of the country}
#' \item{\code{Code}} {Three letter country code}
#' \item{\code{LandArea}} {Size in sq. kilometers}
#' \item{\code{Population}} {Population in millions}
#' \item{\code{Energy}} {Energy usage (kilotons of oil)}
#' \item{\code{Rural}} {Percentage of population living in rural areas}
#' \item{\code{Military}} {Percentage of government expenditures
#' directed toward the military}
#' \item{\code{Health}} {Percentage of government
#' expenditures directed towards healthcare}
#' \item{\code{HIV}} {Percentage of
#' the population with HIV}
#' \item{\code{Internet}} {Percentage of the population
#' with access to the internet}
#' \item{\code{Developed}} {A numeric code for the level of development
#' based on kilowatt hours per capita}
#' \item{\code{kwhPerCap}} {An ordered factor of categories for kilowatt
#' hours per capita, \code{under 2500}, \code{2500 to 5000}, or \code{over 5000}}
#' \item{\code{BirthRate}} {Births per 1000 people}
#' \item{\code{ElderlyPop}} {Percentage of the population at least 65 year old}
#' \item{\code{LifeExpectancy}} {Average life expectancy (years)}
#' \item{\code{CO2}} {CO2 emissions
#' (metric tons per capita)}
#' \item{\code{GDP}} {Gross Domestic Product (per capita)}
#' \item{\code{Cell}} {Cell phone subscriptions (per 100 people)}
#' \item{\code{Electricity}} {Electric power consumption (kWh per capita)}
#' }
#' @source Data collected from the World Bank website, worldbank.org.
#' @keywords datasets
NA





#' AP Multiple Choice
#' 
#' Correct responses on Advanced Placement multiple choice exams
#' 
#' Correct responses from multiple choice sections for a sample of released
#' Advanced Placement exams
#' 
#' @name APMultipleChoice
#' @docType data
#' @format A dataset with 400 observations on the following variable.
#' \itemize{
#\ \item{\code{Answer}} {Correct response: \code{A}, \code{B},
#' \code{C}, \code{D}, or \code{E} }
#' @source Sample exams from several disciplines at
#' \url{http://apcentral.collegeboard.com}
#' @keywords datasets
NA





#' April 14th Temperatures
#' 
#' Temperatures in Des Moines, IA and San Francisco, CA on April 14th
#' 
#' Average temperature for the day of April 14th in each of 16 years from
#' 1995-2010
#' 
#' @name April14Temps
#' @docType data
#' @format A data frame with 16 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Year}} {1995-2010}
#' 
#'    \item{\code{DesMoines}} {Temperature in Des Moines (degrees F)}
#' 
#'    \item{\code{SanFrancisco}} {Temperature in San Francisco (degrees F)} 
#'    }
#' @source The University of Dayton Average Daily Temperature Archive at
#' \url{http://academic.udayton.edu/kissock/http/Weather/citylistUS.htm}
#' @keywords datasets
#' @examples
#' 
#' data(April14Temps)
#' 
NA





#' Baseball Hits
#' 
#' Number of hits and wins for Major League Baseball teams
#' 
#' Data from the 2010 Major League Baseball regular season.
#' 
#' @name BaseballHits
#' @docType data
#' @format A dataset with 30 observations on the following 14 variables.
#' \itemize{
#'   \item{\code{Team}} {Name of baseball team}
#'   \item{\code{League}} {Either \code{AL} or \code{NL}}
#'   \item{\code{Wins}} {Number of wins for the season}
#'   \item{\code{Runs}} {Number of runs scored}
#'   \item{\code{Hits}} {Number of hits}
#'   \item{\code{Doubles}} {Number of doubles}
#'   \item{\code{Triples}} {Number of triples}
#'   \item{\code{HomeRuns}} {Number of home runs}
#'   \item{\code{RBI}} {Number of runs batted in}
#'   \item{\code{StolenBases}} {Number of stolen bases}
#'   \item{\code{CaughtStealing}} {Number of times caught stealing}
#'   \item{\code{Walks}} {Number of walks}
#'   \item{\code{Strikeouts}} {Number of stikeouts}
#'   \item{\code{BattingAvg}} {Team batting average}
#'   }
#' @source
#' \url{http://www.baseball-reference.com/leagues/MLB/2011-standard-batting.shtml}
#' @keywords datasets
NA





#' Baseball Game Times
#' 
#' Information for a sample of 30 Major League Baseball games played during the
#' 2011 season.
#' 
#' Data from a sample of boxscores for Major League Baseball games played in
#' August 2011.
#' 
#' @name BaseballTimes
#' @docType data
#' @format A dataset with 30 observations on the following 9 variables.
#' \itemize{
#'   \item{\code{Away} } { Away team name}
#'   \item{\code{Home}} { Home team name}
#'   \item{\code{Runs}} { Total runs scored (both teams)}
#'   \item{\code{Margin}} { Margin of victory}
#'   \item{\code{Hits}} { Total number of hits (both teams)}
#'   \item{\code{Errors}} { Total number of errors (both teams)}
#' \item{\code{Pitchers}} { Total number of pitchers used (both teams)}
#' \item{\code{Walks}} { Total number of walks (both teams)}
#'   \item{\code{Time}} { Elapsed time for game (in minutes) }
#'   }
#' @source \url{http://www.baseball-reference.com/boxes/2011.shtml}
#' @keywords datasets
NA





#' Benford data
#' 
#' Two examples to test Benford's Law
#' 
#' Leading digits from 1188 addresses sampled from a phone book and 7273
#' amounts from invoices sampled at a company.
#' 
#' @name Benford
#' @docType data
#' @format A dataset with 9 observations on the following 4 variables.
#' \itemize{
#' \item{\code{Digit}} {Leading digit (1-9)}
#' \item{\code{BenfordP}} {Expected proportion according to Benford's law}
#' \item{\code{Address}} {Frequency as a first digit in an address}
#' \item{\code{Invoices}} {Frequency as the first digit in invoice amounts}
#' }
#' @source Thanks to Prof. Richard Cleary for providing the data
#' @keywords datasets
NA





#' Bike Commute
#' 
#' Commute times for two kinds of bicycle
#' 
#' Data from a personal experiment to compare commuting time based on a
#' randomized selection between two bicycles made of different materials.
#' 
#' @name BikeCommute
#' @docType data
#' @format A dataset with 56 observations on the following 9 variables.
#' \itemize{
#' \item{\code{Bike}} {Type of material (\code{Carbon} or \code{Steel})}
#' \item{\code{Date}} {Date of the bike commute as a string}
#' \item{\code{DMY}}  {Date of the bike commute as a date object}
#' \item{\code{Distance}} {Length of commute (in miles)}
#' \item{\code{Time}} {Total commute time (hours:minutes:seconds)}
#' \item{\code{Minutes}} {Time converted to minutes}
#' \item{\code{AvgSpeed}} {Average speed during the ride (miles per hour)}
#' \item{\code{TopSpeed}} {Maximum speed (miles per hour)}
#' \item{\code{Seconds}} {Time converted to seconds}
#' \item{\code{MonthStr}} {Month} 
#' \code{MonthNum} {Numeric value of month}
#' \item{\code{Month}} {Month coded as one of \code{1Jan},
#' \code{2Feb}, \code{3Mar}, \code{4Apr}, \code{5May}, 
#' \code{6June}, or \code{7July}},
#' }
#' @references "Bicycle weight and commuting time: randomised trial," 
#' in British Medical Journal, BMJ 2010;341:c6801.
#' @source Thanks to Dr. Groves for providing his data.
#' @keywords datasets
NA





#' Body Measurements
#' 
#' 
#' @name BodyFat
#' @docType data
#' @format A data frame with 100 observations on the following 10 variables.
#' \itemize{ 
#'    \item{\code{Bodyfat}} {a numeric vector} 
#'    \item{\code{Age}} {a numeric vector} 
#'    \item{\code{Weight}} {a numeric vector}
#'    \item{\code{Height}} {a numeric vector} 
#'    \item{\code{Neck}} {a numeric vector} 
#'    \item{\code{Chest}} {a numeric vector} 
#'    \item{\code{Abdomen}} {a numeric vector} 
#'    \item{\code{Ankle}} {a numeric vector}
#'    \item{\code{Biceps}} {a numeric vector} 
#'    \item{\code{Wrist}} {a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(BodyFat)
NA




#' Body Temperatures
#' 
#' Sample of 50 body temperatures
#' 
#' Body temperatures and pulse rates for a sample of 50 healthy adults.
#' 
#' @name BodyTemp50
#' @docType data
#' @format A dataset with 50 observations on the following 3 variables.
#' \itemize{
#' \item{\code{BodyTemp}} {Body temperatures in degrees F}
#' \item{\code{Pulse}} {Pulse rates (beats per minute)}
#' \item{\code{Sex}} {\code{Female} or \code{Male}}
#' }
#' @source Shoemaker, "What's Normal: Temperature, Gender and Heartrate", 
#' Journal of Statistics Education, Vol. 4, No. 2 (1996)
#' @keywords datasets
NA


#' Bootstrap Correlations for Atlanta Commutes
#' 
#' Boostrap correlations between Time and Distance for 500 commuters in Atlanta
#' 
#' Correlations for bootstrap samples of Time vs. Distance for the data on
#' Atlanta commuters in CommuteAtlanta.
#' 
#' @name BootAtlantaCorr
#' @docType data
#' @format A data frame with 1000 observations on the following variable.
#' \itemize{ 
#'    \item{\code{CorrTimeDist}} {Correlation between Time and Distance
#' for a bootstrap sample of Atlanta commuters} }
#' @source Computer simulation
#' @keywords datasets
#' @examples
#' 
#' data(BootAtlantaCorr)
#' 
NA



#' Caffeine Taps
#' 
#' Finger tap rates with and without caffeine
#' 
#' Results from a double-blind experiment where a sample of male college
#' students to tap their fingers at a rapid rate. The sample was then divided
#' at random into two groups of ten students each.  Each student drank the
#' equivalent of about two cups of coffee, which included about 200 mg of
#' caffeine for the students in one group but was decaffeinated coffee for the
#' second group. After a two hour period, each student was tested to measure
#' finger tapping rate (taps per minute). The goal of the experiment was to
#' determine whether caffeine produces an increase in the average tap rate.
#' 
#' @name CaffeineTaps
#' @docType data
#' @format A dataset with 20 observations on the following 2 variables.
#' \itemize{
#' \item{\code{Taps}} {Number of finger taps in one minute}
#' \item{\code{Group}} {Treatment with levels \code{Caffeine} 
#' and \code{No Caffeine}}
#' \item{\code{Caffeine}} {a recoding of \code{Group} with levels \code{"Yes"} and \code{"No"}}
#' }
#' @source Hand, Daly, Lund, McConway and Ostrowski, Handbook of Small Data
#' Sets, Chapman and Hall, London (1994), pp. 40
#' @keywords datasets
NA





#' CAOS Exam Scores
#' 
#' Scores on a pre-test and post-test of basic statistics concepts
#' 
#' The CAOS (Comprehensive Assessment of Outcomes in First Statistics Course)
#' exam is designed to measure comprehension basic statistical ideas in an
#' introductory statistics course.  This dataset has scores for ten students
#' who took the CAOS pre-test at the start of a course and the post-test during
#' the course itself.  Each exam consists of 40 multiple choice questions and
#' the score is the percentage correct.
#' 
#' @name CAOSExam
#' @docType data
#' @format A datset with 10 observations on the following 3 variables.
#' \itemize{
#' \item{\code{Student}} {ID code for student}
#' \item{\code{Pretest}} {CAOS Pretest score}
#' \item{\code{Posttest}} {CAOS Posttest score}
#' }
#' 
#' @references Find out more about the CAOS exam at
#' \url{http://app.gen.umn.edu/artist/caos.html}
#' @source A sample of 10 students from an introductory statisics course.
#' @keywords datasets
NA





#' Atmospheric carbon dioxide levels by year
#' 
#' Carbon dioxide levels in the atmosphere over a 50 year span from 1960-2010.
#' 
#' @name AtmosphericCO2
#' @docType data
#' @format A dataset with 11 observations on the following 2 variables.
#' \itemize{ 
#' \item{\code{Year}} {Every five years from 1960-2010}
#' \item{\code{CO2}} {Carbon dioxide level in parts per million}
#' }
#' @source Dr. Pieter Tans, NOAA/ESRL (www.esrl.noaa.gov/gmd/ccgg/trends/).
#' Values recorded at the Mauna Loa Observatory in Hawaii.
#' @keywords datasets
NA





#' Breakfast Cereals
#' 
#' Nutrition information for a sample of 30 breakfast cereals
#' 
#' Nutrition contents for a sample of breakfast cererals, derived from
#' nutrition lables.  Values are per cup of cereal (rather than per serving).
#' 
#' @name Cereal
#' @docType data
#' @format A dataset with 30 observations on the following 10 variables.
#' \itemize{
#' \item{\code{Name}} {Brand name of cereal}
#'    \item{ \code{Company}} {Manufacturer coded as 
#'			\code{G}=General Mills, \code{K}=Kellog's or \code{Q}=Quaker}
#'    \item{\code{Serving}} {Serving size (in cups)}
#'    \item{\code{Calories}} {Calories (per cup)}
#'    \item{ \code{Fat}} {Fat (grams per cup)}
#'    \item{ \code{Sodium}} {Sodium (mg per cup)}
#'    \item{ \code{Carbs}} {Carbohydrates (grams per cup)}
#'    \item{ \code{Fiber}} {Dietary Fiber (grams per  cup)}
#'    \item{ \code{Sugars}} {Sugars (grams per cup)}
#'    \item{ \code{Protein}} {Protein (grams per cup)}
#' }
#' @source Cereal data obtained from nutrition labels at
#' \url{http://www.nutritionresource.com/foodcomp2.cfm?id=0800}
#' @keywords datasets
NA



#' Coacaine Treatment
#' 
#' Relapse/no relapse responses to three different treatments for cocaine
#' addiction
#' 
#' Data from an experiment to investigate the effectiveness of the two drugs,
#' desipramine and lithium, in the treatment of cocaine addiction. Subjects
#' (cocaine addicts seeking treatment) were randomly assigned to take one of
#' the treatment drugs or a placebo.  The response variable is whether or not
#' the subject relapsed (went back to using cocaine) after the treatment.
#' 
#' @name CocaineTreatment
#' @docType data
#' @format A dataset with 72 observations on the following 2 variables.
#' \itemize{
#'  \item{\code{Drug}} {Treatment drug: \code{Desipramine},
#' 			\code{Lithium}, or \code{Placebo}}
#' 	\item{\code{Relapse}} { Did the patient relapse? \code{no} or \code{yes}}
#' }
#' @source Gawin, F., et.al., "Desipramine Facilitation of Initial Cocaine
#' Abstinence", Archives of General Psychiatry, 1989; 46(2): 117 - 121.
#' @keywords datasets
NA





#' Cola Calcium
#' 
#' Calcium excretion with diet cola and water
#' 
#' A sample of 16 healthy women aged 18 - 40 were randomly assigned to drink 24
#' ounces of either diet cola or water. Their urine was collected for three
#' hours after ingestion of the beverage and calcium excretion (in mg.) was
#' measured . The researchers were investigating whether diet cola leaches
#' calcium out of the system, which would increase the amount of calcium in the
#' urine for diet cola drinkers.
#' 
#' @name ColaCalcium
#' @docType data
#' @format A dataset with 16 observations on the following 2 variables.
#' \itemize{
#'  	\item{\code{Drink}} {Type of drink: 
#'			\code{Diet cola} or \code{Water}}
#' 		\item{\code{Calcium}} {Amount of calcium excreted (in mg.)}
#' }
#' @source Larson, Amin, Olsen, and Poth, Effect of Diet Cola on Urine Calcium
#' Excretion, Endocrine Reviews, 31[3]: S1070, June 2010. These data are
#' recreated from the published summary statistics, and are estimates of the
#' actual data.
#' @keywords datasets
NA





#' Commute Atlanta
#' 
#' Commute times and distance for a sample of 500 people in Atlanta
#' 
#' Data from the US Census Bureau's American Housing Survey (AHS) which
#' contains information about housing and living conditions for samples from
#' certain metropolitan areas.  These data were extracted respondents in the
#' Atlanta metropolitan area. They include only cases where the respondent
#' worked somewhere other than home.  Values show the time (in minutes) and
#' distance (in miles) that respondents typically traveled on their commute to
#' work each day as well as age and sex.
#' 
#' @name CommuteAtlanta
#' @docType data
#' @format A data frame with 500 observations on the following 5 variables.
#' \itemize{ St. Louis }
#' @source Sample chosen using DataFerret at
#' \url{http://www.thedataweb.org/index.html}.
#' @keywords datasets
#' @examples
#' 
#' data(CommuteAtlanta)
#' 
NA





#' Comute times in St. Louis
#' 
#' Commute times and distance for a sample of 500 people in St. Louis
#' 
#' Data from the US Census Bureau's American Housing Survey (AHS) which
#' contains information about housing and living conditions for samples from
#' certain metropolitan areas.  These data were extracted respondents in the
#' St. Louis metropolitan area. They include only cases where the respondent
#' worked somewhere other than home.  Values show the time (in minutes) and
#' distance (in miles) that respondents typically traveled on their commute to
#' work each day as well as age and sex.
#' 
#' @name CommuteStLouis
#' @docType data
#' @format A dataset with 500 observations on the following 5 variables.
#' \itemize{
#'   	\item{\code{City}} {\code{St. Louis}}
#'    	\item{\code{Age}} {Age of the respondent (in years)}
#' 		\item{\code{Distance}} {Commute distance (in miles)}
#'		\item{\code{Time}} {Commute time (in minutes)}
#'		\item{\code{Sex}} {\code{F} or \code{M}}
#' }
#' @source Sample chosen using DataFerret at
#' \url{http://www.thedataweb.org/index.html}.
#' @keywords datasets
NA





#' Compassionate Rats
#' 
#' Would a rat attempt to free a trapped rat?
#' 
#' In a recent study, some rats showed compassion by freeing another trapped
#' rat, even when chocolate served as a distraction and even when the rats
#' would then have to share the chocolate with their freed companion.
#' 
#' @name CompassionateRats
#' @docType data
#' @format A dataset with 30 observations on the following 2 variables.
#' \itemize{
#'		\item{\code{Sex}} {Sex of the rat: coded as \code{F} or
#' \code{M}}
#'		\item{\code{Empathy}} {Freed the trapped rat?  \code{no} or \code{yes}}
#' }
#' @source Bartal I.B., Decety J., and Mason P., "Empathy and Pro-Social
#' Behavior in Rats," Science, 2011; 224(6061):1427-1430.
#' @keywords datasets
NA





#' Cricket Chirps
#' 
#' Cricket chirp rate and temperature
#' 
#' The data were collected by E.A. Bessey and C.A. Bessey who measured chirp
#' rates for crickets and temperatures during the summer of 1898.
#' 
#' @name CricketChirps
#' @docType data
#' @format A dataset with 7 observations on the following 2 variables.
#' \itemize{
#'		\item{\code{Temperature}} {Air temperature in degrees F}
#' 		\item{\code{Chirps}} {Cricket chirp rate (chirps per minute)}
#' }
#' @source From E.A Bessey and C.A Bessey, Further Notes on Thermometer
#' Crickets, American Naturalist, (1898) 32, 263-264.
#' @keywords datasets
NA





#' Digit counts
#' 
#' Digits from social security numbers and student selected random numbers
#' 
#' A sample of students were asked to pick a random four digit number.  The
#' numbers are given in the dataset, along with separate columns for each of
#' the four digits.  The data also show the last two digits of each student's
#' social security number (SSN).
#' 
#' @name Digits
#' @docType data
#' @format A data frame with 150 observations on the following 7 variables.
#' \itemize{ 
#'    \item{\code{Random}} {Four digit random numbers picked by a
#' sample of students} 
#'    \item{\code{RND1}} {First digit}
#' 
#'    \item{\code{RND2}} {Second digit} 
#'    \item{\code{RND3}} {Third digit}
#' 
#'    \item{\code{RND4}} {Fourth digit} 
#'    \item{\code{SSN8}} {Eighth digit of social
#' security number} 
#'    \item{\code{SSN9}} {Last digit of social security number} }
#' @source In class-student surveys from several classes.
#' @keywords datasets
#' @examples
#' 
#' data(Digits)
#' 
NA





#' Dog/Owner matches
#' 
#' Experiment to match dogs with owners
#' 
#' Pictures were taken of 25 owners and their purebred dogs, selected at random
#' from dog parks. Study participants were shown a picture of an owner together
#' with pictures of two dogs (the owner's dog and another random dog from the
#' study) and asked to choose which dog most resembled the owner. Each
#' dog-owner pair was viewed by 28 naive undergraduate judges, and the pairing
#' was deemed "correct" (yes) if the majority of judges (more than 14) chose
#' the correct dog to go with the owner.
#' 
#' @name DogOwner
#' @docType data
#' @format A data frame with 25 observations on the following variable.
#' \itemize{ 
#'    \item{\code{Match}} {Was the dog correctly paired with it's
#' owner? \code{no} or \code{yes}} }
#' @source Roy and Christenfeld, Do Dogs Resemble their Owners?, Psychological
#' Science, Vol. 15, No. 5, 2004, pp. 361 - 363.
#' @keywords datasets
#' @examples
#' 
#' data(DogOwner)
#' 
NA





#' Election Margin
#' 
#' Approval rating and election margin for recent presidential elections
#' 
#' Data include US Presidential elections since 1940 in which an incumbent was
#' running for president.  The approval rating for the sitting president is
#' comapred to the margin of victory/defeat in the election.
#' 
#' @name ElectionMargin
#' @docType data
#' @format A data frame with 11 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{Year}} {Certain election years from 1940-2004}
#' 
#'    \item{\code{Candidate}} {Incumbent US president}
#' 
#'    \item{\code{Approval}} {Presidential approval rating at time of election}
#' 
#'    \item{\code{Margin}} {margin of victory/defeat (as a percentage)}
#' 
#'    \item{\code{Result}} {Outcome of the election for the imbumbent: \code{Lost}
#' or \code{Won}} }
#' @source Silver, Nate, "Approval Ratings and Re-Election Odds",
#' fivethirtyeight.com, posted January 28, 2011.
#' @keywords datasets
#' @examples
#' 
#' data(ElectionMargin)
#' 
NA





#' Employed in American Community Survey
#' 
#' Employed individiuals from the American Community Survey (ACS) dataset
#' 
#' This is a subset of the ACS dataset including only 431 individuals who were
#' employed.
#' 
#' @name EmployedACS
#' @docType data
#' @format A dataset with 431 observations on the following 9 variables.
#' \itemize{
#'		\item{\code{sex}} {\code{Female} or \code{Male} }
#'		\item{\code{Sex}} {\code{0} or \code{1} }
#'		\item{\code{Age}} {Age (years)}
#'		\item{\code{married}} {\code{0} or \code{1}}
#'		\item{\code{Married}} {\code{Not Married} or \code{Married}}
#'		\item{\code{Income}} {Wages and salary for the past 12 months (in \$1,000's)}
#'		\item{\code{HoursWk}} {Hours of work per week}
#'		\item{\code{Race}} {\code{asian}, \code{black}, \code{white}, or \code{other}}
#'		\item{\code{UScitizen}} {\code{Citizen} or \code{Noncitizen}}
#'		\item{\code{USCitizen}} {\code{0} or \code{1}}
#' 		\item{\code{healthInsurance}} {\code{Insured} or \code{Uninsured}}
#' 		\item{\code{HealthInsurance}} {\code{0} or \code{1}}
#'		\item{code{language}} {Native language: \code{English} or \code{Other}}
#'		\item{code{Language}} {Native language: \code{0} or \code{1}}
#' }
# @source The full public dataset can be downloaded at 
# \url{http://www.census.gov/acs/www/data_documentation/pums_data/}
# and the full list of variables is at 
# \url{http://www.census.gov/acs/www/Downloads/data/documentation/pums/DataDict/PUMSDataDict10.pdf}.
#' @source The full public dataset can be downloaded at 
#' \url{http://www.census.gov/programs-surveys/acs/data/pums.html}, and the
#' full list of variables are at 
#' \url{http://www.census.gov/programs-surveys/acs/guidance.html}.
#' @details
#' Several variables in this data set are included in two encodings.  (Watch your 
#' capitalization.)  The lowercase versions have more intuitive codings and can 
#' be used to interpret the numerical codes.
#' @keywords datasets
NA





#' Exercise Hours
#' 
#' Amount of exercise per week for students (and other variables)
#' 
#' Data from an in-class survey of statistics students asking about amount of
#' exercise, TV viewing, handedness, gender, pulse rate, and number of body
#' piercings.
#' 
#' @name ExerciseHours
#' @docType data
#' @format A data frame with 50 observations on the following 7 variables.
#' \itemize{ 
#'    \item{\code{Year}} {Year in school (1=First year,..., 4=Senior)}
#' 
#'    \item{\code{Gender}} { \code{F} or \code{M}} 
#'    \item{\code{Hand}} {Left
#' (\code{l}) or Right (\code{r}) handed?} 
#'    \item{\code{Exercise}} {Hours of
#' exercise per week} 
#'    \item{\code{TV}} {Hours of TV viewing per week}
#' 
#'    \item{\code{Pulse}} {Resting pulse rate (beats per minute)}
#' 
#'    \item{\code{Pierces}} {Number of body piercings} }
#' @source In-class student survey.
#' @keywords datasets
#' @examples
#' 
#' data(ExerciseHours)
#' 
NA





#' Facebook Friends
#' 
#' Data on number of Facebook friends and grey matter density in brain regions
#' related to social perception and associative memory.
#' 
#' A recent study in Great Britain examines the relationship between the number
#' of friends an individual has on Facebook and grey matter density in the
#' areas of the brain associated with social perception and associative memory.
#' The study included 40 students at City University London.
#' 
#' @name FacebookFriends
#' @docType data
#' @format A dataset with 40 observations on the following 2 variables.
#' \itemize{
#'  	\item{\code{GMdensity}} {Normalized z-scores of grey matter
#' 				density in certain brain regions}
#' 		\item{\code{FBfriends}} {Number of friends on Facebook}
#' }
#' @source Kanai, R., Bahrami, B., Roylance, R., and Rees, G., "Online social
#' network size is reflected in human brain structure," Proceedings of the
#' Royal Society, 7 April 2012; 279(1732): 1327-1334. Data approximated from
#' information in the article.
#' @keywords datasets
NA



#' Fat Mice 18
#' 
#' Weight gain for mice with different nighttime light conditions
#' 
#' This is a subset of the LightatNight dataset, showing body mass gain in mice
#' after 4 weeks for two of the treatment conditions: a normal light/dark cycle
#' (LD) or a bright light on at night (LL).
#' 
#' @name FatMice18
#' @docType data
#' @format A dataset with 18 observations on the following 2 variables.
#' \itemize{
#' 	\item{\code{Light}} {Light treatment: \code{LD}= normal
#' 		light/dark cycle or \code{LL}=bright light at night}
#' 	\item{\code{WgtGain4}} {Weight gain (grams over a four week period)}
#' }
#' @source Fonken, L., et. al., "Light at night increases body mass by shifting
#' time of food intake," Proceedings of the National Academy of Sciences,
#' October 26, 2010; 107(43): 18664-18669.
#' @keywords datasets
NA





#' FishGills12
#' 
#' An experiment to look at fish repiration rates in water withdiffernet levels
#' of calcium.
#' 
#' Fish were randomly assigned to twelve tanks with different levels (measured
#' in mg/L) levels of calcium.  Respiration rate was measured as number of gill
#' beats per minute.
#' 
#' @name FishGills12
#' @docType data
#' @format A data frame with 360 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Calcium}} {Amount of calcium in the water (mg/L)}
#' 
#'    \item{\code{GillRate}} {Respiration rate (beats per minute)} }
#' @source Thanks to Pro. Brad Baldwin for supplying the data
#' @keywords datasets
#' @examples
#' 
#' data(FishGills12)
#' 
NA





#' FishGills3
#' 
#' Respiration rate for fish in three levels of calcium.
#' 
#' Fish were randomly assigned to three tanks with different levels (low,
#' medium and high) levels of calcium.  Respiration rate was measured as number
#' of gill beats per minute.
#' 
#' @name FishGills3
#' @docType data
#' @format A data frame with 360 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Calcium}} {Level of Calcium \code{Low} 0.71 mg.L
#' \code{Medium} 5.24 mg/L \code{High} 18.24 mg/L}
#' 
#'    \item{\code{GillRate}} {Respiration rate (beats per minute)} }
#' @source Thanks to Pro. Brad Baldwin for supplying the data
#' @keywords datasets
#' @examples
#' 
#' data(FishGills3)
NA





#' Flight times
#' 
#' Flight times for Flight 179 (Boston-SF) and Flight 180 (SF-Boston)
#' 
#' United Airlines Flight 179 is a daily flight from Boston to San Francisco.
#' Flight 180 goes in the other direction (SF to Boston). The data show the
#' airborn flying times for each flight on the three dates each month (5th,
#' 15th and 25th) in 2010.
#' 
#' @name Flight179
#' @docType data
#' @format A data frame with 36 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Date}} {Date of the flight (5th, 15th and 25th of
#' each month in 2010) as a factor} 
#'    \item{\code{MDY}} {Date as a date object}
#'    \item{\code{Flight179}} {Flying time (Boston-SF) in
#' minutes} 
#'    \item{\code{Flight180}} {Fllying time (SF-Boston) in minutes} }
#' @source Data collected from the Bureau of Transportation Statistics website
#' at
#' \url{http://www.bts.gov/xml/ontimesummarystatistics/src/dstat/OntimeSummaryAirtime.xml}
#' @keywords datasets
#' @examples
#' 
#' data(Flight179)
#' 
NA





#' Florida Lakes
#' 
#' Water quality for a sample of lakes in Florida
#' 
#' This dataset describes characteristics of water and fish samples from 53
#' Florida lakes.  Some variables (e.g. Alkalinity, pH, and Calcium) reflect
#' the chemistry of the water samples. Mercury levels were recorded for a
#' sample of large mouth bass selected at each lake.
#' 
#' @name FloridaLakes
#' @docType data
#' @format A data frame with 53 observations on the following 12 variables.
#' \itemize{ 
#'    \item{\code{ID}} {An identifying number for each lake}
#' 
#'    \item{\code{Lake}} {Name of the lake}
#' 
#'    \item{\code{Alkalinity}} {Concentration of calcium carbonate (in mg/L)}
#' 
#'    \item{\code{pH}} {Acidity} 
#'    \item{\code{Calcium}} {Amount of calcium in
#' water} 
#'    \item{\code{Chlorophyll}} {Amount of chlorophyll in water}
#' 
#'    \item{\code{AvgMercury}} {Average mercury level for a sample of fish (large
#' mouth bass) from each lake} 
#'    \item{\code{NumSamples}} {Number of fish sampled
#' at each lake} 
#'    \item{\code{MinMercury}} {Minimum mercury level in a sampled
#' fish} 
#'    \item{\code{MaxMercury}} {Maximum mercury level in a sampled fish}
#' 
#'    \item{\code{ThreeYrStdMercury}} {Adjusted mercury level to account for the
#' age of the fish} 
#'    \item{\code{AgeData}} {Mean age of fish in each sample} }
#' @source Lange, Royals, and Connor, Transactions of the American Fisheries
#' Society (1993)
#' @keywords datasets
#' @examples
#' 
#' data(FloridaLakes)
#' 
NA





#' Global Internet Usage
#' 
#' Internet usage for several countries
#' 
#' The Nielsen Company measured connection speeds on home computers in nine
#' different countries.  Variables include the percent of internet users with a
#' fast connection (defined as 2Mb/sec or faster) and the average amount of
#' time spent online, defined as total hours connected to the web from a home
#' computer during the month of February 2011.
#' 
#' @name GlobalInternet
#' @docType data
#' @format A data frame with 9 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Country}} {Name of country}
#'    \item{\code{PercentFastConnection}} {Percent of internet users with a fast
#' 			connection} 
#'    \item{\code{HoursOnline}} {Average number of hours online in
#' 			February 2011} 
#' }
#' @source NielsenWire, "Swiss Lead in Speed: Comparing Global Internet
#' Connections", April 1, 2011
#' @keywords datasets
#' @examples
#' 
#' data(GlobalInternet)
#' 
NA





#' GPA and Gender
#' 
#' Data from a survey of introductory statistics students.
#' 
#' This is a subset of the StudentSurvey dataset where cases with missing
#' values have been dropped and gender is coded as a 0/1 indicator variable.
#'
#' 
#' @name GPAGender
#' @docType data
#' @format A dataset with 343 observations on the following 6 variables.
#' \itemize{
#'		\item{\code{Exercise}} {Hours of exercise (per week)}
#' 		\item{\code{SAT}} {Combined SAT scores (out of 1600)}
#'		\item{\code{GPA}} {Grade Point Average (0.00-4.00 scale)}
#'		\item{\code{Pulse}} {Pulse rate (beats per minute)}
#'		\item{\code{Piercings}} {Number of body piercings}
#' 		\item{\code{Sex}} {\code{Female} or \code{Male}}
#' }
#' @source A first day survey over several different introductory statistics
#' classes.
#' @keywords datasets
NA





#' Happy Planet Index
#' 
#' Measurements related to happiness and well-being for 143 countries.
#' 
#' Data for 143 countries from the Happy Planet Index Project that works to
#' quantify indicators of happiness, well-being, and ecological footprint at a
#' country level.
#' 
#' @name HappyPlanetIndex
#' @docType data
#' @format A dataset with 143 observations on the following 11 variables.
#' \itemize{
#'	   \item{\code{Country}} {Name of country}
#'     \item{\code{Region}} {Three-digit country code}
#'     \item{\code{Happiness}} {Score on a 0-10 scale for  
#'           average level of happiness (10 is happiest)}
#'     \item{\code{LifeExpectancy}} {Average life expectancy (in years)}
#'     \item{\code{Footprint}} {Ecological footprint -- a measure of the 
#'				(per capita) ecological impact}
#'     \item{\code{HLY}} {Happy Life Years - combines life expectancy with well-being}
#'     \item{\code{HPI}} {Happy Planet Index (0-100 scale)}
#'     \item{\code{HPIRank}} {HPI rank for the country}
#'     \item{\code{GDPperCapita}} {Gross Domestic Product (per capita)}
#'     \item{\code{HDI}} {Human Development Index}
#'     \item{\code{Population}} {Population (in millions)}
#' }
#' @references Marks, N., "The Happy Planet Index", www.TED.com/talks, August
#' 29, 2010.
#' @source Data downloaded from \url{http://www.happyplanetindex.org/data/}
#' @keywords datasets
NA





#' Hockey Penalties
#' 
#' Penalty minutes (per game) for NHL teams in 2010-11
#' 
#' Data give the average numeber of penalty minutes for each of the 30 National
#' HOckey League (NHL) teams during the 2010-11 regualar season.
#' 
#' @name HockeyPenalties
#' @docType data
#' @format A data frame with 30 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Team}} {Name of the team}
#' 
#'    \item{\code{PIMperG}} {Average penaly minutes per game} }
#' @source Data obtained online at www.nhl.com
#' @keywords datasets
#' @examples
#' 
#' data(HockeyPenalties)
#' 
NA





#' Hollywood Movies in 2011
#' 
#' 
#' @name HollywoodMovies2011
#' @docType data
#' @format A data frame with 136 observations on the following 14 variables.
#' \itemize{ 
#'    \item{\code{Movie}} {a factor with many levels}
#' 
#'    \item{\code{LeadStudio}} {a factor with many levels}
#' 
#'    \item{\code{RottenTomatoes}} {a numeric vector}
#' 
#'    \item{\code{AudienceScore}} {a numeric vector} 
#'    \item{\code{Story}} {a factor
#' with many levels} 
#'    \item{\code{Genre}} {a factor with levels \code{Action}
#' \code{Adventure} \code{Animation} \code{Comedy} \code{Drama} \code{Fantasy}
#' \code{Horror} \code{Romance} \code{Thriller}}
#' 
#'    \item{\code{TheatersOpenWeek}} {a numeric vector}
#' 
#'    \item{\code{BOAverageOpenWeek}} {a numeric vector}
#' 
#'    \item{\code{DomesticGross}} {a numeric vector} 
#'    \item{\code{ForeignGross}} {a
#' numeric vector} 
#'    \item{\code{WorldGross}} {a numeric vector}
#' 
#'    \item{\code{Budget}} {a numeric vector} 
#'    \item{\code{Profitability}} {a
#' numeric vector} 
#'    \item{\code{OpeningWeekend}} {a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(HollywoodMovies2011)
#' ## maybe str(HollywoodMovies2011) ; plot(HollywoodMovies2011) ...
#' 
NA





#' Home for Sale
#' 
#' Data on homes for sale in four states
#' 
#' Data for samples of homes for sale in each state, selected from zillow.com.
#' 
#' @name HomesForSale
#' @docType data
#' @format A data frame with 120 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{State}} {Location of the home: \code{CA} \code{NJ}
#' \code{NY} \code{PA}} 
#'    \item{\code{Price}} {Asking price (in \$1,000's)}
#' 
#'    \item{\code{Size}} {Area of all rooms (in 1,000's sq. ft.)}
#' 
#'    \item{\code{Beds}} {Number of bedrooms} 
#'    \item{\code{Baths}} {Number of
#' bathrooms} }
#' @references Data collected from www.zillow.com in 2010.
#' @keywords datasets
#' @examples
#' 
#' data(HomesForSale)
#' 
NA





#' Home for Sale in California
#' 
#' Data for a sample of homes offered for sale in California
#' 
#' This is a subset of the HomesForSale data with just information from homes
#' in California (CA). Data were collected offerings listed on an online site.
#' 
#' @name HomesForSaleCA
#' @docType data
#' @format A data frame with 30 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{State}} {Location of the home: \code{CA}}
#' 
#'    \item{\code{Price}} {Asking price (in \$1,000's)} 
#'    \item{\code{Size}} {Area of
#' all rooms (in 1,000's sq. ft.)} 
#'    \item{\code{Beds}} {Number of bedrooms}
#' 
#'    \item{\code{Baths}} {Number of bathrooms} }
#' @source Data collected from www.zillow.com in 2010.
#' @keywords datasets
#' @examples
#' 
#' data(HomesForSaleCA)
#' 
NA





#' Homes for sale in Canton, NY
#' 
#' Prices of homes for sale in Canton, NY
#' 
#' Data for samples of homes for sale in Canton, NY, selected from zillow.com.
#' 
#' @name HomesForSaleCanton
#' @docType data
#' @format A data frame with 10 observations on the following variable.
#' \itemize{ 
#'    \item{\code{Price}} {Asking price for the home (in \$1,000's)} }
#' @source Data collected from www.zillow.com in 2010.
#' @keywords datasets
#' @examples
#' 
#' data(HomesForSaleCanton)
#' 
NA





#' Home for Sale in New York
#' 
#' Data for a sample of homes offered for sale in New York State
#' 
#' This is a subset of the HomesForSale data with just information from homes
#' in New York State (NY). Data were collected offerings listed on an online
#' site.
#' 
#' @name HomesForSaleNY
#' @docType data
#' @format A data frame with 30 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{State}} {Location of the home: \code{NY}}
#' 
#'    \item{\code{Price}} {Asking price (in \$1,000's)} 
#'    \item{\code{Size}} {Area of
#' all rooms (in 1,000's sq. ft.)} 
#'    \item{\code{Beds}} {Number of bedrooms}
#' 
#'    \item{\code{Baths}} {Number of bathrooms} }
#' @source Data collected from www.zillow.com in 2010.
#' @keywords datasets
#' @examples
#' 
#' data(HomesForSaleNY)
#' 
NA





#' Honeybee Circuits
#' 
#' Number of circuits for honeybee dances and nest quality
#' 
#' When honeybees are looking for a new home, they send out scouts to explore
#' options. When a scout returns, she does a "waggle dance" with multiple
#' circuit repetitions to tell the swarm about the option she found. The bees
#' then decide between the options and pick the best one. Scientists wanted to
#' find out how honeybees decide which is the best option, so they took a swarm
#' of honeybees to an island with only two possible options for new homes: one
#' of very high honeybee quality and one of low quality. They then kept track
#' of the scouts who visited each option and counted the number of waggle dance
#' circuits each scout bee did when describing the option.
#' 
#' @name HoneybeeCircuits
#' @docType data
#' @format A dataset with 78 observations on the following 2 variables.
#' \itemize{
#'		\item{\code{Circuits}} {Number of waggle dance circuits for a
#' 					returning scout bee}
#'		\item{\code{Quality}} {Quality of the nest site: \code{High} or \code{Low}}
#' }
#' @source Seeley, T., Honeybee Democracy, Princeton University Press,
#' Princeton, NJ, 2010, p. 128
#' @keywords datasets
NA





#' Honeybee Waggle
#' 
#' Honeybee dance duration and distance to nesting site
#' 
#' When honeybee scouts find a food source or a nice site for a new home, they
#' communicate the location to the rest of the swarm by doing a "waggle dance."
#' They point in the direction of the site and dance longer for sites farther
#' away. The rest of the bees use the duration of the dance to predict distance
#' to the site.
#' 
#' @name HoneybeeWaggle
#' @docType data
#' @format A dataset with 7 observations on the following 2 variables.
#' \itemize{
#'		\item{\code{Distance}}{Distance to the potential nest site (in meters)} 
#'		\item{\code{Duration}} {Duration of the waggle dance (in seconds)}
#' }
#' @source Seeley, T., Honeybee Democracy, Princeton University Press,
#' Princeton, NJ, 2010, p. 128
#' @keywords datasets
NA





#' Hot Dog Eating Contest
#' 
#' Winning number of hot dogs consumed in an eating contest
#' 
#' Every Fourth of July, Nathan's Famous in New York City holds a hot dog
#' eating contest, in which contestants try to eat as many hot dogs (with buns)
#' as possible in ten minutes.  The winning number of hot dogs are given for
#' each year from 2002-2011.
#' 
#' @name HotDogs
#' @docType data
#' @format A dataset with 10 observations on the following 2 variables.
#' \itemize{
#'  	\item{\code{Year}} {Year of the contest: 2002-2011}
#' 		\item{\code{HotDogs}} {Winning number of hot dogs consumed }
#' }
#' @source Nathan's Famous webste at
#' \url{http://nathansfamous.com/contest/hall_of_fame}
#' @keywords datasets
NA





#' Intensive Care Unit Admissions
#' 
#' data from patients admitted to an intensive care unit
#' 
#' Data from a sample of 200 patients following admission to an adult intensive
#' care unit (ICU).
#' 
#' @name ICUAdmissions
#' @docType data
#' @format A data frame with 200 observations on the following 21 variables.
#' \itemize{ 
#'    \item{\code{ID}} {Patient ID number}
#' 
#'    \item{\code{status}} {Patient status: \code{Lived} or \code{Died}}
#'    \item{\code{Status}} {numerical code for \code{Status}}
#' 
#'    \item{\code{Age}} {Patient's age (in years)}
#' 
#'    \item{\code{sex}} {\code{Male} or \code{Female}}
#'    \item{\code{Sex}} {numerical code for \code{sex}}
#' 
#'    \item{\code{race}} {Patient's race: \code{White}, \code{Black}, or
#' \code{Other} }
#'    \item{\code{Race}} {numerical code for \code{race}}
#'    \item{\code{service}} {Type of service: \code{Medical} or
#' \code{Surgical} }
#'    \item{\code{Service}} {numerical code for \code{service}}
#'    \item{\code{cancer}} {Is cancer involved?  \code{No} or \code{Yes}}
#'    \item{\code{Cancer}} {Is cancer involved?  \code{0} or \code{1}}
#'    \item{\code{renal}} {Is chronic renal failure involved? \code{No} or \code{Yes}} 
#'    \item{\code{Renal}} {Is chronic renal failure involved? \code{0} or \code{1}} 
#'    \item{\code{infection}} {Is infection involved? \code{No} or \code{Yes}}
#'    \item{\code{Infection}} {Is infection involved? \code{0} or \code{1}}
#'    \item{\code{cpr}} {Patient gets CPR prior to admission? \code{No} or \code{Yes}}
#'    \item{\code{CPR}} {Patient gets CPR prior to admission? \code{0} or \code{1}}
#'    \item{\code{Systolic}} {Systolic blood pressure (in mm of Hg)} 
#'    \item{\code{HeartRate}} {Pulse rate (beats per minute)} 
#'    \item{\code{previous}} {Previous admission to ICU wihtin 6 months? \code{No} or \code{Yes}}
#'    \item{\code{Previous}} {Previous admission to ICU wihtin 6 months? \code{0} or \code{1}}
#'    \item{\code{type}} {Admission type: \code{Elective} or \code{Emergency}}
#'    \item{\code{Type}} {Admission type: \code{0} or \code{1}}
#'    \item{\code{fracture}} {Fractured bone involved? \code{No} or \code{Yes}}
#'    \item{\code{Fracture}} {Fractured bone involved? \code{0} or \code{1}}
#'    \item{\code{pO2}} {Partial oxygen level from blood gases under 60? 
#'                      \code{No} or  \code{Yes}}
#'    \item{\code{PO2}} {Partial oxygen level from blood gases under 60? 
#'                      \code{0} or  \code{1}}
#'    \item{\code{pHlow}} {pH from blood gas under 7.25? \code{No} or \code{Yes}}
#'    \item{\code{pH}} {pH from blood gas under or over 7.25? \code{Low} or \code{Hi}}
#'    \item{\code{PH}} {pH from blood gas under or over 7.25? \code{0} or \code{1}}
#'    \item{\code{pCO2hi}} {Partial carbon dioxide level from blood gas over 45?
#'                \code{No} or \code{Yes}}
#'    \item{\code{pCO2}} {Partial carbon dioxide level from blood gas over or under 45?
#'                \code{Low} or \code{Hi}}
#'    \item{\code{PCO2}} {Partial carbon dioxide level from blood gas over or under 45?
#'                \code{0} or \code{1}}
#'    \item{\code{bicarbonateLow}} {Bicarbonate from
#' blood gas under 18? \code{No} or \code{Yes}}
#'    \item{\code{bicarbonate}} {Bicarbonate from
#' blood gas under or over 18? \code{Low} or \code{Hi}}
#'    \item{\code{Bicarbonate}} {Bicarbonate from
#' blood gas under or over 18? \code{0} or \code{1}}
#'    \item{\code{creatinineHi}} {Creatinine from blood gas over 2.0? 
#'          \code{No} or \code{Yes}}
#'    \item{\code{creatinine}} {Creatinine from blood gas over or under 2.0? 
#'          \code{Low} or \code{Hi}}
#'    \item{\code{Creatinine}} {Creatinine from blood gas over or under 2.0? 
#'          \code{0} or \code{1}}
#'    \item{\code{consciousness}} {Levels: \code{Conscious}, \code{Deep Stupor}, or \code{Coma}}
#'    \item{\code{Consciousness}} {Levels: \code{0}, \code{1}, or \code{2}}
#'    }
#' @source DASL dataset downloaded from
#' \url{http://lib.stat.cmu.edu/DASL/Datafiles/ICU.html}
#' @keywords datasets
#' @examples
#' 
#' data(ICUAdmissions)
#' 
NA





#' Immune Tea
#' 
#' Interferon gamma production and tea drinking
#' 
#' Eleven healthy non-tea-drinking individuals were asked to drink five or six
#' cups of tea a day, while ten healthy non-tea and non-coffee-drinkers were
#' asked to drink the same amount of coffee, which has caffeine but not the
#' L-theanine that is in tea. The groups were randomly assigned. After two
#' weeks, blood samples were exposed to an antigen and production of interferon
#' gamma was measured.
#' 
#' @name ImmuneTea
#' @docType data
#' @format A data frame with 21 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{InterferonGamma}} {Measure of interferon gamma
#' production} 
#'    \item{\code{Drink}} {Type of drink: \code{Coffee} or \code{Tea}}
#' }
#' @source Adapted from Kamath, et.al., "Antigens in tea-Beverage prime human V
#' 2V2 T cells in vitro and in vivo for memory and non-memory antibacterial
#' cytokine responses", Proceedings of the National Academy of Sciences, May
#' 13, 2003.
#' @keywords datasets
#' @examples
#' 
#' data(ImmuneTea)
#' 
NA





#' Inkjet Printers
#' 
#' Data from online reviews of inkjet printers
#' 
#' Information from reviews of inkjet printers at PCMag.com in August 2011.
#' 
#' @name InkjetPrinters
#' @docType data
#' @format A dataset with 20 observations on the following 6 variables.
#' \itemize{
#'		\item{\code{Model}} {Model name of printer}
#'		\item{\code{PPM}} {Printing rate (pages per minute) for a 
#'			benchmark set of print jobs}
#' 		\item{\code{PhotoTime}} {Time (in seconds) to print 4x6 color photos}
#' 		\item{\code{Price}} {Typical retail price (in dollars)}
#'		\item{\code{CostBW}} {Cost per page (in cents) for printing in 
#'			black & white}
#'		\item{\code{CostColor}} {Cost per page (in cents) for printing in color}
#' }
#' @source Inkjet printer reviews found at
#' \url{http://www.pcmag.com/reviews/printers}, August 2011.
#' @keywords datasets
NA





#' Life Expectancy and Vehicle Registrations
#' 
#' 
#' @name LifeExpectancyVehicles
#' @docType data
#' @format A data frame with 40 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Year}} {a numeric vector}
#' 
#'    \item{\code{LifeExpectancy}} {a numeric vector} 
#'    \item{\code{Vehicles}} {a
#' numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(LifeExpectancyVehicles)
#' ## maybe str(LifeExpectancyVehicles) ; plot(LifeExpectancyVehicles) ...
#' 
NA





#' Light at Night for Mice
#' 
#' Data from an experiment with mice having different nighttime light
#' conditions
#' 
#' In this study, 27 mice were randomly split into three groups. One group was
#' on a normal light/dark cycle (LD), one group had bright light on all the
#' time (LL), and one group had light during the day and dim light at night
#' (DM). The dim light was equivalent to having a television set on in a room.
#' The mice in darkness ate most of their food during their active (nighttime)
#' period, matching the behavior of mice in the wild. The mice in both dim
#' light and bright light, however, consumed more than half of their food
#' during the well-lit rest period, when most mice are sleeping.
#' 
#' @name LightatNight
#' @docType data
#' @format A dataset with 27 observations on the following 9 variables.
#' \itemize{
#'		\item{\code{Light}} {\code{DM}=dim light at night, 
#'				\code{LD}=dark at night, or \code{LL}=bright light at night}
#'		\item{\code{BMGain}} {Body mass gain (in grams over a four week period)}
#'		\item{\code{Corticosterone}} {Blood corticosterene level (a measure of stress)}
#'		\item{\code{DayPct}} {Percent of calories eaten during the day}
#'		\item{\code{Consumption}} {Daily food consumption (grams)}
#'		\item{\code{GlucoseInt}} {Glucose intolerant? \code{No} or \code{Yes}}
#'		\item{\code{GTT15}} {Glucose level in the blood 15 minutes
#' 			after a glucose injection}
#'		\item{\code{GTT120}} {Glucose level in the blood 120 minutes after a 
#'				glucose injection}
#'		\item{\code{Activity}} {A measure of physical activity level}
#' }
#' @source Fonken, L., et. al., "Light at night increases body mass by shifting
#' time of food intake," Proceedings of the National Academy of Sciences,
#' October 26, 2010; 107(43): 18664-18669.
#' @keywords datasets
NA





#' Malevolent Uniforms NFL
#' 
#' Perceived malevolence of uniforms and penalies for National Football League
#' (NFL) teams
#' 
#' Participants with no knowledge of the teams rated the jerseys on
#' characteristics such as timid/aggressive, nice/mean and good/bad. The
#' averages of these responses produced a "malevolence" index with higher
#' scores signifying impressions of more malevolent uniforms. To measure
#' aggressiveness, the authors used the amount of penalty yards converted to
#' z-scores and averaged for each team over the seasons from 1970-1986.
#' 
#' @name MalevolentUniformsNFL
#' @docType data
#' @format A data frame with 28 observations on the following 7 variables.
#' \itemize{ 
#'    \item{\code{NFLTeam}} {Team name}
#' 
#'    \item{\code{NFL_Malevolence}} {Score reflecting the "malevolence" of a
#' team's uniform} 
#'    \item{\code{ZPenYds}} {Z-score for penalty yards} }
#' @source Frank and Gilovich, "The Dark Side of Self- and Social Perception:
#' Black Uniforms and Aggression in Professional Sports", Journal of
#' Personality and Social Psychology, Vol. 54, No. 1, 1988, p. 74-85.
#' @keywords datasets
#' @examples
#' 
#' data(MalevolentUniformsNFL)
#' 
NA





#' Malevolent Uniforms NHL
#' 
#' Perceived malevolence of uniforms and penalies for National Hockey League
#' (NHL) teams
#' 
#' Participants with no knowledge of the teams rated the jerseys on
#' characteristics such as timid/aggressive, nice/mean and good/bad. The
#' averages of these responses produced a "malevolence" index with higher
#' scores signifying impressions of more malevolent uniforms. To measure
#' aggressiveness, the authors used the amount of penalty minutes converted to
#' z-scores and averaged for each team over the seasons from 1970-1986.
#' 
#' @name MalevolentUniformsNHL
#' @docType data
#' @format A data frame with 28 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{NHLTeam}} {Team name}
#' 
#'    \item{\code{NHL_Malevolence}} {Score reflecting the "malevolence" of a
#' team's uniform} 
#'    \item{\code{ZPenMin}} {Z-score for penalty minutes} }
#' @source Frank and Gilovich, "The Dark Side of Self- and Social Perception:
#' Black Uniforms and Aggression in Professional Sports", Journal of
#' Personality and Social Psychology, Vol. 54, No. 1, 1988, p. 74-85.
#' @keywords datasets
#' @examples
#' 
#' data(MalevolentUniformsNHL)
#' 
NA





#' Mammal Longevity
#' 
#' Longevity and gestation period for mammals
#' 
#' Dataset with average lifespan (in years) and typical gestation period (in
#' days) for 40 different species of mammals.
#' 
#' @name MammalLongevity
#' @docType data
#' @format A data frame with 40 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Animal}} {Species of mammal}
#' 
#'    \item{\code{Gestation}} {Time from fertilization until birth (in days)}
#' 
#'    \item{\code{Longevity}} {Average lifespan (in years)} }
#' @source 2010 World Almanac, pg. 292.
#' @keywords datasets
#' @examples
#' 
#' data(MammalLongevity)
#' 
NA





#' Manhattan Apartment Prices
#' 
#' Monthly rent for one-bedroom apartments in Manhattan, NY
#' 
#' Monthly rents for a sample of 20 one-bedroom apartments in Manhattan, NY
#' that were advertised on Craig's List in July, 2011.
#' 
#' @name ManhattanApartments
#' @docType data
#' @format A dataset with 20 observations on the following variable.
#' \itemize{
#'		\item{\code{Rent}} {Montly rent in dollars}
#' }
#' @source Apartments advertised on Craig's List at \url{http://newyork.craigslist.org},
#' July 5, 2011.
#' @keywords datasets
NA





#' Marriage Ages
#' 
#' Ages for husbands and wives from marriage licenses
#' 
#' Data from a sample of 100 marriage licences in St. Lawrence County, NY gives
#' the ages of husbands and wives for newly married couples.
#' 
#' @name MarriageAges
#' @docType data
#' @format A data frame with 100 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Husband}} {Age of husband at marriage}
#' 
#'    \item{\code{Wife}} {Age of wife at marriage} }
#' @source Thanks to Linda Casserly, St. Lawrence County Clerk's Office
#' @keywords datasets
#' @examples
#' 
#' data(MarriageAges)
#' 
NA





#' Masters Golf Scores
#' 
#' Scores from the 2011 Masters golf tournament
#' 
#' Data for a random sample of 20 golfers who made the cut at the 2011 Masters
#' golf tournament.
#' 
#' @name MastersGolf
#' @docType data
#' @format A dataset with 20 observations on the following 2 variables.
#' \itemize{
#' 		\item{\code{First}} {First round score (in relation to par)}
#' 		\item{\code{Final}} {Final four round score (in relation to par)}
#' }
# This URL is dead now, so removing.
# @source 2011 Masters tournament results at
# \url{http://www.masters.com/en_US/discover/past_winners.html}
# @keywords datasets
NA





#' Mental Muscle
#' 
#' Comparing actual movements to mental imaging movements
#' 
#' In this study, participants were asked to either perform actual arm pointing
#' motions or to mentally imagine equivalent arm pointing motions. Participants
#' then developed muscle fatigue by holding a heavy weight out horizontally as
#' long as they could. After becoming fatigued, they were asked to repeat the
#' previous mental or actual motions. Eight participants were assigned to each
#' group, and the time in seconds to complete the motions was measured before
#' and after fatigue.
#' 
#' @name MentalMuscle
#' @docType data
#' @format A dataset with 32 observations on the following 3 variables.
#' \itemize{
#' 		\item{\code{Action}} {Treatment: \code{Actual} motions or
#' 				\code{Mental} imaging motions}
#'		\item{\code{PreFatigue}} {Time (in seconds) to  complete motions before fatigue}
#'		\item{\code{PostFatigue}} {Time (in seconds) to complete motions after fatigue}
#' }
#' @source Data approximated from summary statistics in: Demougeot L. and
#' Papaxanthis C., "Muscle Fatigue Affects Mental Simulation of Action," The
#' Journal of Neuroscience, July 20, 2011, 31(29):10712-10720.
#' @keywords datasets
NA





#' Miami Heat basketball
#' 
#' Game log data for the Miami Heat basketball team in 2010-11
#' 
#' Information from online boxscores for all 82 regular season games payed by
#' the Miami Heat basketball team during the 2010-11 regular season.
#' 
#' @name MiamiHeat
#' @docType data
#' @format A data frame with 82 observations on the following 33 variables.
#' \itemize{ 
#'    \item{\code{Game}} {ID number for each game}
#' 
#'    \item{\code{MDY}} {Data the game was played as a date object}
#'    \item{\code{Date}} {Data the game was played as a character string}
#' 
#'    \item{\code{Location}} {\code{Away} or \code{Home}}
#' 
#'    \item{\code{Opp}} {Opponent tream} 
#'    \item{\code{Win}} {Game result: \code{L}
#' or \code{W}} 
#'    \item{\code{FG}} {Field goals made} 
#'    \item{\code{FGA}} {Field
#' goals attempted} 
#'    \item{\code{FG3}} {Three-point field goals made}
#' 
#'    \item{\code{FG3A}} {Three-point field goals attempted}
#' 
#'    \item{\code{FT}} {Free throws made} 
#'    \item{\code{FTA}} {Free throws
#' attempted} 
#'    \item{\code{Rebounds}} {Total rebounds}
#' 
#'    \item{\code{OffReb}} {Offensive rebounds} 
#'    \item{\code{Assists}} {Number of
#' assists} 
#'    \item{\code{Steals}} {Number of steals}
#' 
#'    \item{\code{Blocks}} {Number of shots blocked}
#' 
#'    \item{\code{Trunovers}} {Number of turnovers} 
#'    \item{\code{Fouls}} {Number of
#' fouls} 
#'    \item{\code{Points}} {Number of points scored}
#' 
#'    \item{\code{OppFG}} {Opponet's field goals made}
#' 
#'    \item{\code{OppFGA}} {Opponent's Field goals attempted}
#' 
#'    \item{\code{OppFG3}} {Opponent's Three-point field goals made}
#' 
#'    \item{\code{OppFG3A}} {Opponent's Three-point field goals attempted}
#' 
#'    \item{\code{OppFT}} {Opponent's Free throws made}
#' 
#'    \item{\code{OppFTA}} {Opponent's Free throws attempted}
#' 
#'    \item{\code{OppOffReb}} {Opponent's Total rebounds}
#' 
#'    \item{\code{OppRebounds}} {Opponent's Offensive rebounds}
#' 
#'    \item{\code{OppAssists}} {Opponent's assists}
#' 
#'    \item{\code{OppSteals}} {Opponent's steals}
#' 
#'    \item{\code{OppBlocks}} {Opponent's shots blocked}
#' 
#'    \item{\code{OppTurnovers}} {Opponent's turnovers}
#' 
#'    \item{\code{OppFouls}} {Opponent's fouls}
#' 
#'    \item{\code{OppPoints}} {Opponent's points scored} }
#' @source Data for the 2010-11 Miami games downloaded from
#' \url{http://www.basketball-reference.com/teams/MIA/2011/gamelog/}
#' @keywords datasets
#' @examples
#' 
#' data(MiamiHeat)
#' 
NA





#' Mindset Matters
#' 
#' Data from a study of perceived exercise with maids
#' 
#' In 2007 a Harvard psychologist recruited 75 female maids working in
#' different hotels to participate in a study. She informed 41 maids (randomly
#' chosen) that the work they do satisfies the Surgeon General's
#' recommendations for an active lifestyle (which is true), giving them
#' examples for how why their work is good exercise. The other 34 maids were
#' told nothing (uninformed). Various chacteristics (weight, body mass index,
#' ...) were recorded for each subject at the start of the experiment and again
#' four weeks later. Maids with missing values for weight change have been
#' removed.
#' 
#' @name MindsetMatters
#' @docType data
#' @format A data frame with 75 observations on the following 14 variables.
#' \itemize{ 
#'    \item{\code{Condition}} {Treatment condition: \code{uninformed} or
#' \code{informed}} 
#'    \item{\code{Cond}} {Treatment condition: \code{0}=uninformed or
#' \code{1}=informed} 
#'    \item{\code{Age}} {Age (in years)}
#' 
#'    \item{\code{Wt}} {Original weight (in pounds)} 
#'    \item{\code{Wt2}} {Weight
#' after 4 weeks (in pounds} 
#'    \item{\code{BMI}} {Original body mass index}
#' 
#'    \item{\code{BMI2}} {Body mass index after 4 weeks}
#' 
#'    \item{\code{Fat}} {Original body fat percentage} 
#'    \item{\code{Fat2}} {Body
#' fat percentage after 4 weeks} 
#'    \item{\code{WHR}} {Original waist to hip
#' ratio} 
#'    \item{\code{WHR2}} {Waist to hip ratio} 
#'    \item{\code{Syst}} {Original
#' systolic blood pressure} 
#'    \item{\code{Syst2}} {Systolic blood pressure after
#' 4 weeks} 
#'    \item{\code{Diast}} {Original diastolic blood pressure}
#' 
#'    \item{\code{Diast2}} {Diastolic blood pressure after 4 weeks} }
#' @source Crum, A.J. and Langer, E.J. (2007). Mind-Set Matters: Exercise and
#' the Placebo Effect, Psychological Science, 18:165-171. Thanks to the authors
#' for supplying the data.
#' @keywords datasets
#' @examples
#' 
#' data(MindsetMatters)
#' 
NA





#' Ministers and Rum
#' 
#' Number of ministers and rum imports to New England between 1860 and 1940
#' 
#' Data every five years from 1860 to 1940 on he number of Methodist ministers
#' working in New England and the annual rum imports into Boston.
#' 
#' @name MinistersRum
#' @docType data
#' @format A data frame with 17 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Year}} {Every five years from 1860-1940}
#' 
#'    \item{\code{Ministers}} {Number of Methodist ministers in New England}
#' 
#'    \item{\code{Rum}} {Number of barrels of rum imported into Boston} }
#' @keywords datasets
#' @examples
#' 
#' data(MinistersRum)
#' 
NA





#' Mustang Prices
#' 
#' Price, age, and mileage for used Mustangs at an internet website
#' 
#' A statistics student, Gabe McBride, was interested in prices for used
#' Mustang cars being offered for sale on an internet site.  He sampled 25 cars
#' from the website and recorded the age (in years), mileage (in thousands of
#' miles) and asking price (in \$1,000's) for each car in his sample.
#' 
#' @name MustangPrice
#' @docType data
#' @format A data frame with 25 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Age}} {Age of the car (in years)}
#' 
#'    \item{\code{Miles}} {Mileage on the car (in 1,000's)}
#' 
#'    \item{\code{Price}} {Asking price in (\$1,000's)} }
#' @source Student project with data collected from autotrader.com in 2008.
#' @keywords datasets
#' @examples
#' 
#' data(MustangPrice)
#' 
NA





#' NBA Players data for 2010-11 Season
#' 
#' Data from the 2010-2011 regular season for 176 NBA basketball players.
#' 
#' 
#' @name NBAPlayers2011
#' @docType data
#' @format A data frame with 176 observations on the following 25 variables.
#' \itemize{ 
#'    \item{\code{Player}} {Name of player} 
#'    \item{\code{Age}} {Age in
#' years)} 
#'    \item{\code{Team}} {Team name} 
#'    \item{\code{Games}} {Games played
#' (out of 82)} 
#'    \item{\code{Starts}} {Games started}
#' 
#'    \item{\code{Mins}} {Minutes played} 
#'    \item{\code{MinPerGame}} {Minutes per
#' game} 
#'    \item{\code{FGMade}} {Field goals made} 
#'    \item{\code{FGAttempt}} {Field
#' goals attempted} 
#'    \item{\code{FGPct}} {Three-point field goal percentage}
#' 
#'    \item{\code{FG3Made}} {Three-point field goals made}
#' 
#'    \item{\code{FG3Attempt}} {Three-point field goals attempted}
#' 
#'    \item{\code{FG3Pct}} {Field goal percentage} 
#'    \item{\code{FTMade}} {Free
#' throws made} 
#'    \item{\code{FTAttempt}} {Free throws attempted}
#' 
#'    \item{\code{FTPct}} {Free throw percentage}
#' 
#'    \item{\code{OffRebound}} {Offensive rebounds}
#' 
#'    \item{\code{DefRebound}} {Defensive rebounds} 
#'    \item{\code{Rebounds}} {Total
#' rebounds} 
#'    \item{\code{Assists}} {Number of assists}
#' 
#'    \item{\code{Steals}} {Number of steals} 
#'    \item{\code{Blocks}} {Number of
#' blocked shots} 
#'    \item{\code{Turnovers}} {Number of turnovers}
#' 
#'    \item{\code{Fouls}} {Number of personal fouls} 
#'    \item{\code{Points}} {Number
#' of points scored} }
#' @source \url{http://www.basketball-reference.com/leagues/NBA_2011_stats.html}
#' @keywords datasets
#' @examples
#' 
#' data(NBAPlayers2011)
#' 
NA





#' NBA 2010-11 Regular Season Standings
#' 
#' Won-Loss record and statistics for NBA Teams
#' 
#' Won-Loss record and regular season statistics for 30 teams in the National
#' Basketball Association for the 2010-2011 season
#' 
#' @name NBAStandings
#' @docType data
#' @format A data frame with 30 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{Team}} {Team name} 
#'    \item{\code{Wins}} {Number of wins
#' in an 82 game regular season} 
#'    \item{\code{Losses}} {Number of losses}
#' 
#'    \item{\code{WinPct}} {Proportion of games won} 
#'    \item{\code{PtsFor}} {Average
#' points scored per game} 
#'    \item{\code{PtsAgainst}} {Average points allowed per
#' game} }
#' @source Data downloaded from
#' \url{http://www.basketball-reference.com/leagues/NBA_2011_games.html}
#' @keywords datasets
#' @examples
#' 
#' data(NBAStandings)
#' 
NA





#' NFL Game Scores in 2011
#' 
#' Results for all NFL games for the 2011 regular season
#' 
#' Data for all 256 regular season games in the National Football League (NFL)
#' for the 2011 season.
#' 
#' @name NFLScores2011
#' @docType data
#' @format A dataset with 256 observations on the following 11 variables.
#' \itemize{
#' \item{\code{Week}} {a numeric vector}
#' \item{\code{HomeTeam}} {Home team name}
#' \item{\code{AwayTeam}} {Visiting team name}
#' \item{\code{HomeScore}} {Points scored by the home team}
#' \item{\code{AwayScore}} {Points scored by the visiting team}
#' \item{\code{HomeYards}} {Yards gained by the home team}
#' \item{\code{AwayYards}} {Yards gained by the visiting team}
#' \item{\code{HomeTO}} {Turnovers lost by the home team}
#' \item{\code{AwayTO}} {Turnovers lost by the visiting team}
#' \item{\code{Date}}  {Date of the game (as a character string)}
#' \item{\code{YDM}}  {Date of the game (as a date object)}
#' \item{\code{Day}} {Day of the week: \code{Mon}, \code{Sat}, \code{Sun}, or \code{Thu}}
#' }
#' @source NFL scores and game statistics found at
#' \url{http://www.pro-football-reference.com/years/2011/games.htm}.
#' @keywords datasets
#' 

NA





#' Nutrition Study
#' 
#' Variables related to nutrition and health for 315 individuals
#' 
#' Data from a cross-sectional study to investigate the relationship between
#' personal characteristics and dietary factors, and plasma concentrations of
#' retinol, beta-carotene and other carotenoids. Study subjects were patients
#' who had an elective surgical procedure during a three-year period to biopsy
#' or remove a lesion of the lung, colon, breast, skin, ovary or uterus that
#' was found to be non-cancerous.
#' 
#' @name NutritionStudy
#' @docType data
#' @format A data frame with 315 observations on the following 17 variables.
#' \itemize{ 
#'    \item{\code{ID}} {ID number for each subject in this sample}
#' 
#'    \item{\code{Age}} {Subject's age (in years)} 
#'    \item{\code{Smoke}} {a factor
#' with levels \code{No} \code{Yes}} 
#'    \item{\code{Quetelet}} {Weight/(Height^2)}
#' 
#'    \item{\code{Vitamin}} {Vitamin use: \code{1}=Regular,
#' \code{2}=Occasional, or \code{3}=No} 
#'    \item{\code{Calories}} {Number of
#' calories consumed per day} 
#'    \item{\code{Fat}} {Grams of fat consumed per day}
#' 
#'    \item{\code{Fiber}} {Grams of fiber consumed per day}
#' 
#'    \item{\code{Alcohol}} {Number of alcoholic drinks consumed per week}
#' 
#'    \item{\code{Cholesterol}} {Cholesterol consumed (mg per day)}
#' 
#'    \item{\code{BetaDiet}} {Dietary beta-carotene consumed (mcg per day)}
#' 
#'    \item{\code{RetinolDiet}} {Dietary retinol consumed (mcg per day)}
#' 
#'    \item{\code{BetaPlasma}} {Plasma beta-carotene (ng/ml)}
#' 
#'    \item{\code{RetinolPlasma}} {Plasma retinol (ng/ml)}
#' 
#'    \item{\code{Sex}} {Cosed as \code{Female} or \code{Male}}
#' 
#'    \item{\code{VitaminUse}} {Coded as \code{No} \code{Occasional}
#' \code{Regular}} 
#'    \item{\code{EverSmoke}} {Smoking status: \code{Never},
#' \code{Former}, or \code{Current} }
#'    \item{\code{PriorSmoke}} {Smoking status: \code{1},
#' \code{2}, or \code{3} }
#' }
#' @references Data downloaded from
#' \url{http://lib.stat.cmu.edu/DASL/}.
#' @source Nierenberg, Stukel, Baron, Dain, and Greenberg, "Determinants of
#' plasma levels of beta-carotene and retinol", American Journal of
#' Epidemiology (1989).
#' @keywords datasets
#' @examples
#' 
#' data(NutritionStudy)
#' 
NA





#' 2008 Olympic Men's Marathon
#' 
#' Times for all finishers in the men's marathon at the 2008 Olympics
#' 
#' Results for all finishers in the 2008 Men's Olympic marathon in Beijing,
#' China.
#' 
#' @name OlympicMarathon
#' @docType data
#' @format A data frame with 76 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{Rank}} {Order of finish} 
#'    \item{\code{Athlete}} {Name
#' of marathoner} 
#'    \item{\code{Nationality}} {Country of marathoner}
#' 
#'    \item{\code{Time}} {Time as H:MM:SS} 
#'    \item{\code{Minutes}} {Time in minutes}
#' }
#' @source
#' \url{http://2008olympics.runnersworld.com/2008/08/mens-marathon-results.html}
#' @keywords datasets
#' @examples
#' 
#' data(OlympicMarathon)
#' 
NA





#' Ottaw Senators hockey team
#' 
#' Data for 24 players on the 2009-10 Ottawa Senators
#' 
#' Points scored and penaly minutes for 24 players (excluding goalies) playing
#' ice hockey for the Ottawa Senators during the 2009-10 NHL regular season.
#' 
#' @name OttawaSenators
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Points}} {Number of points (goals + assists) scored}
#' 
#'    \item{\code{PenMins}} {Number of penalty minutes} }
#' @source Data obtained from \url{http://senators.nhl.com/club/stats.htm}.
#' @keywords datasets
#' @examples
#' 
#' data(OttawaSenators)
#' 
NA





#' Pizza Girl Tips
#' 
#' Data on tips for pizza deliveries
#' 
#' "Pizza Girl" collected data on her deliveries and tips over three different
#' evening shifts.
#' 
#' @name PizzaGirl
#' @docType data
#' @format A dataset with 24 observations on the following 2 variables.
#' \itemize{
#'		\item{\code{Tip}} {Amount of tip (in dollars)}
#'		\item{\code{Shift}} {Which of three different shifts}
#' }
#' @source Pizza Girl: Statistical Analysis at
#' \url{http://slice.seriouseats.com/archives/2010/04/statistical-analysis-of-a-pizza-delivery-shift-20100429.html}.
#' @keywords datasets
NA





#' Quiz vs Lecture Pulse Rates
#' 
#' Paired data with pulse rates in a lecture and during a quiz for 10 students
#' 
#' Ten students in an introductory statistics class measured their pulse rate
#' in two settings: in the middle of a regular class lecture and again while
#' taking a quiz.
#' 
#' @name QuizPulse10
#' @docType data
#' @format A data frame with 10 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Student}} {ID number for the student}
#' 
#'    \item{\code{Quiz}} {Pulse rate (beats per minute) during a quiz}
#' 
#'    \item{\code{Lecture}} {Pulse rate (beats per minute) during a lecture} }
#' @source In-class data collection
#' @keywords datasets
#' @examples
#' 
#' data(QuizPulse10)
#' 
NA





#' Simulated proportions
#' 
#' Counts and proportions for 5000 simulated samples with n=200 and p=0.50
#' 
#' Results from 5000 simulations of samples of size n=200 from a population
#' with proportoin of "yes" responses at p=0.50.
#' 
#' @name RandomP50N200
#' @docType data
#' @format A data frame with 5000 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Count}} {Number of simulated "yes" responses in 200
#' trials} 
#'    \item{\code{Phat}} {Sample proportion (Count/200)} }
#' @source Computer simulation
#' @keywords datasets
#' @examples
#' 
#' data(RandomP50N200)
#' 
NA





#' Restaurant Tips
#' 
#' Tip data from the First Crush Bistro
#' 
#' The owner of a bistro called First Crush in Potsdam, NY was interested in
#' studying the tipping patterns of his customers.  He collected restaurant
#' bills over a two week period that he believes provide a good sample of his
#' customers. The data recorded from 157 bills include the amount of the bill,
#' size of the tip, percentage tip, number of customers in the group, whether
#' or not a credit card was used, day of the week, and a coded identity of the
#' server.
#' 
#' @name RestaurantTips
#' @docType data
#' @format A data frame with 157 observations on the following 7 variables.
#' \itemize{ 
#'    \item{\code{Bill}} {Size of the bill (in dollars)}
#' 
#'    \item{\code{Tip}} {Size of the tip (in dollars)} 
#'    \item{\code{CreditCard}} {Paid
#' with a credit card?  \code{No} or \code{Yes}} 
#'    \item{\code{Credit}} {Paid
#' with a credit card?  \code{n} or \code{y}} 
#'    \item{\code{Guests}} {Number of
#' peole in the group} 
#'    \item{\code{Day}} {Day of the week: \code{m}=Monday,
#' \code{t}=Tuesday, \code{w}=Wednesday, \code{th}=Thursday, or \code{f}=Friday
#' } 
#'    \item{\code{Server}} {Code for waiter/waitress: \code{A}, \code{B}, or
#' \code{C}} 
#'    \item{\code{PctTip}} {Tip as a percentage of the bill} }
#' @source Thanks to Tom DeRosa for providing the tipping data.
#' @keywords datasets
#' @examples
#' 
#' data(RestaurantTips)
#' 
NA





#' Retail Sales
#' 
#' Monthly U.S. Retail Sales (in billions)
#' 
#' Data show the monthly retail sales (in billions) for the U.S. economy in
#' each month from 2002 through 2011.
#' 
#' @name RetailSales
#' @docType data
#' @format A data frame with 144 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Month}} {Month of the year} 
#'    \item{\code{Year}} {Year (from 2002 to 2011) } 
#'    \item{\code{Date}} {Date in date format (day of month is meaningless)}
#'    \item{\code{Sales}} {U.S. retail sales (in billions of dollars)} }
#' @source \url{http://www.census.gov/retail/}
#' @keywords datasets
#' @examples
#' 
#' data(RetailSales)
#' if (require(lattice)) {
#'   xyplot(Sales ~ Date, RetailSales, type='l')
#'   xyplot(Sales ~ Date, RetailSales, type='l', groups=Month)
#' }
#' 
NA





#' Rock & Roll Hall of fame
#' 
#' Groups and Individuals in the Rock and Roll hall of Fame
#' 
#' All inductees of the Rock & Roll Hall of Fame as of 2012.
#' 
#' @name RockandRoll
#' @docType data
#' @format A data frame with 273 observations on the following 4 variables.
#' \itemize{ 
#'    \item{\code{Inductee}} {Name of he group or individual}
#' 
#'    \item{\code{FemaleMembers}} {\code{Yes} if individual or membr of the group
#' is female, otherwise \code{No}} 
#'    \item{\code{Category}} {Type of indiviudal
#' or group: \code{Performer}, \code{Non-performer}, \code{Early INfluence},
#' \code{Lifetime Achievement}, \code{Sideman}} 
#'    \item{\code{People}} {Number of
#' people in the group} }
#' @source Rock & Roll Hall of Fame website,
#' \url{http://rockhall.com/inductees/alphabetical/}
#' @keywords datasets
#' @examples
#' 
#' data(RockandRoll)
#' 
NA





#' Salary and Gender
#' 
#' Salaries for college teachers
#' 
#' A random sample of college teachers taken from the 2010 American Community
#' Survey (ACS) 1-year Public Use Microdata Sample (PUMS).
#' 
#' @name SalaryGender
#' @docType data
#' @format A dataset with 100 observations on the following 4 variables.
#' \itemize{
#' \item{\code{Salary}} {Annual salary in \$1,000's}
#' \item{\code{Gender}} {0=female or 1=male}
#' \item{\code{Sex}} {\code{Female} or \code{Male}}
#' \item{\code{Age}} {Age in years}
#' \item{\code{phD}} {\code{No} or \code{Yes}}
#' \item{\code{PhD}} {\code{0} or \code{1}}
#' }
# This original source URL is no longer valid
# @source \url{http://www.census.gov/acs/www/data_documentation/public_use_microdata_sample/}
# Here is an updated URL that seems to provide PUMS data.
#' @source \url{https://www.census.gov/main/www/pums.html}
#' @keywords datasets
NA





#' AllCountries
#' 
#' Data on a sample of countries of the world
#' 
#' A subset of data from AllCountries for a random sample of 50 countries. Data
#' for 2008 to avoid many missing values in more recent years.
#' 
#' @name SampCountries
#' @docType data
#' @format A data frame with 50 observations on the following 13 variables.
#' \itemize{ 
#'    \item{\code{Country}} {Name of the country}
#'    \item{\code{LandArea}} {Size in sq. kilometers}
#'    \item{\code{Population}} {Population in millions}
#'    \item{\code{Energy}} {Energy usage (kilotons of oil)}
#'    \item{\code{Rural}} {Percentage of population living in rural areas}
#'    \item{\code{Military}} {Percentage of government expenditures directed
#' 				toward the military} 
#'    \item{\code{Health}} {Percentage of government
#' 				expenditures directed towards healthcare} 
#'    \item{\code{HIV}} {Percentage of
#' 				the population with HIV} 
#'    \item{\code{Internet}} {Percentage of the
#' 				population with access to the internet} 
#'    \item{\code{kwhPerCap}} {An ordered factor of categories for kilowatt
#' 				hours per capita, \code{under 2500}, \code{2500 to 5000}, or \code{over 5000}}
#'    \item{\code{Developed}} {A numerical code for \code{kwhPerCap}}
#'    \item{\code{BirthRate}} {Births per 1000 people}
#'    \item{\code{ElderlyPop}} {Percentage of the population at least 65 year old}
#'    \item{\code{LifeExpectancy}} {Average life expectancy (in years)} 
#' }
#' @source Data collected from the World Bank website, worldbank.org.
#' @keywords datasets
#' @examples
#' 
#' data(SampCountries)
#' 
NA





#' S \& P 500 Prices
#' 
#' Daily data for S\& P 500 Stock Index
#' 
#' Daily prices for the S\& P 500 Stock Index for trading days in 2010.
#' 
#' @name SandP500
#' @docType data
#' @format A data frame with 252 observations on the following 6 variables.
#' \itemize{ 
#'    \item{\code{Date}} {Date as a character string} 
#'    \item{\code{MDY}} {Date as a date object} 
#'    \item{\code{Open}} {Opening value}
#' 
#'    \item{\code{High}} {High point for the day} 
#'    \item{\code{Low}} {Low point for
#' the day} 
#'    \item{\code{Close}} {Closing value} 
#'    \item{\code{Volume}} {Shares
#' traded (in millions)} }
#' @source Downladed from
#' \url{http://finance.yahoo.com/q/hp?s=^GSPC+Historical+Prices}
#' @keywords datasets
#' @examples
#' 
#' data(SandP500)
#' if (require(lattice)) {
#'   xyplot( High + Low ~ Date, data=SandP500, type="l", 
#'     main="S and P 500",
#'     auto.key=list(lines=TRUE, points=FALSE))
#'  }
#' 
NA





#' Sandwich Ants
#' 
#' Ant Counts on samples of different sandwiches
#' 
#' As young students, Dominic Kelly and his friends enjoyed watching ants
#' gather on pieces of sandwiches.  Later, as a university student, Dominic
#' decided to study this with a more formal experiment. He chose three types of
#' sandwich fillings (vegemite, peanut butter, and ham \& pickles), four types
#' of bread (multigrain, rye, white, and wholemeal), and put butter on some of
#' the sandwiches. \cr To conduct the experiment he randomly chose a sandwich,
#' broke off a piece, and left it on the ground near an ant hill.  After
#' several minutes he placed a jar over the sandwich bit and counteed the
#' number of ants. He repeated the process, allowing time for ants to return to
#' the hill after each trial, until he had two samples for each combination of
#' the factors. \cr This dataset has only sandwiches with no butter. The data
#' in SandwichAnts2 adds information for samples with butter.
#' 
#' @name SandwichAnts
#' @docType data
#' @format A data frame with 24 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{Butter}} {Butter on the sandwich? \code{no}}
#' 
#'    \item{\code{Filling}} {Type of filling: \code{Ham \& Pickles}, \code{Peanut
#' Butter}, or \code{Vegemite}} 
#'    \item{\code{Bread}} {Type of bread:
#' \code{Multigrain}, \code{Rye}, \code{White}, or \code{Wholemeal}}
#' 
#'    \item{\code{Ants}} {Number of ants on the sandwich}
#' 
#'    \item{\code{Order}} {Trial number} }
#' @source Margaret Mackisack, ``Favourite Experiments: An Addendum to What is
#' the Use of Experiments Conducted by Statistics Students?", Journal of
#' Statistics Education (1994)
#' \url{http://www.amstat.org/publications/jse/v2n1/mackisack.supp.html}
#' @keywords datasets
#' @examples
#' 
#' data(SandwichAnts)
#' 
NA





#' Sandwich Ants - Part 2
#' 
#' Ant counts on samples of different sandwiches
#' 
#' As young students, Dominic Kelly and his friends enjoyed watching ants
#' gather on pieces of sandwiches.  Later, as a university student, Dominic
#' decided to study this with a more formal experiment. He chose three types of
#' sandwich fillings (vegemite, peanut butter, and ham \& pickles), four types
#' of bread (multigrain, rye, white, and wholemeal), and put butter on some of
#' the sandwiches. \cr To conduct the experiment he randomly chose a sandwich,
#' broke off a piece, and left it on the ground near an ant hill.  After
#' several minutes he placed a jar over the sandwich bit and counteed the
#' number of ants. He repeated the process, allowing time for ants to return to
#' the hill after each trial, until he had two samples for each combination of
#' the three factors.
#' 
#' @name SandwichAnts2
#' @docType data
#' @format A data frame with 48 observations on the following 5 variables.
#' \itemize{ 
#'    \item{\code{Butter}} {Butter on the sandwich? \code{no}}
#' 
#'    \item{\code{Filling}} {Type of filling: \code{Ham \& Pickles}, \code{Peanut
#' Butter}, or \code{Vegemite}} 
#'    \item{\code{Bread}} {Type of bread:
#' \code{Multigrain}, \code{Rye}, \code{White}, or \code{Wholemeal}}
#' 
#'    \item{\code{Ants}} {Number of ants on the sandwich}
#' 
#'    \item{\code{Order}} {Trial number} }
#' @source Margaret Mackisack, ``Favourite Experiments: An Addendum to What is
#' the Use of Experiments Conducted by Statistics Students?", Journal of
#' Statistics Education (1994)
#' \url{http://www.amstat.org/publications/jse/v2n1/mackisack.supp.html}
#' @keywords datasets
#' @examples
#' 
#' data(SandwichAnts2)
#' 
NA





#' Skateboard Prices
#' 
#' Prices of skateboards for sale online
#' 
#' Prices for skateboards offered for sale on eBay.
#' 
#' @name SkateboardPrices
#' @docType data
#' @format A dataset with 20 observations on the following variable.
#' \itemize{
#'		\item{\code{Price}} {Selling price in dollars}
#' }
#' @source Random sample taken from all skateboards available for sale on eBay
#' on February 12, 2012.
#' @keywords datasets
NA





#' Sleep Caffeine
#' 
#' Experimentn to compare word recall after sleep or caffeine
#' 
#' A random sample of 24 adults were divided equally into two groups and given
#' a list of 24 words to memorize. During a break, one group takes a 90 minute
#' nap while another group is given a caffeine pill.  The response variable is
#' the number of words participants are able to recall following the break.
#' 
#' @name SleepCaffeine
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Group}} {Treatment: \code{Caffeine} \code{Sleep}}
#'    \item{\code{Words}} {Number of words recalled} 
#' }
#' @source Mednick, Cai, Kanady, and Drummond, "Comparing the benefits of
#' caffeine, naps and placebo on verbal, motor and perceptual memory",
#' Behavioural Brain Research, 193 (2008), 79-86.
#' @keywords datasets
#' @examples
#' 
#' data(SleepCaffeine)
#' 
NA





#' Sleep Study
#' 
#' Data from a study of sleep patterns for college students.
#' 
#' The data were obtained from a sample of students who did skills tests to
#' measure cognitive function, completed a survey that asked many questions
#' about attitudes and habits, and kept a sleep diary to record time and
#' quality of sleep over a two week period.
#' 
#' @name SleepStudy
#' @docType data
#' @format A dataset with 253 observations on the following 27 variables.
#' \itemize{
#' \item{\code{Gender}} {1=male, 0=female}
#' \item{\code{Sex}} {\code{Female} or \code{Male}}
#' \item{\code{ClassYear}} {Year in school, 1=first year, ..., 4=senior}
#' \item{\code{LarkOwl}} {Early riser or night owl?
#' 		\code{Lark}, \code{Neither}, or \code{Owl}}
#' \item{\code{NumEarlyClass}} {Number of classes per week before 9 am}
#' \item{\code{earlyClass}} {Indicator for any early classes}
#' \item{\code{EarlyClass}} {Indicator for any early classes}
#' \item{\code{GPA}} {Grade point average (0-4 scale)}
#' \item{\code{ClassesMissed}} {Number of classes missed in a semester}
#' \item{\code{CognitionZscore}} {Z-score on a test of cognitive skills}
#' \item{\code{PoorSleepQuality}} {Measure of sleep quality (higher values are poorer sleep)}
#' \item{\code{DepressionScore}} {Measure of degree of depression}
#' \item{\code{AnxietyScore}} {Measure of amount of anxiety}
#' \item{\code{StressScore}} {Measure of amount of stress}
#' \item{\code{DepressionStatus}} {Coded depression score: \code{normal},
#' 		\code{moderate}, or \code{severe}}
#' \item{\code{AnxietyStatus}} {Coded anxiety score: 
#'      \code{normal}, \code{moderate}, or \code{severe}}
#' \item{\code{Stress}} {Coded stress score: \code{normal} or \code{high}}
#' \item{\code{DASScore}} {Combined score for depression, anxiety and stress}
#' \item{\code{Happiness}} {Measure of degree of happiness}
#' \item{\code{AlcoholUse}} {Self-reported:
#' 		\code{Abstain}, \code{Light}, \code{Moderate}, or \code{Heavy}}
#' \item{\code{Drinks}} {Number of alcoholic drinks per week}
#' \item{\code{WeekdayBed}} {Average weekday bedtime (24.0=midnight)}
#' \item{\code{WeekdayRise}} {Average weekday rise time (8.0=8 am)}
#' \item{\code{WeekdaySleep}} {Average hours of sleep on weekdays}
#' \item{\code{WeekendBed}} {Average weekend bedtime (24.0=midnight)}
#' \item{\code{WeekendRise}} {Average weekend rise time (8.0=8 am)}
#' \item{\code{WeekendSleep}} {Average weekend bedtime (24.0=midnight)}
#' \item{\code{AverageSleep}} {Average hours of sleep for all days}
#' \item{\code{allNighter}} {Had an all-nighter this semester?  \code{Yes} or \code{No}}
#' \item{\code{AllNighter}} {Had an all-nighter this semester?  \code{0} or \code{1}}
#' }
#' @source Onyper, S., Thacher, P., Gilbert, J., Gradess, S., "Class Start
#' Times, Sleep, and Academic Performance in College: A Path Analysis," April
#' 2012; 29(3): 318-335.  Thanks to the authors for supplying the data.
#' @keywords datasets
NA





#' Smiles
#' 
#' Experiment to study effect of smiling on leniency in judicial matters
#' 
#' Hecht and LeFrance conducted a study examining the effect of a smile on the
#' leniency of disciplinary action for wrongdoers. Participants in the
#' experiment took on the role of members of a college disciplinary panel
#' judging students accused of cheating.  For each suspect, along with a
#' description of the offense, a picture was provided with either a smile or
#' neutral facial expression.  A leniency score was calculated based on the
#' disciplinary decisions made by the participants.
#' 
#' @name Smiles
#' @docType data
#' @format A data frame with 68 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Leniency}} {Score assigned by a judgement panel
#' (higher is more lenient)} 
#'    \item{\code{Group}} {Treatment group:
#' \code{neutral} or \code{smile}} }
#' @source LaFrance, M., & Hecht, M. A., "Why smiles generate leniency",
#' Personality and Social Psychology Bulletin, 21, 1995, 207-214.
#' @keywords datasets
#' @examples
#' 
#' data(Smiles)
#' 
NA





#' Speed Dating
#' 
#' Data from a sample of four minute speed dates.
#' 
#' Participants were students at Columbia's graduate and professional schools,
#' recruited by mass email, posted fliers, and fliers handed out by research
#' assistants. Each participant attended one speed dating session, in which
#' they met with each participant of the opposite sex for four minutes. Order
#' and session assignments were randomly determined. After each four minute
#' "speed date," participants filled out a form rating their date on a scale of
#' 1-10 on various attributes.  Only data from the first date in each session
#' is recorded here.
#' 
#' @name SpeedDating
#' @docType data
#' @format A dataset with 276 observations on the following 22 variables.
#' \itemize{ 
#' \item{\code{DecisionMale}} {Would the male like another date?  \code{Yes} or \code{No}}
#' \item{\code{DecisionM}} {Would the male like another date?  \code{0} or \code{1}}
#'    \item{ \code{DecisionFemale}} {Would the female like another date? \code{Yes} or \code{No}}
#'    \item{ \code{DecisionF}} {Would the female like another date? \code{0} or \code{1}}
#'    \item{ \code{LikeM}} {How much the male likes his partner (1-10 scale)}
#'    \item{ \code{LikeF}} {How much the female likes her partner (1-10 scale)}
#'    \item{ \code{PartnerYesM}} {Male's estimate of chance the female wants
#' 		another date (1-10 scale)}
#'    \item{ \code{PartnerYesF}} {Female's estimate of chance the male wants another 
#'        date (1-10 scale)}
#'    \item{ \code{AgeM}} {Male's age (in years)}
#'    \item{ \code{AgeF}} {Females age (in years)}
#'    \item{ \code{RaceM}} {Male's race: 
#' 		\code{Asian}, \code{Black}, \code{Caucasian}, \code{Latino}, or \code{Other}}
#'    \item{ \code{RaceF}} {Female's race: 
#' 		\code{Asian}, \code{Black}, \code{Caucasian}, \code{Latino}, or \code{Other}}
#'    \item{\code{AttractiveM}} {Male's rating of female's attractiveness (1-10 scale)}
#'    \item{ \code{AttractiveF}} {Female's rating of male's attractiveness (1-10 scale)}
#'    \item{ \code{SincereM}} {Male's rating of female's sincerity (1-10 scale)}
#'    \item{ \code{SincereF}} {Female's rating of male's sincerity (1-10 scale)}
#'    \item{ \code{IntelligentM}} {Male's rating of female's intelligence (1-10 scale)}
#'    \item{ \code{IntelligentF}} {Female's rating of male's intelligence (1-10 scale)}
#'    \item{ \code{FunM}} {Male's rating of female as fun (1-10 scale)}
#'    \item{ \code{FunF}} {Female's rating of male as fun (1-10 scale)}
#'    \item{ \code{AmbitiousM}} {Male's rating of female's ambition (1-10 scale)}
#'    \item{ \code{AmbitiousF}} {Female's rating of male's ambition (1-10 scale)}
#'    \item{ \code{SharedInterestsM}} {Male's rating of female's shared interests 
#'             (1-10 scale)}
#'    \item{ \code{SharedInterestsF}} {Female's rating of male's shared interests 
#'             (1-10 scale)}
#' }
#' @source Gelman, A. and Hill, J., Data analysis using regression and
#' multilevel/hierarchical models, Cambridge University Press: New York, 2007
#' @keywords datasets
NA





#' Statistics Exam Grades
#' 
#' Grades on statistics exams
#' 
#' Exam scores for a sample of students who completed a course using
#' Statistics: Unlocking the Power of Data as a text.  The dataset contains
#' scores on Exam 1 (Chapters 1 to 4), Exam 2 (Chapters 5 to 8), and the Final
#' exam (entire book).
#' 
#' @name StatGrades
#' @docType data
#' @format A dataset with 50 observations on the following 3 variables.
#' \itemize{
#'		\item{\code{Exam1}} {Score (out of 100 points) on the first exam}
#'		\item{\code{Exam2}} {Score (out of 100 points) on the first exam}
#' 		\item{\code{Final}} {Score (out of 100 points) on the final exam}
#' }
#' @source Random selection of students in an introductory statistics course.
#' @keywords datasets
NA





#' Statistics PhD Programs
#' 
#' Enrollments in Statistics PhD Programs
#' 
#' Graduate student enrollments in Statistics and Biostatistics departments in
#' 2009. The list does not include combined departments of mathematics and
#' statistics and does not include departments that did not reply to the AMS
#' survey.
#' 
#' @name StatisticsPhD
#' @docType data
#' @format A data frame with 82 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{University}} {Name of the school}
#' 
#'    \item{\code{Department}} {Type of department: \code{Biostatistics} or
#' \code{Statistics}} 
#'    \item{\code{FTGradEnrollment}} {Full time graduate
#' student enrollment} }
#' @source The full list of the 82 Group IV departments was obtained at
#' \url{http://www.ams.org/profession/data/annual-survey/group_iv}. Data on
#' enrollment were obtained primarily from Assistantships and Graduate
#' Fellowships in the Mathematical Sciences, 2009, American Mathematical
#' Society.
#' @keywords datasets
#' @examples
#' 
#' data(StatisticsPhD)
#' 
NA





#' Stock Changes
#' 
#' Stock price change for a smaple of stocks from the S \& P 500 (August 2-6,
#' 2010)
#' 
#' A random sample of 50 companies from Standard \& Poor's index of 500
#' companies was selected. The change in the price of the stock (in dollars)
#' over the 5-day period from August 2 - 6, 2010 was recorded for each company
#' in the sample.
#' 
#' @name StockChanges
#' @docType data
#' @format A data frame with 50 observations on the following variable.
#' \itemize{ 
#'    \item{\code{SPChange}} {Change in sock price (in dollars)} }
#' @keywords datasets
#' @examples
#' 
#' data(StockChanges)
#' 
NA





#' Story Spoilers
#' 
#' Raitngs for stories with and without spoilers
#' 
#' This study investigated whether a story spoiler that gives away the ending
#' early diminishes suspense and hurts enjoyment?  For twelve different short
#' stories, the study's authors created a second version in which a spoiler
#' paragraph at the beginning discussed the story and revealed the outcome.
#' Each version of the twelve stories was read by at least 30 people and rated
#' on a 1 to 10 scale to create an overall rating for the story, with higher
#' ratings indicating greater enjoyment of the story.  Stories 1 to 4 were
#' ironic twist stories, stories 5 to 8 were mysteries, and stories 9 to 12
#' were literary stories.
#' 
#' @name StorySpoilers
#' @docType data
#' @format A dataset with 12 observations on the following 3 variables.
#' \itemize{
#'		\item{\code{Story}} {ID for story}
#'		\item{\code{Spoiler}} {Average (0-10) rating for spoiler version}
#'		\item{\code{Original}} {Average (0-10) rating for original version}
#' }
#' @source Leavitt, J. and Christenfeld, N., "Story Spoilers Don't Spoil
#' Stories," Psychological Science, published OnlineFirst, August 12, 2011.
#' @keywords datasets
NA





#' Stressed Mice
#' 
#' Time in darkness for mice in different environments
#' 
#' In the recent study, mice were randomly assigned to either an enriched
#' environment where there was an exercise wheel available, or a standard
#' environment with no exercise options.  After three weeks in the specified
#' environment, for five minutes a day for two weeks, the mice were each
#' exposed to a "mouse bully" - a mouse who was very strong, aggressive, and
#' territorial.  one measure of mouse anxiety is amount of time hiding in a
#' dark compartment, with mice who are more anxious spending more time in
#' darkness. The amount of time spent in darkness is recorded for each of the
#' mice.
#' 
#' @name StressedMice
#' @docType data
#' @format A data frame with 14 observations on the following 2 variables.
#' \itemize{ 
#'    \item{\code{Time}} {Time spent in darkness (in seconds)}
#' 
#'    \item{\code{Environment}} {Type of environment: \code{Enriched} or
#' \code{Standard}} }
#' @source Data approximated from summary statistics in: Lehmann and Herkenham,
#' "Environmental Enrichment Confers Stress Resiliency to Social Defeat through
#' an Infralimbic Cortex-Dependent Neuroanatomical Pathway", The Journal of
#' Neuroscience, April 20, 2011, 31(16):61596173.
#' @keywords datasets
#' @examples
#' 
#' data(StressedMice)
#' 
NA





#' Student Survey Data
#' 
#' Data from a survey of students in introductory statistics courses
#' 
#' Data from an in-class survey given to introductory statistics students over
#' several years.
#' 
#' @name StudentSurvey
#' @docType data
#' @format A data frame with 362 observations on the following 17 variables.
#' \itemize{ 
#'    \item{\code{Year}} {Year in school}
#' 
#'    \item{\code{Gender}} {Student's gender: \code{F} or \code{M}}
#' 
#'    \item{\code{Smoke}} {Smoker? \code{No} or \code{Yes}}
#' 
#'    \item{\code{Award}} {Prefered award: \code{Academy} \code{Nobel}
#' \code{Olympic}} 
#'    \item{\code{HigherSAT}} {Which SAT is higher?  \code{Math}
#' or \code{Verbal}} 
#'    \item{\code{Exercise}} {Hours of exercsie per week}
#' 
#'    \item{\code{TV}} {Hours of TV viewing per week} 
#'    \item{\code{Height}} {Height
#' (in inches)} 
#'    \item{\code{Weight}} {Weight (in pounds)}
#' 
#'    \item{\code{Siblings}} {Number of siblings} 
#'    \item{\code{BirthOrder}} {Birth
#' order, 1=oldest} 
#'    \item{\code{VerbalSAT}} {Verbal SAT score}
#' 
#'    \item{\code{MathSAT}} {Math SAT score} 
#'    \item{\code{SAT}} {Combined Verbal +
#' Math SAT} 
#'    \item{\code{GPA}} {College grade point average}
#' 
#'    \item{\code{Pulse}} {Pulse rate (beats per minute)}
#' 
#'    \item{\code{Piercings}} {Number of body piercings} }
#' @source In-class student survey
#' @keywords datasets
#' @examples
#' 
#' data(StudentSurvey)
#' 
NA





#' Ten Countries
#' 
#' A subset of the ALLCountries data for a random sample of ten countries
#' 
#' 
#' @name TenCountries
#' @docType data
#' @format A data frame with 10 observations on the following 4 variables.
#' \itemize{ 
#'    \item{\code{Country}} {a factor with levels \code{Armenia}
#' \code{Bahamas} \code{Lebanon} \code{Macedonia} \code{North Korea}
#' \code{Romania} \code{Serbia} \code{Slovenia} \code{Tunisia}
#' \code{Uzbekistan}} 
#'    \item{\code{Code}} {a factor with levels \code{ARM}
#' \code{BHS} \code{LBN} \code{MKD} \code{PRK} \code{ROU} \code{SRB} \code{SVN}
#' \code{TUN} \code{UZB}} 
#'    \item{\code{Area}} {a numeric vector}
#' 
#'    \item{\code{PctRural}} {a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(TenCountries)
#' ## maybe str(TenCountries) ; plot(TenCountries) ...
#' 
NA





#' Textbook Costs
#' 
#' Prices for textbooks for different courses
#' 
#' Data are from samples of ten courses in each of four disciplines at a
#' liberal arts college.  For each course the bookstore's website lists the
#' required texts(s) and costs. Data were collected for the Fall 2011 semester.
#' 
#' @name TextbookCosts
#' @docType data
#' @format A data frame with 40 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Field}} {General discipline of the course:
#' \code{Arts}, \code{Humanities}, \code{NaturalScience}, or
#' \code{SocialScience}} 
#'    \item{\code{Books}} {Number of books requiired}
#' 
#'    \item{\code{Cost}} {Total cost (in dollars) for required books} }
#' @source Bookstore online site
#' @keywords datasets
#' @examples
#' 
#' data(TextbookCosts)
#' 
NA





#' Toenail Arsenic
#' 
#' Arsenic in toenails of 19 people using private wells in New Hampshire
#' 
#' Level of arsenic was measured in toenails of 19 subjects from New Hampshire,
#' all with private wells as their main water source.
#' 
#' @name ToenailArsenic
#' @docType data
#' @format A data frame with 19 observations on the following variable.
#' \itemize{ 
#'    \item{\code{Arsenic}} {Level of arsenic found in toenails} }
#' @source Adapted from Karagas, et.al.,"Toenail Samples as an Indicator of
#' Drinking Water Arsenic Exposure", Cancer Epidemiology, Biomarkers and
#' Prevention 1996;5:849-852.
#' @keywords datasets
#' @examples
#' 
#' data(ToenailArsenic)
#' 
NA





#' Traffic Flow
#' 
#' Traffic flow times from a simulation with timed and flexible traffic lights
#' 
#' Engineers in Dresden, Germany were looking at ways to improve traffic flow
#' by enabling traffic lights to communicate information about traffic flow
#' with nearby traffic lights. The data in show results of one experiment where
#' they simulated buses moving along a street and recorded the delay time (in
#' seconds) for both a fixed time and a flexible system of lights. The process
#' was repeated under both conditions for a sample of 24 simulated scenarios.
#' 
#' @name TrafficFlow
#' @docType data
#' @format A data frame with 24 observations on the following 3 variables.
#' \itemize{ 
#'    \item{\code{Timed}} {Delay time (in minutes) for fixed timed
#' lights} 
#'    \item{\code{Flexible}} {Delay time (in minutes) for flexible
#' communicating lights.} 
#'    \item{\code{Difference}} {Difference (Timed-Flexible)
#' for each simulation} }
#' @source Lammer and Helbing, "Self-Stabilizing decentralized signal control
#' of realistic, saturated network traffic", Santa Fe Institute working paper
#' \# 10-09-019, September 2010.
#' @keywords datasets
#' @examples
#' 
#' data(TrafficFlow)
#' 
NA





#' US State Data
#' 
#' Various data for all 50 US States
#' 
#' 
#' @name USStates
#' @docType data
#' @format A data frame with 50 observations on the following 17 variables.
#' \itemize{ 
#'    \item{\code{State}} {Name of state}
#' 
#'    \item{\code{HouseholdIncome}} {Mean household income (in dollars)}
#' 
#'    \item{\code{IQ}} {Mean IQ score of residents}
#' 
#'    \item{\code{McCainVote}} {Percentage of votes for John McCain in 2008
#' Presidential election} 
#'    \item{\code{Region}} {Area of the country:
#' \code{MW}=Midwest, \code{NE}=Northeast, \code{S}=South, or \code{W}=West}
#' 
#'    \item{\code{ObamaMcCain}} {Which 2008 Presidential candidate won state?
#' \code{M}=McCain or \code{O}=Obama} 
#'    \item{\code{Pres2008}} {Which 2008 Presidential candidate won state?
#' \code{M}=McCain or \code{O}=Obama} 
#'    \item{\code{Population}} {Number of
#' residents (in millions)} 
#'    \item{\code{EighthGradeMath}} {a numeric vector}
#' 
#'    \item{\code{HighSchool}} {Percentage of high school graduates}
#' 
#'    \item{\code{GSP}} {Gross State Product (dollars per capita)}
#' 
#'    \item{\code{FiveVegetables}} {Percentage of residents who eat at least five
#' servings of fruits/vegetables per day} 
#'    \item{\code{Smokers}} {Percentage of
#' residents who smoke} 
#'    \item{\code{PhysicalActivity}} {Percentage of residents
#' who have competed in a physical activity in past month}
#' 
#'    \item{\code{Obese}} {Percentage of residents classified as obese}
#' 
#'    \item{\code{College}} {Percentage of residents with college degrees}
#' 
#'    \item{\code{NonWhite}} {Percentage of residents who are not white}
#' 
#'    \item{\code{HeavyDrinkers}} {Percentage of residents who drink heavily} }
#' @source Various online sources, mostly at www.census.gov
#' @keywords datasets
#' @examples
#' 
#' data(USStates)
#' 
NA





#' Water Striders
#' 
#' Mating activity for water striders
#' 
#' Water striders are common bugs that skate across the surface of water. Water
#' striders have different personalities and some of the males are
#' hyper-aggressive, meaning they jump on and wrestle with any other water
#' strider near them. Individually, because hyper-aggressive males are much
#' more active, they tend to have better mating success than more inactive
#' striders. This study examined the effect they have on a group. Four males
#' and three females were put in each of ten pools of water. Half of the groups
#' had a hyper-aggressive male as one of the males and half did not. The
#' proportion of time females are in hiding was measured for each of the 10
#' groups, and a measure of mean mating activity was also measured with higher
#' numbers meaning more mating.
#' 
#' @name WaterStriders
#' @docType data
#' @format A dataset with 10 observations on the following 3 variables.
#' \itemize{
#'		\item{\code{AggressiveMale}} {Hyper-aggressive male in group?
#' 			\code{No} or \code{Yes}}
#'		\item{\code{FemalesHiding}}  {Proportion of time the
#' 				female water striders were in hiding}
#'		\item{\code{MatingActivity}} {Measure of mean mating activity 
#'				(higher numbers meaning more mating)}
#' }
#' @source Sih, A. and Watters, J., "The mix matters: behavioural types and
#' group dynamics in water striders," Behaviour, 2005; 142(9-10): 1423.
#' @keywords datasets
NA





#' WaterTaste
#' 
#' Blind taste test to compare brands of bottled water
#' 
#' Result from a blind taste test comparing different four different types of
#' water (Sam's Choice, Aqufina, Fiji, and tap water).  Participants rank
#' ordered waters when presented in a random order.
#' 
#' @name WaterTaste
#' @docType data
#' @format A data frame with 100 observations on the following 10 variables.
#' \itemize{ 
#'    \item{\code{Gender}} {Gender of respondent: \code{F}=Female
#' \code{M}=Male} 
#'    \item{\code{Age}} {Age (in years)} 
#'    \item{\code{Class}} {Year
#' in school \code{F}=First year \code{J}=Junior \code{O}=Other \code{P}
#' \code{SO}=Sophomore \code{SR}=Senior} 
#'    \item{\code{UsuallyDrink}} {Usual
#' source of drinking water: \code{Bottled}, \code{Filtered}, or \code{Tap}}
#' 
#'    \item{\code{FavBotWatBrand}} {Favorite brand of bottled water }
#' 
#'    \item{\code{Preference}} {Order of perference: \code{A}=Sams Choice,
#' \code{B}=Aquafina, \code{C}=Fiji, and \code{D}=Tap water }
#' 
#'    \item{\code{First}} {Top choice among \code{Aquafina}, \code{Fiji},
#' \code{SamsChoice}, or \code{Tap}} 
#'    \item{\code{Second}} {Second choice}
#' 
#'    \item{\code{Third}} {Third choice } 
#'    \item{\code{Fourth}} {Fourth choice} }
#' @references \url{http://www.amstat.org/publications/jse/v18n1/lunsford.pdf}
#' @source "Water Taste Test Data" by M. Leigh Lunsford and Alix D. Dowling
#' Finch in the Journal of Statistics Education (Vol 18, No, 1) 2010
#' @keywords datasets
#' @examples
#' 
#' data(WaterTaste)
#' 
NA





#' Wetsuits
#' 
#' Swim velocity (for 1500 meters) with and withut wearing a wetsuit
#' 
#' A study tested whether wearing wetsuits influences swimming velocity.
#' Twelve competitive swimmers and triathletes swam 1500m at maximum speed
#' twice each; once wearing a wetsuit and once wearing a regular bathing suit.
#' The order of the trials was randomized.  Each time, the maximum velocity in
#' meters/sec of the swimmer was recorded.
#' 
#' @name Wetsuits
#' @docType data
#' @format A data frame with 12 observations on the following 4 variables.
#' \itemize{ 
#'    \item{\code{Wetsuit}} {Maximum swim velocity (m/sec) when wearing
#' a wetsuit} 
#'    \item{\code{NoWetsuit}} {Maximum swim velocity (m/sec) when
#' wearing a regular bathing suit} 
#'    \item{\code{Gender}} {Gender of swimmer:
#' \code{F} or \code{M}} 
#'    \item{\code{Type}} {Type of athlete: \code{swimmer} or
#' \code{triathlete}} }
#' @source de Lucas, R.D., Balildan, P., Neiva, C.M., Greco, C.C., Denadai,
#' B.S. (2000). "The effects of wetsuits on physiological and biomechanical
#' indices during swimming," Journal of Science and Medicine in Sport, 3 (1):
#' 1-8.
#' @keywords datasets
#' @examples
#' 
#' data(Wetsuits)
#' 
NA



