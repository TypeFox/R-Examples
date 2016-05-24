#' Alcohol at Georgetown College
#' 
#' Alcohol policy violations on the Georgetown College campus over several years.
#' 
#' 
#' @name alcohol
#' @docType data
#' @format A data frame with 10 observations on the following 4 variables.
#' \describe{ \item{year}{Academic year ending with Spring of the given
#' year.} 
#' \item{enrollment}{Full-time equivalent enrollment.}
#' \item{writeups}{Number of write-ups for alcohol violations.}
#' \item{writeups.per.100}{Number of writeups per 100 students.} }
#' @source Collected by MAT 111 students as a project.
#' @keywords datasets
#' 

NA


#' 
#' Attitudes Experiment 2001
#' 
#' Study conducted in November 2001 by students in MAT 111.  Subjects were 267
#' Georgetown College students.  Not all subjects got the same survey form.
#' 
#' @name attitudes
#' @docType data
#' @format A data frame with 268 observations on the following 8 variables.
#' \describe{ 
#' \item{def.race}{Suggested race of the defendant in the survey form.}
#' \item{vic.race}{Suggested race of the victim in the survey form.}
#' \item{conc.situation}{Scenario described in the
#' in the "rock concert" question on the survey form.}
#' \item{sentence}{Sentence, in years, reccommended for the defendant.}
#' \item{conc.decision}{Whether or not the subject chose to buy a ticket (or buy another
#' ticket).}
#' \item{year}{Class rank of the subject.}
#' \item{sex}{a factor with level Sex of the survey participant.}
#' \item{major}{possible values: \code{humanities} \code{math.sci} \code{pre.prof}
#' \code{social.sci} Type of major the subject intends.}
#' }
#' 
#' @details
#' Here is a sample survey form, with variants noted.
#' 
#' Attitudes Survey
#' 
#' Crime: You are on a jury for a manslaughter case in Lewistown, PA.  The
#' defendant has been found guilty, and in Pennsylvania it is part of the job
#' of the jury to recommend a sentence to the judge.  The facts of the case are
#' as follows.  The defendant, Tyrone Marcus Watson, a 35-year old native of
#' Lewistown, was driving under the influence of alchohol on the evening of
#' Tuesday July 17, 2001.  At approximately 11:00 PM Watson drove through a red
#' light, striking a pedestrian, Betsy Brockenheimer, a 20-year old resident of
#' Lewistown.  Brockenheimer was taken unconscious to the hospital and died of
#' her injuries about one hour later.  Watson did not flee the scene, nor did
#' he resist arrest.
#' 
#' The prior police record for Mr. Watson is as follows: two minor traffic
#' violations, and one previous arrest, five years ago, for DUI.  No one was
#' hurt in that incident.
#' 
#' Watson has now been convicted of DUI and manslaughter.  The minimum jail
#' term for this combination of offenses is two years; the maximum term is
#' fifty years.  In the blank below, write a number from 2 to 50 as your
#' recommended length of sentence for Tyrone Marcus Watson.  _______________
#' 
#' [In the question above, name of defedant could vary: either William Shane
#' Winchester or Tyrone Marcus Watson.  The name of the victim could also vary:
#' either Betsy Brockenheimer or Latisha Dawes.]
#' 
#' Spending Habits
#' 
#' You have purchased a $30 ticket to see a rock concert in Rupp Arena.  When
#' you arrive at the Arena on the night of the performance, you find that you
#' have lost the ticket.  You have no receipt, so it will not be possible to
#' see the concert unless you purchase another ticket.  Would you purchase
#' another ticket?  Circle below.
#' 
#' YES NO
#' 
#' [In other forms, the question above could have been: You plan to see a rock
#' concert in Rupp Arena.  Tickets for the performance are $30. When you arrive
#' at the Arena on the night of the performance, you find that you have lost
#' two bills from your purse or wallet: a ten and a twenty.  Would you buy the
#' ticket anyway?]
#' 
#' Respondent Data
#' 
#' I am (circle one): freshman sophomore junior senior
#' 
#' I am (circle one) male female
#' 
#' (Optional) My intended major is: _____________________
#' 
#' @source Georgetown College
#' @keywords datasets


NA


#' 
#' Beans
#' 
#' Experiment performed at UC-Davis; fifteen students participated.  Each
#' student was asked to place as many beans into a cup as he/she could, in 15
#' seconds.  Each student performed this task once with the dominant hand, and
#' once with the nondominat hand, but the order of performance was randomized.
#' The purpose of the study was to see whether manual dexterity was better for
#' the dominant hand.  Terminology: your dominant hand is the hand you use the
#' most.
#' 
#' 
#' @name beans
#' @docType data
#' @format A data frame with 15 observations on the following 3 variables.
#' \describe{ \item{Dom}{Number of beans placed into
#' cup with the dominant hand.}
#' \item{NonDom}{Number
#' of beans placed with the nondominant hand.}
#' \item{Diff}{Difference in number of beans placed (dominant hand minus nondominant
#' hand).} }
#' @source Uts and Heckard, Mind on Statistics, 4th Edition.
#' @keywords datasets


NA

#' Miguel Cabrera
#' 
#' PITCHf/x data on nine-time All-Star Miguel Cabrera, who is thought to be one of the best
#' pure hitters in baseball.  Covers seasons 2009 through 2012.
#' 
#' 
#' @name cabrera
#' @docType data
#' @format A data frame with 6265 observations on the following 12 variables.  Each observation
#' is a pitch to Cabrera.
#' \describe{ \item{season}{The year of play}
#' \item{gamedate}{Date of the game in which the pitch was thrown}
#' \item{pitch_type}{Type of pitch thrown, as determined by a computer algorithm.}
#' \item{balls}{Current ball count}
#' \item{strikes}{Current strike count}
#' \item{speed}{speed of pitch (in mph).  (When crossing plate?)}
#' \item{px}{x-coordinate of pitch (in feet, measured from center of plate)}
#' \item{pz}{vertical coordinate of pitch (in feet above plate)}
#' \item{swung}{Whether or not cabrera swung at the ball.  Factor with levels "no", "yes".}
#' \item{hitx}{x-coordinate of landing point of ball (if it was hit).  Relative to park.}
#' \item{hity}{y-coordinate of landing point of ball (if it was hit).  Relative to park.}
#' \item{hit_outcome}{Outcome when ball was hit.  Factor with levels E (error), H (hit), O (batter out).}
#' }
#' @source Marchi and Albert:  Analyziing Baseball Data with R, CRC Press 2014.  For
#' more on the PITCHf/x system, see \url{http://en.wikipedia.org/wiki/PITCHf/x}.
#' 
#' @keywords datasets


NA

#' 
#' Time to Chug
#' 
#' College-aged males chugging a 12-ounce can of a certain beverage.
#' 
#' 
#' @name chugtime
#' @docType data
#' @format A data frame with 13 observations on the following 2 variables.
#' \describe{
#' \item{Weight}{Weight of the subject (in pounds).}
#' \item{ChugTime}{How long (in seconds) the subject requires to drink the
#' beverage.}
#' }
#' @source Utts and Heckard, Mind on Statistics, 4th Edition.
#' @keywords datasets
#' 


NA

#' The Death Penalty and Race
#' 
#' A dataset recreated from summary data that describes relationships 
#' between race of defendant, race of victim, and outcome of trial in a number of capital 
#' cases in Florida in 1976-1977.

#' @name deathpen
#' @docType data
#' @format A data frame with 326 rows and 3 variables
#' 
#' \describe{
#'   \item{defrace}{Race of the defendant in the capital case}  
#'   \item{vicrace}{Race of the victim}
#'   \item{death}{Whether or not the defendant in the case received the death penalty}
#' }
#' 
#' @source Michael J. Radelet:  "Racial Characteristics and 
#' the Imposition of the Death Penalty", \emph{American Sociological Review}, 46 (1981).
#' @keywords datasets


NA


#' Speed and Fuel Efficiency (British Ford Escort)
#' 
#' A British Ford Escort was driven along a prescribed course.  Each drive was
#' done at a different speed, and the fuel efficiency was recorded for each
#' drive.
#' 
#' 
#' @name fuel
#' @docType data
#' @format A data frame with 15 observations on the following 2 variables.
#' \describe{
#' \item{speed}{in kilometers per hour.}
#' \item{efficiency}{fuel efficiency, measured in liters of fuel required to travel 100 kilometers. }
#' }
#' @source The Basic Practice of Statistics, by Moore and McCabe.
#' @keywords datasets

NA

#' Galton's Father-Son Data
#' 
#' Data on father-son pairs.  Collected in 1885 by Francis Galton.
#' 
#' 
#' @name galton
#' @docType data
#' @format A data frame with 1078 observations on the following 2 variables.
#' \describe{
#' \item{fheight}{Height of the father, in inches.}
#' \item{sheight}{Height of the son, in inches.}
#' }
#' @keywords datasets


NA

#' Feelings About Georgetown College
#' 
#' Results of a survey conducted by Georgetown College students on 47 Georgetown College upperclass students.
#'
#' 
#' @name gcfeeling
#' @docType data
#' @format A data frame with 47 observations on the following 6 variables.
#' \describe{
#' \item{rating.fresh}{how happy the subjects remembers being as a first-year student, on a scale of 1 to 10.}
#' \item{rating.js}{how happy the subjects feels now, on a scale of 1 to 10.}
#' \item{greek}{whether or not the subject belongs to a greek organization.}
#' \item{athlete}{whether or not the subject is a varsity athlete}
#' \item{rating.diff}{upper-level happiness rating minus remembered first-year rating}
#' \item{happier}{whether or not subject feels happier now than as a first-year student}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets

NA

#' Georgetown College Students
#' 
#' Data collected by GC students.
#'
#' 
#' @name gcstudents
#' @docType data
#' @format A data frame with 62 observations on the following 4 variables.
#' \describe{ \item{height}{height of the survey participant, in inches}
#' \item{GPA}{grade-point average}
#' \item{enough_Sleep}{Does the particpant feel that he/she gets enough sleep?}
#' \item{sex}{sex of the survey participant}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets

NA

#' General Social Survey, 2002
#' 
#' The General Social Survey (GSS) is a nationwide poll that has been conducted
#' since 1972 (semiannually since 1994).  Most interviews are done
#' face-to-face.  For further information, see below.
#'
#' 
#' @name gss02
#' @docType data
#' @format A data frame with 2765 observations on the following 13 variables.
#' \describe{
#' \item{sex}{a factor with levels \code{Female}
#' \code{Male}}
#' \item{race}{a factor with levels \code{AfrAm}
#' \code{Hispanic} \code{Other} \code{White}}
#' \item{degree}{a factor
#' with levels \code{Bachelor} \code{Graduate} \code{HighSchool} \code{JunColl}
#' \code{NotHs}}
#' \item{relig}{a factor with levels \code{Catholic}
#' \code{Jewish} \code{Other} \code{Protestant}}
#' \item{polparty}{a
#' factor with levels \code{Democrat} \code{Independent} \code{Other}
#' \code{Republican}}
#' \item{cappun}{a factor with levels \code{Favor}
#' \code{Oppose} Whether or not the subject favors capital punishment.}
#' \item{tvhours}{the subject estimates number of
#' hours per day he or she watches TV.}
#' \item{marijuan}{a factor with
#' levels \code{Legal} \code{NotLegal} Whether or not subject believes that
#' marijuana should be legalized.}
#' \item{owngun}{a factor with levels
#' \code{No} \code{Yes}.  Does the subject own a gun?}
#' \item{gunlaw}{a factor with levels
#' \code{Favor} \code{Oppose} Whether or not the subject favors stricter
#' gunlaws.}
#' \item{age}{age of the subject}
#' \item{chldidel}{the ideal number of children the subject would like to have.}
#' \item{emailtime}{estimated number of hours per week subject spends using email.}
#' }
#' @source National Opinion Research Center:  \url{http://www3.norc.org/gss+website/}.  
#' Found in Uts and Heckard: Mind on Statistics, 4th Edition.
#' @keywords datasets


NA

#' General Social Survey, 2008
#'
#' 
#' 
#' @name gss08
#' @docType data
#' @format A data frame with 2023 observations on the following 12 variables.
#' \describe{ \item{sex}{a factor with levels \code{Female}
#' \code{Male}} \item{race}{a factor with levels \code{AfrAm}
#' \code{Other} \code{White}} \item{degree}{a factor with levels
#' \code{Bachelor} \code{Graduate} \code{HighSchool} \code{JunColl}
#' \code{NotHs}} \item{relig}{a factor with levels \code{Catholic}
#' \code{Jewish} \code{None} \code{Other} \code{Protestant}}
#' \item{polparty}{a factor with levels \code{Democrat}
#' \code{Independent} \code{Other} \code{Republican}} \item{cappun}{a
#' factor with levels \code{Favor} \code{Oppose}} \item{tvhours}{a
#' numeric vector} \item{marijuan}{a factor with levels \code{Legal}
#' \code{NotLegal}} \item{owngun}{a factor with levels \code{No}
#' \code{Yes}} \item{gunlaw}{a factor with levels \code{Favor}
#' \code{Oppose}} \item{age}{a numeric vector}
#' \item{chldidel}{a numeric vector} }
#' @references For more information see \code{gss02}
#' @source National Opinion Research Center:  \url{http://www3.norc.org/gss+website/}.  
#' Found in Uts and Heckard: Mind on Statistics, 4th Edition.
#' @keywords datasets

NA

#' @title General Social Survey, 2012
#' 
#' @description A selection of variables from the 2012 General Social Survey.  The variables are as follows:
#' 
#' \itemize{
#'   \item age. Age of the subject.
#'   \item sex. Sex of the subject.
#'   \item race. Race of the subject.
#'   \item polviews. Subject's political views.
#'   \item relig. Religion of the subject.
#'   \item cappun.  Opinion on capital punishment.
#'   \item owngun.  Whether or not one owns a gun.
#'   \item emailhr.  Number of hours per week spent on email.
#'   \item bigbang. Whether or not subject believes the Big Bang theory is true.
#'   \item premarsx.  Opinion on premarital sex.
#'   \item pornlaw.  Should pornography be legal?
#'   \item zodiac.  Sign of the Zodiac under which the subject was born.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source \url{http://www3.norc.org/gss+website/}.
#' @format A data frame with 1976 rows and 12 variables
#' @name gss2012

NA

#' Hair Color and ACT Score
#' 
#' A study performed by MAT 111 students at Georgetown College.
#'
#' 
#' @name hair_and_act
#' @docType data
#' @format A data frame with 100 observations on the following 3 variables.
#' \describe{
#' \item{sex}{a factor with levels \code{female} \code{male}}
#' \item{hair.color}{a factor with levels \code{dark} \code{light}}
#' \item{act}{composite ACT score of
#' subject.}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets


NA


#' Height and Handspan
#' 
#' Height and handspan of a few subjects.
#' 
#' 
#' @name handheight
#' @docType data
#' @format A data frame with 167 observations on the following 3 variables.
#' \describe{
#' \item{Sex}{a factor with levels \code{Female} \code{Male}}
#' \item{Height}{height of subject, in inches.}
#' \item{HandSpan}{handspan of subject, in centimeters.}
#' }
#' @source Uts and Heckard, Mind on Statistics, 4th Edition.
#' @keywords datasets


NA

#' Handford Weather Station, 1984-2010
#' 
#' The station is located in Hanford, WA.
#'
#' 
#' @name hanford1
#' @docType data
#' @format A data frame with 27 observations on the following 2 variables.
#' \describe{
#' \item{year}{calendar year}
#' \item{temp}{average high temperature for that year.}
#' }
#' @source See \url{http://www.hanford.gov/hms/}
#' @keywords datasets


NA


#' Hanford Weather Station, 1945-2010
#' 
#' The weather station is located in Hanford, WA.  Note that this dataset is more complete than \code{\link{hanford1}}
#'
#' 
#' @name hanford2
#' @docType data
#' @format A data frame with 66 observations on the following 2 variables.
#' \describe{
#' \item{year}{calendar year}
#' \item{temp}{average high temperature for that year.}
#' }
#' @source See \url{www.hanford.gov/hms/}
#' @keywords datasets


NA


#' Ricky Henderson's Career Stats
#' 
#' Various year-by-year statistics for the Major League player Ricky Henderson
#' 
#' 
#' @name henderson
#' @docType data
#' @format A data frame with 23 observations on the following 18 variables.
#' \describe{ \item{Season}{season}
#' \item{TM}{team played for}
#' \item{G}{number of games played in}
#' \item{AB}{number of at-bats}
#' \item{R}{runs scored}
#' \item{H}{number of base hits}
#' \item{X2B}{doubles}
#' \item{X3B}{triples}
#' \item{HR}{home runs}
#' \item{RBI}{runs batted in}
#' \item{BB}{bases on balls}
#' \item{SO}{number of times struck out}
#' \item{SB}{Number of stolen bases}
#' \item{CS}{number of times caught stealing}
#' \item{AVG}{batting average}
#' \item{OBP}{on-base percentage}
#' \item{SLG}{slugging rate}
#' \item{OPS}{on-base plus slugging}
#'  }
#' @source unknown (possibly Albert: Teaching Statistics With Baseball)
#' @keywords datasets


NA


#' Hall of Fame Batters
#' 
#' Data on the 147 batters inducted into the Major LeagueBaseball Hall of Fame
#' as of the year 2013.
#' 
#' 
#' @name hofbatting
#' @docType data
#' @format A data frame with 147 observations on the following 29 variables.
#' \describe{ \item{Name}{Player name}
#' \item{Rk}{Unknown}
#' \item{Inducted}{Year inducted into Hall of Fame}
#' \item{Yrs}{Number of years planed in the Majors}
#' \item{From}{First year in the Majors}
#' \item{To}{Last year in the Majors}
#' \item{MidCareer}{Middle year of player's career}
#' \item{Era}{Era of Baseball history in which the player was (for the most part) active.
#' Values are:  19th century  (before 1900), Dead Ball (1900-1919), Lively Ball (1920-1940),
#' Integration, (1941-1959), Expansion (1960-1975), Free Agency (1976-1992), Long Ball (1993 +).}
#' \item{ASG}{Number of All-Star games played}
#' \item{WAR.pos}{Wins above replacement (WAR), as defined for a position player}
#' \item{G}{Games played}
#' \item{PA}{Number of plate appearances}
#' \item{AB}{Number of times at bat}
#' \item{R}{Runs scored}
#' \item{H}{Base hits}
#' \item{X2B}{Doubles}
#' \item{X3B}{Triples}
#' \item{Triple.Rate}{Number of triples divided by number of times at bat}
#' \item{HR}{Home runs}
#' \item{HR.Rate}{Number of home runs divided by number of times at bat}
#' \item{RBI}{Runs batted in}
#' \item{SB}{Number of succesful stolen base attempts}
#' \item{CS}{Number of times thrown out while attempting to steal a base}
#' \item{BB}{Base on Balls (number of times "walked")}
#' \item{SO}{Number of times struck out}
#' \item{BA}{Batting average}
#' \item{OBP}{On base percentage}
#' \item{SLG}{Slugging average}
#' \item{OPS}{OBP plus SLG}
#'  }
#' @source Modified from Marchi and Albert:  Analyziing Baseball Data with R, CRC Press 2014.
#' @keywords datasets


NA


#' Hall of Fame Pitchers
#' 
#' Data on the 70 pitchers inducted into the Major League Baseball Hall of Fame
#' as of the year 2013.
#' 
#' 
#' @name hofpitching
#' @docType data
#' @format A data frame with 147 observations on the following 32 variables.
#' \describe{ \item{Name}{Player name}
#' \item{Rk}{Unknown}
#' \item{Inducted}{Year inducted into Hall of Fame}
#' \item{Yrs}{Number of years planed in the Majors}
#' \item{From}{First year in the Majors}
#' \item{To}{Last year in the Majors}
#' \item{MidCareer}{Middle year of player's career}
#' \item{Era}{Era of Baseball history in which the player was (for the most part) active.
#' Values are:  19th century  (before 1900), Dead Ball (1900-1919), Lively Ball (1920-1940),
#' Integration, (1941-1959), Expansion (1960-1975), Free Agency (1976-1992), Long Ball (1993 +).}
#' \item{ASG}{Number of All-Star games played}
#' \item{WAR}{Wins above replacement}
#' \item{W}{Games won}
#' \item{L}{Games lost}
#' \item{W.L.}{proportion of games won}
#' \item{ERA}{Earned run average}
#' \item{G}{Games played}
#' \item{GS}{Games started}
#' \item{GF}{Games finished}
#' \item{CG}{Complete games}
#' \item{SHO}{Shut-outs}
#' \item{SV}{Saves}
#' \item{IP}{Innings pitched}
#' \item{H}{Hits allowed}
#' \item{R}{Runs allowed}
#' \item{ER}{Earned Runs allowed}
#' \item{HR}{Home Runs allowed}
#' \item{BB}{Bases on Balls (number of "walks")}
#' \item{IBB}{Intentional walks}
#' \item{SO}{Strikeouts}
#' \item{HBP}{Hit batter with pitch (?)}
#' \item{BK}{Balks}
#' \item{WP}{Wild Pitches}
#' \item{BF}{Total batters faced}
#'  }
#' @source Modified from Marchi and Albert:  Analyziing Baseball Data with R, CRC Press 2014.
#' @keywords datasets


NA


#' @title An Imaginary Population
#' 
#' @description An imaginary population, used for instructional purposes.  The variables are as follows:
#' 
#' \itemize{
#'   \item sex. (male, female).  
#'   \item math. Whether or not you were a mathematics major.
#'   \item income.  Your annual income, rounded to the nearest $100.
#'   \item cappun.  Opinion about the death penalty (favor, oppose).
#'   \item height.  Height in inches.
#'   \item idealheight.  The height you would like to be, in inches.
#'   \item diff.  Idealheight - actual height.
#'   \item kkardashtemp.  Your feelings about Kim Kardashian on a 0-100 scale  (0=very cold, 100=very warm).
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 10000 rows and 8 variables
#' @name imagpop


NA


#' IQs of Siblings
#' 
#' IQs of pairs of siblings.
#'
#' 
#' @name iqsiblings
#' @docType data
#' @format A data frame with 80 observations on the following 2 variables.
#' \describe{
#' \item{First}{IQ of the older sibling.}
#' \item{Second}{IQ of the younger sibling.}
#' }
#' @source William Harris, Georgetown College
#' @keywords datasets


NA


#' Knife or Gun?
#' 
#' What will make you yell louder: being killed with a knife or being killed
#' with a gun?  Results of an entirely imaginary experiment performed on very
#' strange volunteers.  Members of the Knife group are killed by a knife, and
#' members of the Gun group are killed by a gun. The volume of the screams of each
#' subject during slaying is recorded.  In order to ensure that the two groups
#' are similar with respect to how loud they can yell in the first place,
#' subjects are blocked by whether or not they have participated in
#' hog-hollering contests.  After blocking, subjects are randomly assigned to
#' groups.
#' 
#' 
#' @name knifeorgunblock
#' @docType data
#' @format A data frame with 20 observations on the following 3 variables.
#' \describe{
#' \item{hoghollerer}{a factor with levels \code{no}
#' \code{yes} whether or not the subject competes in hog-hollerin' contests}
#' \item{means}{a factor with levels \code{gun}
#' \code{knife} means by which subject is slain}
#' \item{volume}{volume of expiring subject's cries.}
#' }
#' @source A morbid imagination.
#' @keywords datasets


NA


#' Labels and Perception of Quality
#' 
#' Students in MAT 111 performed an experiment to see whether the perception of
#' the quality of peanut butter was affected by the labeling on the peanut
#' butter jar.  Each subject tasted from two jars, one of which was labeled
#' Jiff, and the other of which was labeled Great Value (a cheaper brand).
#' Unknown to the subjects, both jars contained Great Value peanut butter.
#' Each subject rated the quality of the peanut butter on a scale of 1 to 10.
#' 
#' 
#' @name labels
#' @docType data
#' @format A data frame with 30 observations on the following 3 variables.
#' \describe{
#' \item{jiffrating}{rating subject gave to the PB in the jar with the Jiff label}
#' \item{greatvaluerating}{rating subject gave to the PB in the jar with the Great Value label}
#' \item{sex}{a factor with levels \code{female} \code{male}}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets


NA


#' @title Crowd Behavior at Ledge-Jumping Incidents
#' 
#' @description A dataset recreated from summary data describing the relationship between weather and crowd behavior during 21 recorded  
#' incidents in England, in which a (suicidal) person was contemplating jumping from a ledge or other high structure and a crowd 
#' gathered to watch.  The variables are as follows:
#' 
#' \itemize{
#'   \item weather. Warm or cool, based on the time of year when the incident occurred.  
#'   \item crowd.behavior. The crowd either baited the would-be jumper, or was polite.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source "The baiting crowd in episodes of threatened suicide",
#' \emph{Journal of Personality and Social Psychology}, 41, 703-709.  See also dataset 59 in
#' \emph{A Handbook of Small Datasets} by Hand et. al.  See also \url{http://www.ncbi.nlm.nih.gov/pubmed/7288565}.
#' @format A data frame with 21 rows and 2 variables
#' @name ledgejump


NA


#' @title MAT 111 Survey
#' 
#' @description Results of a survey of MAT 111 students at Georgetown College.
#' 
#' \itemize{
#'       \item height.  How tall are you, in inches?
#'       \item ideal_ht.  A numeric vector  How tall would you LIKE to be, in inches?
#'       \item sleep.  How much sleep did you get last night?
#'       \item fastest.  What is the highest speed at which you have ever driven a car?
#'       \item weight_feel.  How do you feel about your weight?
#'       \item love_first.  Do you believe in love at first sight?
#'       \item extra_life.  Do you believe in extraterrestrial life?
#'       \item seat.  When you have a choice, where do you prefer to sit in a classroom?
#'       \item GPA.  What is your college GPA?
#'       \item enough_Sleep.  Do you think you get enough sleep?
#'       \item sex.  What sex are you?
#'       \item diff.  Your ideal height minus your actual height.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Georgetown College, MAT 111.
#' @format A data frame with 71 rows and 12 variables
#' @name m111survey


NA


#' MAT 111 Survey, Fall 2012
#' 
#' Results of a survey given at beginning of semester, to all students in MAT 111.
#'
#' 
#' @name m111surveyfa12
#' @docType data
#' @format A data frame with 89 observations on the following 14 variables.
#' \describe{
#' \item{height}{Your height in inches.}
#' \item{ideal_ht}{How tall you would LIKE to be, in
#' inches.}
#' \item{sleep}{How much sleep you got last
#' night, in hours.}
#' \item{fastest}{What is the
#' highest speed at which you have ever driven a car (in mph)?}
#' \item{wt_feel}{a factor with levels \code{1_underweight}
#' \code{2_about_right} \code{3_overweight} How do you feel about your weight?}
#' \item{love_first}{a factor with levels \code{no} \code{yes} Do you
#' believe in love at first sight?}
#' \item{et_life}{a factor with levels
#' \code{no} \code{yes} Do you believe in life on other planets?}
#' \item{seat}{a factor with levels \code{1_front}
#' \code{2_middle} \code{3_back} When you have a choice, where do you prefer to
#' sit in a classroom?}
#' \item{GPA}{What is your
#' current GPA?}
#' \item{engh_slp}{a factor with levels \code{no}
#' \code{yes} Do you think you get enough sleep?}
#' \item{sex}{a factor
#' with levels \code{female} \code{male} What sex are you?}
#' \item{anchor}{a factor with levels \code{australia}
#' \code{united_states} (Anchor for the next question.)  For the next question,
#' either Australia or the US, along with its population, was given in the
#' leadup information to the question.  The "anchor" variable records which version of the question you were
#' given.}
#' \item{canada}{"The population of country
#' XXX is YYY million.  About what is the population of Canada, in millions?"  XXX was either the U.S. or Australia.}
#' \item{diff.ih.ah.}{Your ideal height minus your actual height.}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets


NA


#' MAT 111 Survey, Fall 2013
#' 
#' Results of a survey given at beginning of semester, to all students in MAT 111.
#'
#' 
#' @name m111surveyfa13
#' @docType data
#' @format A data frame with 85 observations on the following 14 variables.
#' \describe{
#' \item{height}{Your height in inches.}
#' \item{ideal_ht}{How tall you would LIKE to be, in
#' inches.}
#' \item{sleep}{How much sleep you got last
#' night, in hours.}
#' \item{fastest}{What is the
#' highest speed at which you have ever driven a car (in mph)?}
#' \item{weight_feel}{a factor with levels \code{1_underweight}
#' \code{2_about_right} \code{3_overweight} How do you feel about your weight?}
#' \item{love_first}{a factor with levels \code{no} \code{yes} Do you
#' believe in love at first sight?}
#' \item{extra_life}{a factor with levels
#' \code{no} \code{yes} Do you believe in life on other planets?}
#' \item{seat}{a factor with levels \code{1_front}
#' \code{2_middle} \code{3_back} When you have a choice, where do you prefer to
#' sit in a classroom?}
#' \item{GPA}{What is your
#' current GPA?}
#' \item{enough_Sleep}{a factor with levels \code{no}
#' \code{yes} Do you think you get enough sleep?}
#' \item{sex}{a factor
#' with levels \code{female} \code{male} What sex are you?}
#' \item{diff}{ideal height minus actual height}
#' \item{symbol}{a factor with levels \code{a}
#' \code{b} (Anchor for the next question.)  For the next question,
#' either Australia or the US, along with its population, was given in the
#' leadup information to the question.  The "anchor" variable records which version of the question you were
#' given.  If "a", the population of Australia was given.  If "b", the U.S. population was given.}
#' \item{pop_Canada}{"The population of country
#' XXX is YYY million.  About what is the population of Canada, in millions?"  XXX was either the U.S. or Australia.}
#' }
#' @source MAT 111 at Georgetown College
#' @keywords datasets


NA


#' @title Music and Reading Comprehension
#' 
#' @description An experiment performed by a student at Georgetown College.  Forty-four subjects were randomized into four
#' groups.  All subjects read an article; one group read in a silent environment, while the other three groups
#' heard each three different genres of music.  Each subject took a reading comprehension test afterward.
#' 
#' \itemize{
#' \item{\code{sex}} {a factor with levels \code{Female} \code{Male}}
#' \item{\code{year}} {class rank of subject}
#' \item{\code{type}} {type of music subject listened to while reading}
#' \item{\code{score}} {number of questions correct on reading comprehension test}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Matt Doolin, MAT 111 at Georgetown College
#' @format A data frame with 44 observations on 4 variables.
#' @name music


NA


#' Napkin Use
#' 
#' Students at GC observed their fellow students in the Cafe at lunch.
#'
#' 
#' @name napkins
#' @docType data
#' @format A data frame with 86 observations on the following 2 variables.
#' \describe{ \item{napkins}{number of napkins used
#' by the subject during the meal.} \item{sex}{a factor with levels
#' \code{female} \code{male} Sex of the person being observed}}
#' @source MAT 111 at Georgetown College
#' @keywords datasets


NA


#' Non-Response to Surveys
#' 
#' Results of a study on non-response to a mail survey.  Subjects were residents of Denmark.
#' 
#' 
#' @name nonresponse
#' @docType data
#' @format A data frame with 4229 observations on the following 3 variables.
#' \describe{ \item{residence}{where the subject resides:  either in Copenhagen, a city outside of Copenhagen, or 
#' in the countryside} 
#' \item{gender}{sex of the subject}
#' \item{response}{Whether or not eh subject responded to the mail survey}}
#' @source Rebuilt from a contingency table in E. B. Andersen (1991), 
#' The Statistical Analysis of Categorical Data. 2nd edition. Springer-Verlag, Berlin.  Table found in
#' package \code{vcd}.
#' @keywords datasets


NA


#' Nicotine Withdrawal and Accidents
#' 
#' Results of study conducted in Great Britain to see if nicotine withdrawal increases the risk of an accident.
#' 
#' 
#' @name nosmokeday
#' @docType data
#' @format A data frame with 10 observations on the following 3 variables.
#' \describe{ \item{year}{calendar year}
#' \item{Injuries.before.NSD}{number of injury accidents on the day one week prior to 
#' National No Smoke Day in the United Kingdom} 
#' \item{Injuries.on.NSD}{number of injury accidents on  
#' National No Smoke Day in the United Kingdom}}
#' @source J. Knwles, "Nicotine withdrawal and road accidents", Science, 400, 128, (8 July 1999).
#' Found in Whitlock and Schluter, The Analysis of Biological Data.
#' @keywords datasets


NA


#' Old Faithful
#' 
#' Old faithful geyser at Yellowstone Park.
#'
#' 
#' @name oldfaithful
#' @docType data
#' @format A data frame with 299 observations on the following 2 variables.
#' \describe{
#' \item{Duration}{duration of eruption, in minutes}
#' \item{TimeNext}{time until the next eruption, in minutes}
#' }
#' @source Unknown
#' @keywords datasets


NA


#' Ovarian Cancer Study
#' 
#' Results of a retrospective study, conducted in 1973,
#' on 299 women who been surgically treated for ovarian cancer 10 years before.
#' 
#' 
#' @name Ovarian
#' @docType data
#' @format A data frame with 299 observations on the following 4 variables.
#' \describe{ \item{stage}{factor indicating the stage of the cancer at the time of operation (early, advanced)} 
#' \item{operation}{factor indicating the amount of tissue removed during surgery (radical,limited)}
#' \item{survival}{whether or not the subject was still alive after ten years (yes,no)}
#' \item{xray}{factor indicating whether or not the subject also received x-ray treatments (yes,no)}}
#' @source Rebuilt from a contingency table in E. B. Andersen (1991), 
#' The Statistical Analysis of Categorical Data. 2nd edition. Springer-Verlag, Berlin.  Table found in
#' package \code{vcd}.
#' @keywords datasets


NA


#' Penn State #1
#' 
#' A study of students at Penn State University.
#'
#' 
#' @name pennstate1
#' @docType data
#' @format A data frame with 190 observations on the following 9 variables.
#' \describe{
#' \item{Sex}{a factor with levels \code{F} \code{M}}
#' \item{HrsSleep}{how many hours of sleep the subject gets per night}
#' \item{SQpick}{a factor
#' with levels \code{Q} \code{S}.  Each subject was presented with two letters
#' (S and Q), and asked to pick one.  This variable indicates which letter the
#' subject picked.}
#' \item{Height}{height in inches}
#' \item{RandNumb}{a numeric vector: Each subject was asked to choose
#' randomly an integer from 1 to 10.}
#' \item{Fastest}{highest speed, in mph, at which subject has ever driven a car}
#' \item{RtSpan}{span of the right hand, in
#' centimeters.}
#' \item{LftSpan}{span of the left
#' hand, in centimeters.}
#' \item{Form}{a factor with levels \code{QorS}
#' \code{SorQ}.  The order of presentation of the S and Q options to the
#' subject varied from one survey form to another.  This variable indicates
#' which letter was presented first on the form.}
#' }
#' @source Uts and Heckard, Mind on Statistics, 4th Edition.
#' @keywords datasets


NA


#' Pushups by Football Players at Georgetown College
#' 
#' Two football players at GC asked their team-mates to do as many pushups as
#' they could in two minutes.
#' 
#' 
#' @name pushups
#' @docType data
#' @format A data frame with 30 observations on the following 3 variables.
#' \describe{
#' \item{weight}{weight of subject in
#' pounds.}
#' \item{pushups}{number of pushups
#' completed.}
#' \item{position}{a factor with levels \code{LINE}
#' \code{SKILL}: type of position played by the subject.  Line positions
#' require high body mass, skill positions require a lot of running around.}
#' }
#' @source MAT 111, Georgetown College
#' @keywords datasets


NA


#'Effect of Soil Salinity on Plant Growth
#' 
#' @description Result of an experiment conducted to investigate the effect of salinity level in soil 
#' on the growth of plants.
#' 
#' @details  From the source (see below):  "Experimental
#' fields of land were located at an agricultural field station, and each field was divided into six smaller
#' plots. Each of the smaller plots was treated with a different amount of salt (measured in ppm) and the
#' biomass at the end of the experiment was recorded."
#' 
#' 
#' @name saltmarsh
#' @docType data
#' @format A data frame with 24 observations on the following 3 variables.
#' \describe{ \item{salt}{amount of salt applied to the plot (in parts per million)}
#' \item{biomass}{total biomass of plot at the end of the study period (units unknown)}
#' \item{block}{field in which the plot was located} }
#' @source The Course Notes of Carl Schwarz, Simon Fraser University:  
#' \url{http://people.stat.sfu.ca/~cschwarz/CourseNotes/}
#' @keywords datasets


NA


#' @title SAT Scores
#' 
#' @description SAT scores by state.  The variables are as follows:
#' 
#' \itemize{
#'   \item state.  A state in the U.S.
#'   \item expend.  Mean annual expenditure per student (in 1000$).
#'   \item ratio.  Mean student-teacher ratio.
#'   \item salary.  Mean annual teacher salary.
#'   \item frac.  Percentage of students in the state who take the SAT.
#'   \item verbal.  Mean SAT Verbal score for the state.
#'   \item math.  Mean SAT Math score for the state.
#'   \item sat. Sum of mean Verbal and mean Math.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Deborah Lynn Guber, "Getting what you pay for: the debate over equity 
#' in public school expenditures" (1999), \emph{Journal of Statistics Education} 7(2).
#' \url{http://www.amstat.org/publications/jse/secure/v7n2/datasets.guber.cfm}.
#' @format A data frame with 50 rows and 8 variables
#' @name sat


NA


#' Weddell Seal Oxygen Consumptions
#' 
#' Results of an experiment conducted on ten Weddell seals.
#' 
#' 
#' @name sealsO2
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{o2.nonfeeding}{Oxygen consumption during recovery time after a dive during which
#'  no plankton was consumed by the seal, in ml of O2 per kilogram of weight} 
#' \item{o2.feeding}{Oxygen consumption during recovery time after a dive during which
#'  plankton was consumed, in ml of O2 per kilogram of weight}}
#' @source Williams, T. M., L. A. Fuiman, M. Horning, and R. W. Davis. 2004. 
#' The Journal of Experimental Biology 207: 973 to 982.  \url{http://jeb.biologists.org/content/207/6/973.full}
#' @keywords datasets


NA

#' A Small Experiment
#' 
#' Subjects in a hypothetical experiment
#' 
#' 
#' @name SmallExp
#' @docType data
#' @format A data frame with 16 observations on the following 3 variables.
#' \describe{ \item{name}{name of the subject}
#' \item{sex}{sex of the subject}
#' \item{athlete}{whether or not the subject is an athlete} }
#' @keywords datasets


NA


#' Larvae on Stumps
#' 
#' Biologists were interested in whether beetles prefer areas where beavers have cut down cottonwood trees.
#' (The tree-stumps produce tender green shoots that beetles are thought to like.)  23 circular plots, all of equal area,
#' were studied.  For each plot the researchers counted the number of cottonwood stumps, and also the number of
#' clusters of beetle larvae found in the plot.
#' 
#' 
#' @name stumps
#' @docType data
#' @format A data frame with 23 observations on the following 2 variables.
#' \describe{
#' \item{stumps}{number of stumps in the plot}
#' \item{larvae}{number of larvae clusters in the plot}
#' }
#' @source Basic Practice of Statistics, by Moore and McCabe.
#' @keywords datasets


NA


#'@title Temperature in U.S. Cities
#' 
#' @description Average temperatures for cities in the United States.
#' 
#' \itemize{
#' \item{\code{city}} {Name of the city}
#' \item{\code{latitude}} {latitude of the city, in degrees north of the Equator}
#' \item{\code{JanTemp}} {mean temperature of the city in January.}
#' \item{\code{AprTemp}} {mean temperature of the city in April.}
#' \item{\code{AugTemp}} {mean temperature of the city in August.}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Mind on Statistics, 4th edition, Uts and Heckard.
#' @format A data frame with 20 observations on 5 variables.
#' @name temperature


NA

#' @title Restaurant Tips
#' 
#' @description A waiter recorded all tips earned during a 2.5 month period in early 1990, along with other
#' information about the customers who gave the tips.  The variables are as follows:
#' 
#' \itemize{
#'   \item {obs} {Observation number of event}
#'   \item {totbill} {The total bill}
#'   \item {tip} {Amount of the tip}
#'   \item {sex} {Sex of tipper (F or M)}
#'   \item {smoker} {Whether the tipper smokes: Yes or No}
#'   \item {day} {Day of the week}
#'   \item {time} {Whether the meal was during the Day or the Night}
#'   \item {size} {Number of people in the dining party}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Bryant, P. G. and Smith, M. A.(1995), Practical Data Analysis: Case
#` Studies in Business Statistics, Richard D. Irwin Publishing, Homewood, IL.
#' Found in Interactive and Dynamic Graphics for Data Analysis: 
#' With Examples Using R and GGobi, Dianne Cook and Deborah F. Swayne.
#' See also \url{http://www.ggobi.org/book/}.
#' @format A data frame with 244 rows and 8 variables
#' @name tips


NA

#'@title Tornado Damage
#' 
#' @description Tornado damage in the U.S., by state.  Also includes Puerto Rico.
#'
#' \itemize{ 
#' \item{\code{state}} {the state}
#' \item{\code{damage}} {mean annual damage from tornados, over a five-year period, in millions of dollars}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Moore and McCabe, The Basic Practice of Statistics.
#' @format A data frame with 51 observations on 2 variables.
#' @name tornado


NA


#' @title UC Davis #1
#' 
#' @description Results of a survey of students at UC-Davis.
#' 
#' \itemize{
#' \item{\code{Sex}}{ a factor with levels \code{Female} \code{Male}}
#' \item{\code{TV}}{ Number of hours spent watching TV per week}
#' \item{\code{computer}}{ number of hours spent on computer per week}
#' \item{\code{Sleep}}{ hours of sleep per night}
#' \item{\code{Seat}}{ a factor with levels \code{Back} \code{Front} \code{Middle} Where do you prefer to sit in class,
#' when you have a choice?}
#' \item{\code{alcohol}}{ number of alcoholic drinks consumed per week}
#' \item{\code{Height}}{ height in inches}
#' \item{\code{momheight}}{ height of mother, in inches}
#' \item{\code{dadheight}}{ height of father, in inches}
#' \item{\code{exercise}}{ number of hours of exercise per week}
#' \item{\code{GPA}}{ grade point average}
#' \item{\code{class}}{ a factor with levels \code{LibArts} \code{NonLib} Student Category:  liberal arts or
#' not}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Mind on Statistics, 4th edition, Uts and Heckard.
#' @format A data frame with 173 observations on 12 variables..
#' @name ucdavis1


NA


#' Justin Verlander
#' 
#' PITCHf/x data on Justin Verlander, winner of the 2011 Cy Young Award.  Covers his 
#' 2009 through 2012 seasons.
#' 
#' 
#' @name verlander
#' @docType data
#' @format A data frame with 15307 observations on the following 12 variables.  Each observation
#' is a single pitch.
#' \describe{ \item{season}{The year of play}
#' \item{gamedate}{Date of the game in which the pitch was thrown}
#' \item{pitch_type}{Type of pitch thrown:  CH (Change-up), CU (Curveball), FF (Four-Seam Fastball),
#' FT (Two-Seam Fastball), SL (Slider).  Pitch type is determined by computer algorithm.}
#' \item{balls}{Current ball count}
#' \item{strikes}{Current strike count}
#' \item{pitches}{number of pitches previously thrown in the game}
#' \item{speed}{speed of pitch (in mph).  (When crossing plate?)}
#' \item{px}{x-coordinate of pitch (in feet, measured from center of plate)}
#' \item{pz}{vertical coordinate of pitch (in feet above plate)}
#' \item{pfx_x}{the horizontal movement, in inches, of the pitch between the release point and 
#' home plate, as compared to a theoretical pitch thrown at the same speed with no spin-induced 
#' movement. Measured at 40 feet from home plate.}
#' \item{pfx_z}{the vertical movement, in inches, of the pitch between the release point and 
#' home plate, as compared to a theoretical pitch thrown at the same speed with no spin-induced 
#' movement. Measured at 40 feet from home plate.}
#' \item{batter_hand}{A factor with two values: L (left) and R (right).}
#'  }
#' @source Marchi and Albert:  Analyziing Baseball Data with R, CRC Press 2014.  For
#' more on the PITCHf/x system, see \url{http://en.wikipedia.org/wiki/PITCHf/x}.
#' 
#' @keywords datasets


NA


#' @title Youth Risk 2003
#' 
#' @description A Study of Risky Behaviors in High School Seniors, from year 2003.
#' 
#' \itemize{
#' \item \code{Sex}
#' \item \code{Grades} Typical grades you earn in school
#' \item \code{WtAction} What do you plan to do about your weight?
#' \item \code{Seatbelt} How often do you wear a seat-belt?
#' \item \code{Sunscreen} How much do you wear sunscreen?
#' \item \code{Grades_1} Same as grades, but with some groups combined
#' \item \code{Sun_1} Same as sunscreen, but with some groups combined
#' }
#' 
#' @docType data
#' @keywords datasets
#' @source Mind on Statistics, 4th Edition, by Uts and Heckard.
#' @format A data frame with 3042 observations on 7 variables.
#' @name youthrisk03

NA
