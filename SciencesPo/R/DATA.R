#' @encoding UTF-8
#' @title Galton's Family Data on Human Stature.
#'
#' @description It is a reproduction of the data set used by Galton in his 1885's paper on correlation between parent's height and their children. However, Galton would only introduce the concept of correlation few years later, in 1888. Galton suggested the use of the regression line and was the first to describe the so-called common phenomenon of regression toward the mean by comparing his experiments on the size of the seeds of successive generations of peas. This dataset contains the following columns:
#'
#' \itemize{
#'   \item parent the parents' average height
#'   \item child the child's height
#' }
#'
#' @details Regression analysis is the statistical method most often used in political science research. The reason is that most scholars are interested in identifying \dQuote{causal} effects from non-experimental data and that regression is the method for doing this. The term \dQuote{regresssion} (1889) was first crafted by Sir Francis Galton upon investigating the relationship between body size of fathers and sons. Thereby he \dQuote{invented} regression analysis by estimating: \eqn{S_s = 85.7 + 0.56S_F} meaning that the size of the son regresses towards the mean.
#'
#' @references  Francis Galton (1886) Regression Towards Mediocrity in Hereditary Stature. \emph{The Journal of the Anthropological Institute of Great Britain and Ireland,} Vol. \bold{15}, pp. 246--263.
#'
#' @docType data
#' @keywords datasets
#' @name galton
#' @usage data(galton)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::galton)} variables and \Sexpr{nrow(SciencesPo::galton)} observations.
NULL




#' @encoding UTF-8
#' @title Same Sex Marriage Public Opinion Data
#'
#' @description Data set fielded by the PEW Research Center on same sex marriage support in US. It covers public opinion on the issue starting from 1996 up to date. This dataset contains the following columns:
#'
#' \itemize{
#'   \item Date The year of the measurement
#'   \item Oppose Percent opposing same-sex marriage
#'   \item Favor Percent favoring same-sex marriage
#'   \item DK Percent of Don't Know
#' }
#'
#' @references PEW Research Center. \emph{Support for same-sex marriage}. \url{http://www.pewresearch.org}
#'
#' @docType data
#' @keywords datasets
#' @name marriage
#' @usage data(marriage)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::marriage)} variables and \Sexpr{nrow(SciencesPo::marriage)} observations.
NULL



#' @encoding UTF-8
#' @title The Penn World Table
#'
#' @description The Penn World Table used in Summers and Heston (1991). This dataset contains the following columns:
#'
#' \itemize{
#' \item year Year
#' \item pop Population (thousands)
#' \item rgdppc Real per capita GDP
#' \item savrat a numeric vector
#' \item country Country
#' \item com Communist regime
#' \item opec OPEC country
#' \item name Country name
#' }
#'
#' @references Summers, R. and Heston, A. (1991) The Penn World Table (Mark 5): an expanded set of international comparisons, 1950--1988. \emph{The Quarterly Journal of Economics,} \bold{106(2),} 327--368.
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#'
#' @docType data
#' @keywords datasets
#' @name sheston91
#' @usage data(sheston91)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::sheston91)} variables and \Sexpr{nrow(SciencesPo::sheston91)} observations.
NULL




#' @encoding UTF-8
#' @title Griliches's (1976) Data
#'
#' @description Data set used by Griliches (1976) on wages of very young men. This dataset contains the following columns:
#'
#' \itemize{
#' \item rns residency in South.
#' \item rns80 a numeric vector.
#' \item mrt marital status = 1 if married.
#' \item mrt80 a numeric vector.
#' \item smsa reside metro area = 1 if urban.
#' \item smsa80 a numeric vector.
#' \item med mother's education, years.
#' \item iq iq score.
#' \item kww score on knowledge in world of work test.
#' \item year Year.
#' \item age a numeric vector.
#' \item age80 a numeric vector.
#' \item s completed years of schooling.
#' \item s80 a numeric vector.
#' \item expr experience, years.
#' \item expr80 a numeric vector.
#' \item tenure tenure, years.
#' \item tenure80 a numeric vector.
#' \item lw log wage.
#' \item lw80 a numeric vector.
#' }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Griliches, Z. (1976) Wages of very young men. \emph{The Journal of Political Economy,} \bold{84(4),} 69--85.
#'
#' @docType data
#' @keywords datasets
#' @name griliches76
#' @usage data(griliches76)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::griliches76)} variables and \Sexpr{nrow(SciencesPo::griliches76)} observations.
NULL



#' @encoding UTF-8
#' @title Marc Nerlove's (1963) data
#'
#' @description Data used by Marc Nerlove (1963) on returns of electricity supply. This dataset contains the following columns:
#'
#' \itemize{
#' \item totcost costs in 1970, MM USD.
#'  \item output output, billion KwH.
#' \item plabor price of labor.
#'  \item pfuel price of fuel.
#'  \item pkap price of capital.
#'  }
#'
#'
#'#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Nerlove, M. (1963) Returns to Scale in Electricity Supply. In \emph{Measurement in Economics-Studies in Mathematical Economics and Econometrics in Memory of Yehuda Grunfeld,} edited by Carl F. Christ. Stanford: Stanford.
#'
#' @docType data
#' @keywords datasets
#' @name nerlove63
#' @usage data(nerlove63)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::nerlove63)} variables and \Sexpr{nrow(SciencesPo::nerlove63)} observations.
NULL




#' @encoding UTF-8
#' @title Stock's and Watson's (1993) Data.

#' @description Data set used by Stock and Watson (1993) to estimate co-integration. This dataset contains the following columns:
#' \itemize{
#' \item lnm1 Log M1.
#' \item lnp Log NNP price deflator.
#' \item lnnnp Log NNP.
#' \item cprate A numeric vector.
#' \item year Commercial paper rate.
#'  }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Stock, J. H., and Watson, M. W. (1993) A simple estimator of cointegrating vectors in higher order integrated systems. \emph{Econometrica: Journal of the Econometric Society,} 783--820.
#'
#' @docType data
#' @keywords datasets
#' @name swatson93
#' @usage data(swatson93)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::swatson93)} variables and \Sexpr{nrow(SciencesPo::swatson93)} observations.
NULL




#' @encoding UTF-8
#' @title Lothian's and Taylor's (1996) Data Set

#' @description Data used by Lothian and Taylor (1996) on the real exchange rate behaviour. This dataset contains the following columns:
#'
#' \itemize{
#' \item year Year
#' \item spot dollar/sterling exchange rate.
#' \item USwpi U.S. wholesale price index, 1914==100.
#' \item UKwpi U.K. wholesale price index, 1914==100.
#' }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#'
#' @references Lothian, J. R., and Taylor, M. P. (1996) Real exchange rate behavior: the recent float from the perspective of the past two centuries. \emph{Journal of Political Economy,} 488--509.
#'
#' @docType data
#' @keywords datasets
#' @name ltaylor96
#' @usage data(ltaylor96)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::ltaylor96)} variables and \Sexpr{nrow(SciencesPo::ltaylor96)} observations.
NULL




#' @encoding UTF-8
#' @title Christensen's and Greene's (1976) Data

#' @description Data set used by Christensen and Greene (1976) on economies of scale in US electric power generation. This dataset contains the following columns:
#'
#' \itemize{
#' \item firmid Observation id.
#' \item costs Costs in 1970, MM USD.
#' \item output Output, million KwH.
#' \item plabor Price of labor.
#' \item pkap Price of capital.
#' \item pfuel Price of fuel.
#' \item labshr Labor's cost share.
#' \item kapshr Capital's cost share.
#' }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Christensen, L. R., and Greene, W. H. (1976) Economies of scale in US electric power generation. \emph{The Journal of Political Economy,} 655--676.
#'
#' @docType data
#' @keywords datasets
#' @name cgreene76
#' @usage data(cgreene76)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::cgreene76)} variables and \Sexpr{nrow(SciencesPo::cgreene76)} observations.
NULL




#' @encoding UTF-8
#' @title Mishkin's (1992) Data
#'
#' @description Data from the Frederic S. Mishkin (1992) paper \dQuote{Is the Fisher Effect for real?}. This dataset contains the following columns:
#'
#' \itemize{
#'  \item year Year
#' \item mon a numeric vector
#' \item inf1mo a numeric vector
#' \item inf3mo a numeric vector
#' \item tbill1mo a numeric vector
#' \item tbill3mo a numeric vector
#' \item cpiu a numeric vector
#' \item quote a numeric vector
#' }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Mishkin, F. S. (1992) Is the Fisher effect for real?: A reexamination of the relationship between inflation and interest rates. \emph{Journal of Monetary Economics,} \bold{30(2),} 195--215.
#'
#' @docType data
#' @keywords datasets
#' @name mishkin92
#' @usage data(mishkin92)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::mishkin92)} variables and \Sexpr{nrow(SciencesPo::mishkin92)} observations.
NULL



#' @encoding UTF-8
#' @title Bekaert's and Hodrick's (1993) Data

#' @description Data set used by Bekaert and Hodrick (1993) on biases in the measurement of foreign exchange risk premiums. This dataset contains the following columns:
#'
#' \itemize{
#'  \item date A character vector for date.
#'  \item jyspot Price of USD in JY, spot.
#'  \item jyfwd Price of USD in JY, 30-day forward.
#'  \item jys30 Price of USD in JY, spot market at 30-day forward deliver/date.
#'  \item dmspot Price of USD in DM, spot.
#'  \item dmfwd Price of USD in DM, 30-day forward
#'  \item dms30 Price of USD in DM, spot market at 30-day forward deliver/date.
#'  \item bpspot Price of USD in BP, spot.
#'  \item bpfwd Price of USD in BP, 30-day forward.
#'  \item bps30 Price of USD in BP, spot market at 30-day forward deliver/date.
#' \item quote A numeric vector.
#' }
#'
#' @source Hayashi, F. (2000). \emph{Econometrics.} Princeton. New Jersey, USA: Princeton University. \url{http://fmwww.bc.edu/ec-p/data/Hayashi/}
#'
#' @references Bekaert, G., and Hodrick, R. J. (1993) On biases in the measurement of foreign exchange risk premiums. \emph{Journal of International Money and Finance,} \bold{12(2),} 115-138.
#'
#' @docType data
#' @keywords datasets
#' @name bhodrick93
#' @usage data(bhodrick93)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::bhodrick93)} variables and \Sexpr{nrow(SciencesPo::bhodrick93)} observations.
NULL



#' @encoding UTF-8
#' @title Titanic

#' @description Population at Risk and Death Rates for an Unusual Episode. For each person on board the fatal maiden voyage of the ocean liner Titanic, this dataset records sex, age [adult/child], economic status [first/second/third class, or crew] and whether or not that person survived. This dataset contains the following columns:
#'
#' \itemize{
#' \item CLASS Class (0 = crew, 1 = first, 2 = second, 3 = third)
#' \item AGE Age   (1 = adult, 0 = child)
#' \item SEX Sex   (1 = male, 0 = female)
#' \item SURVIVED Survived (1 = yes, 0 = no)
#' }
#'
#' @note There is not complete agreement among primary sources as to the exact numbers on board, rescued, or lost. \bold{STORY BEHIND THE DATA:} The sinking of the Titanic is a famous event, and new books are still being published about it.  Many well-known facts--from the proportions of first-class passengers to the "women and children first" policy, and the fact that that policy was not entirely successful in saving the women and children in the third class--are reflected in the survival rates for various classes of passenger.  These data were originally collected by the British Board of Trade in their investigation of the sinking.
#'
#' @source British Board of Trade Inquiry Report (1990). \emph{Report on the Loss of the `Titanic' (S.S.)"}, Gloucester, UK: Allan Sutton Publishing.
#'
#' @references Dawson (1995). "The `Unusual Episode' Data Revisited" in the \emph{Journal of Statistics Education}.
#'
#' @docType data
#' @keywords datasets
#' @name titanic
#' @usage data(titanic)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::titanic)} variables and \Sexpr{nrow(SciencesPo::titanic)} observations.
NULL




#' @encoding UTF-8
#' @title Turnout Data
#'
#' @description Data on voter turnout in the 50 states and D.C. for the 1992
#' Presidential election and 1990 Congressional elections. Per capita income,
#' populations in poverty and populations with no high school degree are also given. This dataset contains the following columns:
#'
#' \itemize{
#' \item v1 state name (alphabetic, 20 characters)
#' \item v2 region of country:
#' 1 Northeast
#' 2 Midwest
#' 3 South
#' 4 West
#' \item v3 division within region:
#' 1 Northeast-New England
#' 2 Northeast-Middle Atlantic
#' 3 Midwest-East North Central
#' 4 Midwest-West North Central
#' 5 South-South Atlantic
#' 6 South-East South Central
#' 7 South-West South Central
#' 8 West-Mountain
#' 9 West-Pacific.
#' \item v4 "Elazar's state political culture assignments:
#' 1 = moralistic
#' 2 = individualistic
#' 3 = traditionalistic
#' \item v5 percent of population below poverty level, 1992.
#' \item v6 per capita personal income, 1993.
#' \item v7 percent casting votes for presidential electors, 1992.
#' \item v8 percent casting votes for U.S. Representatives, 1990.
#' \item v9 population without a high school degree, of those 25 years or older, 1990.
#' \item v10 population 25 years or older, 1990.
#' \item v11 South =1; all others =0.
#'  }
#'
#' @source U.S. Bureau of the Census, Statistical Abstract of the United States, 1994.
#'
#' @docType data
#' @keywords datasets
#' @name turnout
#' @usage data(turnout)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::turnout)} variables and \Sexpr{nrow(SciencesPo::turnout)} observations.
NULL



#' @encoding UTF-8
#' @title Paired Data
#'
#' @description Artificial data for a paired experiment. This dataset contains the following columns:
#'
#' \itemize{
#' \item patient the patient id.
#' \item before_X before treatment.
#' \item after_Y after treatment.
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name paired
#' @usage data(paired)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::paired)} variables and \Sexpr{nrow(SciencesPo::paired)} observations.
NULL



#' @encoding UTF-8
#' @title Word frequencies from Mosteller and Wallace
#'
#' @description The data give the frequencies of words in works from four different sources: the political writings of eighteenth century American political figures Alexander Hamilton, James Madison, and John Jay, and the book Ulysses by twentieth century Irish writer James Joyce. This dataset contains the following columns:
#'
#' \itemize{
#' \item Hamilton Hamilton frequency.
#' \item HamiltonRank Hamilton rank.
#' \item Madison Madison frequency.
#' \item MadisonRank Madison rank.
#' \item Jay Jay frequency.
#' \item JayRank Jay rank.
#' \item Ulysses Word frequency in Ulysses.
#' \item UlyssesRank Word rank in Ulysses.
#'  }
#'
#'  @references Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.

#' @source Mosteller, F. and Wallace, D. (1964). Inference and Disputed Authorship: The Federalist. Reading, MA: Addison-Wesley.
#'
#' @docType data
#' @keywords datasets
#' @name words
#' @usage data(words)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::words)} variables and \Sexpr{nrow(SciencesPo::words)} observations.
NULL




#' @encoding UTF-8
#' @title Cathedrals
#'
#' @description Heights and lengths of Gothic and Romanesque cathedrals. This dataset contains the following columns:
#'
#' \itemize{
#' \item Type Romanesque or Gothic.
#' \item Height Total height, feet.
#' \item Length Total length, feet.
#'  }
#'
#'  @references
#' Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.
#'
#' @docType data
#' @keywords datasets
#' @name cathedrals
#' @usage data(cathedrals)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::cathedrals)} variables and \Sexpr{nrow(SciencesPo::cathedrals)} observations.
NULL



#' @encoding UTF-8
#' @title Burt's twin data
#'
#' @description The given data are IQ scores from identical twins; one raised in a foster home, and the other raised by birth parents.
#' This dataset contains the following columns:
#'
#' \itemize{
#' \item C Social class, C1=high, C2=medium, C3=low, a factor.
#' \item IQb biological.
#' \item IQf foster.
#'  }
#'
#'  @references
#'  Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.

#' @source Burt, C. (1966). The genetic estimation of differences in intelligence:
#'  A study of monozygotic twins reared together and apart. Br. J. Psych., 57, 147-153.
#'
#' @docType data
#' @keywords datasets
#' @name twins
#' @usage data(twins)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::twins)} variables and \Sexpr{nrow(SciencesPo::twins)} observations.
NULL




#' @encoding UTF-8
#' @title The Measure of US Presidents
#'
#' @description The US presidents and their main opponents' heights (cm). This dataset contains the following columns:
#'
#' \itemize{
#' \item election The election year.
#' \item winner The winner candidate.
#' \item winner.height The winner candidate's height in centimeters (cm).
#' \item winner.vote The popular vote support for the winner.
#' \item winner.party The winner's party.
#' \item opponent The main opponent candidate.
#' \item opponent.height The opponent candidate's height in centimeters (cm).
#' \item opponent.vote Popular vote support for the main opponent candidate.
#' \item opponent.party The opponent's party.
#' \item turnout The electorate turnout in percentages.
#' \item winner.bmi The winner Body Mass Index (BMI) estimate \code{(BMI = weight in kg/(height in meter)**2)}
#'  }
#' @source
#' US Presidents: \url{http://www.jimwegryn.com/Names/Presidents.php}
#' Inside Gov. \url{http://www.us-presidents.insidegov.com}.
#' Wikipedia: \url{http://en.wikipedia.org/wiki/United_States_presidential_election,_2012}.
#' Wikipedia:\url{http://en.wikipedia.org/wiki/Heights_of_presidents_and_presidential_candidates_of_the_United_States}.
#'
#' @docType data
#' @keywords datasets
#' @name Presidents
#' @usage data(Presidents)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::Presidents)} variables and \Sexpr{nrow(SciencesPo::Presidents)} observations.
NULL


#' @title Approval Ratings for President George W. Bush
#'
#' @description Approval ratings for George W. Bush.
#'
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::Bush)} variables and \Sexpr{nrow(SciencesPo::Bush)} observations.
#' \itemize{
#' \item start.date. Start date of the survey.
#' \item end.date. End date of the survey.
#' \item approve. Percent which approve of the president.
#' \item disapprove. Percent which disapprove of the president.
#' \item undecided. Percent undecided about the president.
#' }
#' @docType data
#' @keywords datasets
#' @name Bush
#' @usage data(Bush)
#'
NULL



#' @title Measurement System Units
#'
#' @description A dataset with measurement system units.
#'
#' \itemize{
#' \item from A character defining the original unit.
#' \item to A character defining the target unit.
#' \item factor The factor to be applied in conversion.
#' \item description Some details about the measure.
#' }
#' @docType data
#' @keywords datasets
#' @name Units
#' @usage data(Units)
#' @format A \code{data.frame} object with \Sexpr{ncol(SciencesPo::Units)} variables and \Sexpr{nrow(SciencesPo::Units)} observations.
#'
NULL

#' # setwd("~/SciencesPo/data")
#' # system("cp Bush.txt Bush-cp.txt")
#' # system("rm Bush-cp.txt.gz")
#' # system("gzip Bush.txt")
#' # system("rm Bush-cp.txt")



