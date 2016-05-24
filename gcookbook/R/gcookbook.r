#' gcookbook: Data sets for "R Graphics Cookbook"
#'
#' This package contains data sets used in the book "R Graphics Cookbook"
#' by Winston Chang, published by O'Reilly Media.
#'
#' @name gcookbook
#' @docType package
#' @aliases gcookbook package-gcookbook
NULL

#' Apple stock data
#'
#' Weekly stock data for AAPL (Apple, Inc.) from 1984 to 2012.
#'
#' @section Variables:
#' \itemize{
#'   \item date
#'   \item adj_price: Price, adjusted for splits and dividends.
#' }
#'
#' @docType data
#' @name aapl
#' @usage aapl
NULL


#' Homing in desert ants
#'
#' Data from an experiment on the homing performance of a desert ant,
#' \emph{Cataglyphis bicolor}.
#'
#' @section Variables:
#' \itemize{
#'   \item angle: Angle between true home direction and the direction that the
#'     ant went in (positive is clockwise).
#'   \item expt: Number of ants in the experimental condition that went in this
#'     direction.
#'   \item ctrl: Number of ants in the control condition that went in this
#'     direction.
#' }
#'
#' @docType data
#' @name anthoming
#' @usage anthoming
#'
#' @source
#'   Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. & Ostrowski, E. (1994),
#'   \emph{A Handbook of Small Data Sets}, Chapman & Hall.
#'
#'   Duelli, P. and Wehner, R. (1973), The spectal sensitivity of polarized
#'   light orientation in \emph{Cataglyphis bicolor} (Formicidae, Hymenoptera).
#'   \emph{Journal of Comparative Physiology}, \bold{86}, 36-53.
NULL


#' Summary of cabbages data set
#'
#' This data set has groupwise means, standard deviations, counts, and standard
#' error of the mean for the \code{\link{cabbages}} data set from the MASS
#' package. The purpose of this summarized data set is to make it easy to use
#' for example graphs.
#'
#' @docType data
#' @name cabbage_exp
#' @usage cabbage_exp
#'
#' @seealso The source data set in the MASS package, \code{\link{cabbages}}.
NULL


#' Global climate temperature anomaly data from 1800 to 2011 
#'
#' This data set includes estimated global temperature anamoly data for the
#' years 1800 through 2011. The anomaly is the difference from the baseline
#' temperature, which is the mean of the yearly temperatures from 1951-1980.
#'
#' @section Variables:
#' \itemize{
#'   \item Source: Data source (Berkeley, CRUTEM3, NASA).
#'   \item Year: Year for the estimate.
#'   \item Anomaly1y: Temperature anomaly in Celcius, smoothed over one year.
#'   \item Anomaly5y: Temperature anomaly in Celcius, smoothed over five years.
#'   \item Anomaly10y: Temperature anomaly in Celcius, smoothed over ten years.
#'   \item Unc10y: Uncertainty for 10-year-smoothed anomaly.
#' }
#'
#' @docType data
#' @name climate
#' @usage climate
#'
#' @source
#'   Berkeley Earth Project: \url{http://berkeleyearth.org/dataset/}
#'
#'   Climatic Research Unit (CRUTEM3): \url{http://www.cru.uea.ac.uk/cru/data/temperature/}
#'
#'   NASA: \url{http://data.giss.nasa.gov/gistemp/}
#'
NULL


#' Corneal thickness of eyes
#'
#' Corneal thickness of eight people who had glaucoma in one eye.
#'
#' @section Variables:
#' \itemize{
#'   \item affected Corneal thickness (in microns) of eye affected by glaucoma.
#'   \item notaffected Corneal thickness (in microns) of eye not affected by glaucoma.
#' }
#'
#' @docType data
#' @name corneas
#' @usage corneas
#'
#' @source
#'   Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. & Ostrowski, E. (1994),
#'   \emph{A Handbook of Small Data Sets}, Chapman & Hall.
#'
#'   Ehlers, N. On corneal thickness and introcular pressure, II. (1970).
#'   \emph{Acta Opthalmologica}, \bold{48}, 1107-1112.
NULL


#' Health and economic data about countries around the world from 1960-2010
#'
#' Health and economic data about countries around the world from 1960-2010,
#' from the World Bank.
#'
#' @section Variables:
#' \itemize{
#'   \item Name: Name of country
#'   \item Code: Short country code
#'   \item Year
#'   \item GDP: Per capita Gross Domestic Product, in adjusted 2011 U.S. Dollars
#'   \item laborrate: Labor rate.
#'   \item healthexp: Health expenditures in U.S. Dollars.
#'   \item infmortality: Infant mortality per 1000 live births.
#' }
#'
#' @docType data
#' @name countries
#' @usage countries
#'
#' @source
#'   World Bank: \url{http://data.worldbank.org/}
NULL


#' Convictions for drunkenness
#'
#' Number of people convicted for drunkenness at Tower Bridge and Lambeth
#' Magistrates' Courts from January 1 to June 27, 1970, classified by age
#' and sex.
#'
#' @docType data
#' @name drunk
#' @usage drunk
#'
#' @source
#'   Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. & Ostrowski, E. (1994),
#'   \emph{A Handbook of Small Data Sets}, Chapman & Hall.
#'
#'   Cook, T. (1971). \emph{New Society}, 20 May, 1971.
NULL


#' Height and weight of schoolchildren
#'
#' @section Variables:
#' \itemize{
#'   \item sex
#'   \item ageYear: Age in years.
#'   \item ageMonth: Age in months.
#'   \item heightIn: Height in inches.
#'   \item weightLb: Weight in pounds.
#' }
#'
#' @docType data
#' @name heightweight
#' @usage heightweight
#'
#' @source
#'   Lewis, T., & Taylor, L.R. (1967),
#'   \emph{Introduction to Experimental Ecology}, Academic Press.
NULL


#' Data from simulation of hurricane Isabel
#'
#' This data is from a simulation of hurricane Isabel in 2003. It includes
#' temperature and wind data for a 2139km (east-west) x 2004km (north-south) x 
#' 19.8km (vertical) volume. The simluation data is from the National Center
#' for Atmospheric Research, and it was used in the IEEE Visualization 2004
#' Contest.
#'
#' @section Variables:
#' \itemize{
#'   \item x: Latitude (x coordinate).
#'   \item y: Longitude (y coordinate).
#'   \item z: Height in km (z coordinate).
#'   \item vx: x wind component in m/s
#'   \item vy: y wind component in m/s
#'   \item vz: z wind component in m/s
#'   \item t: Temperature in Celcius
#'   \item speed: wind speed, sqrt(vx^2 + vy^2 + vz^2)
#' }
#'
#' @docType data
#' @name isabel
#' @usage isabel
#'
#' @source
#'   \url{http://vis.computer.org/vis2004contest/data.html}
#'
#'   \url{http://ncar.ucar.edu/}
#'
NULL


#' Successful sexual relations in Mad Men (TV show)
#'
#' Each row of this data frame represents a pair of characters who had a sexual
#' relationship on the TV show Mad Men, as of the end of season 4. This data
#' can be displayed with an undirected graph.
#'
#' @section Variables:
#' The placement of names in column Name1 as opposed to Name2 is arbitrary, and
#' not meaningful. In other words, for any row, you can swap the values of Name1
#' and Name2, and it will represent the same information.
#'
#' \itemize{
#'   \item Name1: Name of one sexual partner.
#'   \item Name2: Name of another sexual partner.
#' }
#'
#' @docType data
#' @name madmen
#' @usage madmen
#'
#' @seealso For a list of attempted sexual pairings, see \code{\link{madmen2}}.
#'
#' @source
#'   Wired Magazine 20.02, February 2012
NULL


#' Attempted sexual relations in Mad Men (TV show)
#'
#' Each row of this data frame represents a pair of characters on the TV show
#' Mad Men, as of the end of season 4. Each row represents an attempted sexual
#' relation: the character in the first column, Name1, attempted to have sex
#' with the character in the second column, Name2. If the relationship goes in
#' both directions (the characters had sex with each other), then there will be
#' two rows, representing each direction. This data can be displayed with a
#' directed graph.
#'
#' @section Variables:
#'
#' \itemize{
#'   \item Name1: Character who made sexual advances.
#'   \item Name2: Character who was the target of sexual advances.
#' }
#'
#' @docType data
#' @name madmen2
#' @usage madmen2
#'
#' @seealso For a list of successful sexual pairings, see \code{\link{madmen}}.
#'
#' @source
#'   Wired Magazine 20.02, February 2012
NULL


#' Marathon and half-marathon times
#'
#' This data set contains mrathon and half-marathon running times for 520
#' people. Each row represents one person's times.
#'
#' @section Variables:
#' \itemize{
#'   \item Half: Time in minutes, for half marathon.
#'   \item Full: Time in minutes, for full marathon.
#' }
#'
#' @docType data
#' @name marathon
#' @usage marathon
#'
#' @source
#'   Downey, A.B. (2011), \emph{Think Stats}, O'Reilly Media.
NULL


#' Means of results from an experiment on plant growth
#'
#' This data set simply has groupwise means of the \code{\link{PlantGrowth}}
#' data set. The purpose of this summarized data set is to make it easy to
#' use for example graphs.
#'
#' @docType data
#' @name pg_mean
#' @usage pg_mean
#'
#' @seealso The source data set, \code{\link{PlantGrowth}}.
#'
#' @source
#'   Dobson, A. J. (1983)
#'     \emph{An Introduction to Statistical Modelling}, Chapman & Hall.
NULL


#' Plum root cuttings (long format)
#'
#' This is data from an experiment to investigate the effect of cutting length
#' and planting time on the survival of plum root cuttings.
#'
#' @section Variables:
#' \itemize{
#'   \item length: Cutting length.
#'   \item time: Planting time.
#'   \item survival: Survival status.
#'   \item count: Number of plants.
#' }
#'
#' @docType data
#' @name plum
#' @usage plum
#'
#' @seealso This data frame is in "long" format. See \code{\link{plum_wide}} for
#'   the same data in "wide" format.
#'
#' @source
#'   Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. & Ostrowski, E. (1994),
#'   \emph{A Handbook of Small Data Sets}, Chapman & Hall.
#'
#'  Bartlett, M.S. (1935), Contingency table interactions,
#'  \emph{Journal of the Royal Statistical Society Supplement}, \bold{2}, 248-252.
NULL


#' Plum root cuttings (wide format)
#'
#' This is data from an experiment to investigate the effect of cutting length
#' and planting time on the survival of plum root cuttings.
#'
#' @section Variables:
#' \itemize{
#'   \item length: Cutting length.
#'   \item time: Planting time.
#'   \item dead: Number of dead plants in this condition.
#'   \item alive: Number of alive plants in this condition.
#' }
#'
#' @docType data
#' @name plum_wide
#' @usage plum_wide
#'
#' @seealso This data frame is in "wide" format. See \code{\link{plum}} for
#'   the same data in "long" format.
#'
#' @source
#'   Hand, D.J., Daly, F., Lunn, A.D., McConway, K.J. & Ostrowski, E. (1994),
#'   \emph{A Handbook of Small Data Sets}, Chapman & Hall.
#'
#'  Bartlett, M.S. (1935), Contingency table interactions,
#'  \emph{Journal of the Royal Statistical Society Supplement}, \bold{2}, 248-252.
NULL


#' Simple example data set
#'
#' This data set is for examples of R graphics.
#'
#' @seealso This data frame is in "wide" format. See \code{\link{simpledat_long}}
#'   for the same data in "long" format.
#' 
#' @docType data
#' @name simpledat
#' @usage simpledat
NULL


#' Simple example data set (long format)
#'
#' This data set is for examples of R graphics.
#'
#' @seealso This data frame is in "long" format. See \code{\link{simpledat}}
#'   for the same data in "wide" format.
#' 
#' @docType data
#' @name simpledat_long
#' @usage simpledat_long
NULL


#' Batting averages of the top hitters in Major League Baseball in 2001
#'
#' Batting statistics for the top 144 hitters in Major League Baseball in 2001.
#'
#' @section Variables:
#' Variables:
#' \itemize{
#'  \item id: Unique player id
#'  \item first: First name
#'  \item last: Last name
#'  \item name: Full name (first and last)
#'  \item year: Year of data
#'  \item stint
#'  \item team: Abbreviation of team played for
#'  \item lg: League (American League or National League)
#'  \item g: Number of games
#'  \item ab: Number of times at bat
#'  \item r: Number of runs
#'  \item h: Number of hits (times reached base because of a batted, fair ball 
#'    without error by the defense)
#'  \item 2b: Hits on which the batter reached second base safely
#'  \item 3b: Hits on which the batter reached third base safely
#'  \item hr: Number of home runs
#'  \item rbi: Runs batted in
#'  \item sb: Stolen bases
#'  \item cs: Caught stealing
#'  \item bb: Base on balls (walk)
#'  \item so: Strike outs
#'  \item ibb: Intentional base on balls
#'  \item hbp: Hits by pitch
#'  \item sh: Sacrifice hits
#'  \item sf: Sacrifice flies
#'  \item gidp: Ground into double-play
#'  \item avg: Batting average (hits divided by at-bats)
#' }
#'
#' @docType data
#' @name tophitters2001
#' @usage tophitters2001
#'
#' @source \url{http://www.baseball-databank.org/}.
NULL


#' Age distribution of population in the United States, 1900-2002
#'
#' These are the estimated (not counted) values by the U.S. Census.
#'
#' @section Variables:
#' \itemize{
#'   \item Year
#'   \item AgeGroup
#'   \item Thousands: Number of people, in thousands.
#' }
#'
#' @docType data
#' @name uspopage
#' @usage uspopage
#'
#' @source
#'   U.S. Census Bureau, Statistical Abstract of the United States, 2003, HS-3:
#'     \url{http://www.census.gov/statab/hist/HS-03.pdf} and
#'     \url{www.census.gov/statab/hist/02HS0003.xls}.
NULL


#' Change in population of states in the U.S. between 2000 and 2010
#'
#' This data set represents the percent change in population of states in the
#' U.S. from 2000 to 2010.
#'
#' @section Variables:
#' \itemize{
#'   \item State
#'   \item Abb: Abbreviated state name.
#'   \item Region: Region of country that the state is in.
#'   \item Change: Percent change in population.
#' }
#'
#' @docType data
#' @name uspopchange
#' @usage uspopchange
#'
#' @source
#'  U.S. Census Bureau, Statistical Abstract of the United States, 2012, Table 14.
#'    \url{http://www.census.gov/compendia/statab/2012/tables/12s0014.pdf}
NULL


#' Wind speed and direction
#'
#' This data set contains the wind speed and direction over the course of a
#' single day in Chicago (February 29, 2012).
#'
#' @section Variables:
#' \itemize{
#'   \item TimeUTC: Time of day in minutes; 0 is midnight.
#'   \item Temp: Temperature in Celcius.
#'   \item WindAvg: Average wind speed in m/s in this time block.
#'   \item WindMax: Maximum wind speed in m/s in this time block.
#'   \item WindDir: Average direction that wind comes from (0=north, 90=east).
#'   \item SpeedCat: Average wind speed, categorized in 5 m/s groups.
#'   \item DirCat: Average direction, categorized in 15-degree groups.
#' }
#'
#' @docType data
#' @name wind
#' @usage wind
#'
#' @source
#'   Great Lakes Environmental Research Laboratory:
#'     \url{http://www.glerl.noaa.gov/metdata/chi/}
NULL


#' World population estimates from 10,000 B.C. to 2,000 A.D.
#'
#' @section Variables:
#' \itemize{
#'   \item Year
#'   \item Population: Estimated population, in thousands
#' }
#'
#' @docType data
#' @name worldpop
#' @usage worldpop
NULL
