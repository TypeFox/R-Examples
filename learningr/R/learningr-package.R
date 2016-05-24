#' Datasets to accompany the book Learning R
#'
#' \code{learningr} contains datasets that are used in examples in the book 
#' ``Learning R''.
#'
#' @docType package
#' @name learningr
#' @aliases learningr learningr-package
#' @references
#'   Cotton, R (2013) \emph{Learning R}
#'   O Reilly.  ISBN 978-1-4493-5710-8.
#' @author Richard Cotton \email{richierocks@@gmail.com}
NULL

#' Alpe d'Huez
#' 
#' Fastest times for the Alpe d'Huez stage of the Tour de France cycle race, 
#' plus some contextual information on year and drug use allegations.
#' 
#' @note
#' The data was kindly compiled by William Hogan.
#'
#' 
#' The dataset is not guaranteed to be error free.  Please double check the data
#' if you want to use it for something serious.
#' 
#' @docType data
#' @name alpe_d_huez
#' @aliases alpe_d_huez2
#' @format A data frame with the following columns.
#' \describe{
#' \item{Time}{Character time of ride in the form M' S".} 
#' \item{NumericTime}{Numeric time of ride in minutes.} 
#' \item{Name}{Name of rider.} 
#' \item{Year}{Year of race.} 
#' \item{Nationality}{Nationality of rider.} 
#' \item{DrugUse}{Have allegations of drug use been made against the
#' rider?  In \code{alpe_d_huez} the values are "Y" or "N"; in
#' \code{alpe_d_huez2} this is a logical vector.} 
#' }
#' @references William Hogan \email{william.hogan@@amd.com} compiled the data
#' from 
#' \url{http://en.wikipedia.org/wiki/Alpe_d\%27Huez}.
#' Richard Cotton \email{richierocks@@gmail.com} made some modifications 
#' while importing it into R.
NULL

#' Crab tag
#' 
#' Depth and temperature data from a sensor tag attached to an edible crab 
#' (Cancer Pangurus) in the North Sea in 2008 and 2009.
#'
#' 
#' @note
#' This data was kindly supplied by Ewan Hunter of CEFAS. 
#' It is part of a larger dataset consisting of many crabs.
#' 
#' @docType data
#' @name crab_tag
#' @format A list with 5 elements, as follows.
#' \code{id_block} is a list with 2 elements
#' \describe{
#' \item{Firmware Version No}{Version of the firmware used in the crab tag.} 
#' \item{Firmware Build Level}{Build level of the firmware used in the crab tag.} 
#' }
#' \code{tag_notebook} is a list with 5 elements.
#' \describe{
#' \item{Mission Day}{Number of days of data.} 
#' \item{Last Deployment Date}{Date and time that the tag was released into the sea.} 
#' \item{Deployed by Host Version}{UNKNOWN} 
#' \item{Downloaded by Host Version}{UNKNOWN} 
#' \item{Last Clock Set Date}{UNKNOWN} 
#' }
#' 
#' \code{lifetime_notebook} is a list with 3 elements.
#' \describe{
#' \item{Tag ID}{The unique ID of the tag.} 
#' \item{Pressure Range}{UNKNOWN} 
#' \item{No of sensors}{The number of sensors on the tag.} 
#' }
#' 
#' \code{deployment_notebook} is a data frame with X columns.
#' \describe{
#' \item{Start}{Start date and time of recording.} 
#' \item{Stop}{Stop date and time of recording.} 
#' \item{Logging Rate}{UNKNOWN} 
#' \item{Pointer}{UNKNOWN} 
#' \item{PA inc}{UNKNOWN} 
#' \item{sensors}{UNKNOWN} 
#' \item{Flags}{UNKNOWN} 
#' \item{Resolution}{UNKNOWN} 
#' \item{Fast Rate}{UNKNOWN} 
#' \item{V1}{UNKNOWN} 
#' \item{V2}{UNKNOWN} 
#' }
#' 
#' \code{daylog} is a data frame with X columns.
#' \describe{
#' \item{Mission.Day}{Integer number of days since the start of recording.} 
#' \item{Date}{Date of record.} 
#' \item{Max.Temp}{Maximum temperature (Celcius) recorded that day.} 
#' \item{Min.Temp}{Minimum temperature (Celcius) recorded that day.} 
#' \item{Max.Depth}{Maximum depth (m) recorded that day.} 
#' \item{Min.Depth}{Minimum depth (m) recorded that day.} 
#' \item{Batt.Volts}{Voltage of tag battery.} 
#' }
#' @references Ewan Hunter \email{ewan.hunter@@cefas.gov.uk} ran the 
#' project where the data was collected. A full analysis is in
#' Hunter, E, Eaton, D, Stewart, C, Lawler, A & Smith, M. 2013.
#' Edible crabs "go west": migrations and incubation cycle revealed by 
#' electronic tags. \url{https://www.ncbi.nlm.nih.gov/pubmed/23734180}. 
#' Richard Cotton \email{richierocks@@gmail.com} made some modifications 
#' while importing it into R.
NULL

#' Deer Endocranial Volume
#' 
#' The dataset contains the endocranial volume of 33 red deer  (Cervus 
#' elaphus), using four different methods of measurement: cardiac
#' tomography, filling the skull with glass beads (yes, the skulls are
#' from dead deer), simply measuring the length, width and height and 
#' multiplying the numbers, and using Finarelli's equation.  "Endocranial
#' volume" is a proxy for brain size.
#' 
#' @note
#' The data was kindly provided by Corina Logan.
#' Second measurements are provided for several of the deer.
#' Finarelli's equation is used for estimating the brain volume of 
#' non-bovid ruminant Artiodactylid species (say that 10 times fast).
#' \deqn{\ln(volume(skull)) = 2.6616 * \ln(width(skull)) - 6.2722}{%
#'       \ln(volume(skull)) = 2.6616 * \ln(width(skull)) - 6.2722}
#' 
#' 
#' @docType data
#' @name deer_endocranial_volume
#' @format A data frame with the following columns.
#' \describe{
#' \item{SkullID}{A unique identifier for each red deer.} 
#' \item{VolCT}{Endocranial volume calculated by cardiac tomography.} 
#' \item{VolBead}{Endocranial volume calculated by glass beads.} 
#' \item{VolLWH}{Endocranial volume calculated by length*width*height.} 
#' \item{VolFinarelli}{Endocranial volume calculated by Finarelli's equation.} 
#' \item{VolCT2}{A second measurement via cardiac tomography.} 
#' \item{VolBead2}{A second measurement via glass beads.} 
#' \item{VolLWH2}{A second measurement via l*w*h.} 
#' }
#' @references The dataset was collected by Corina Logan 
#' \email{itsme@@corinalogan.com}.  It is stored in the Dryad Digital  
#' Repository, \url{doi:10.5061/dryad.4t7h2}.   
#' A more full analysis is given in the paper
#' Logan CJ, Clutton-Brock TH. 2013. Validating methods for estimating
#' endocranial volume in individual red deer (Cervus elaphus). Behavioural
#' Processes 92:143-146. doi:10.1016/j.beproc.2012.10.015.
#' \url{http://www.sciencedirect.com/science/article/pii/S037663571200232X}
#' Richard Cotton \email{richierocks@@gmail.com} made some modifications 
#' while importing it into R.
NULL

#' English Monarchs
#' 
#' Names, dates and houses of English kings and queens from post-Roman rule
#' (the fifth century) until England invaded Ireland in the Early 13th century.
#' 
#' @note
#' This dataset is a bit messy and ambiguous in places, because history
#' is like that. In fact, the messy parts of the dataset are in general
#' a good indicator that soemthing interesting was happening at the time.
#' (See, for example, missing or multiple rulers, starts and ends of 
#' reigns in the same year, and rulers that appear several times with 
#' different territories.)  Evan defining a monarch of England is tricky.
#' Most of the monarchs in this dataset were around before England 
#' existed (it consisted of seven territories called the heptarchy).
#' The data stops before John I (the bad guy from the Robin Hood stories)
#' because he proclaimed himself King of Ireland, although some people
#' consider monarchs up to Anne, five hundred years later, to be English 
#' Monarchs even though they ruled over Ireland, Wales and Scotland to
#' varying degrees.
#' 
#' The heptarchy consisted of East Anglia, Essex, Kent, Mercia, 
#' Northumbria, Sussex and Wessex.  Northubria was originally divided into
#' Deria and Bernicia.  There are also periods of Norse and Danish rule.
#'
#' The dataset was compiled from Wikipedia and thus is not guaranteed to 
#' be error free.  Please double check the data if you want to use it for 
#' something serious.
#' 
#' @docType data
#' @name english_monarchs
#' @format A data frame with the following columns.
#' \describe{
#' \item{name}{Name of monarch(s).} 
#' \item{house}{Royal house of monarch(s).} 
#' \item{start.of.reign}{Year they rose to power.} 
#' \item{end.of.reign}{Year they left power.} 
#' \item{domain}{Region of England ruled over.} 
#' }
#' @references Richard Cotton \email{richierocks@@gmail.com} compiled
#' the dataset from various Wikipedia pages. 
#' \url{http://en.wikipedia.org/wiki/Kings_of_england}
#' \url{http://en.wikipedia.org/wiki/Kings_of_East_Anglia}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Essex}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Kent}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Mercia}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Northumbria}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Sussex}
#' \url{http://en.wikipedia.org/wiki/Kings_of_Wessex}
NULL

#' Gonorrhoea
#' 
#' Rates of gonorrhoea infection in the US by year, with contextual 
#' information about age, ethnicity and gender.
#' @docType data
#' @name gonorrhoea
#' @format A data frame with the following columns.
#' \describe{
#' \item{Year}{Year that infected people visited the clinic.} 
#' \item{Age.Group}{Age group of person infected.} 
#' \item{Ethnicity}{Ethnicity of person infected.} 
#' \item{Gender}{Gender of person infected.} 
#' \item{Rate}{Number of infections per 100000 people.} 
#' }
#' @references Compiled by Richard Cotton \email{richierocks@@gmail.com} 
#' from \url{http://www.cdc.gov/std/stats11/tables/22b.htm} 
NULL

#' Hafu
#' 
#' Half-caste manga characters.
#'
#' 
#' @note
#' The dataset was kindly provided by Gwern Branwen.
#'
#' \code{hafu2} is a lightly cleaned up version of \code{hafu}.
#' 
#' Gwern's notes: The following list includes manga, light novel, anime, 
#' and video game characters (there being little point in keeping the 
#' mediums separate). It also includes characters who are not hafu  
#' themselves but a quarter-foreign inasmuch as they imply a hafu at some  
#' point. Characters are treated separately even if they are in the same  
#' work (eg. siblings). Classification is based on in-universe or out-of- 
#' universe information, since appearance can be highly misleading in anime 
#' (blue eyes may indicate heroic status, rather than being Caucasian; hair  
#' color may be chosen for contrast against other characters or signal 
#' stereotypes like red hair indicating a fiery personality), and different 
#' groups will identify the same anime character as belonging to their own 
#' race (Lu 2009), perhaps due to minimalistic drawings intended to save  
#' money or enable viewers to project themselves onto a character. 
#' 
#' @docType data
#' @name hafu
#' @aliases hafu2
#' @format Both data frames have the following columns.
#' \describe{
#' \item{Year}{Integer year that the manga was made.} 
#' \item{Series}{Name of series.} 
#' \item{Character}{Name of character.} 
#' \item{Gender}{Gender of character.} 
#' \item{Father}{Nationality of character's father.} 
#' \item{Mother}{Nationality of character's mother.} 
#' \item{Eyes}{Character's eye colour.} 
#' \item{Hair}{Character's hair colour.} 
#' \item{Notes}{Notes on data collection or ambiguity.} 
#' }
#'
#' \code{hafu2} has these additional columns.
#'
#' \describe{
#' \item{FathersNationalityIsUncertain}{} 
#' \item{MothersNationalityIsUncertain}{} 
#' }
#' @references The dataset was compiled by Gwern Branwen
#' \email{gwern0@@gmail.com}.  The original is available from 
#' \url{http://www.gwern.net/hafu#list}.
NULL

#' Obama vs. McCain
#' 
#' State-by-state voting information in the 2008 US presidential election, 
#' along with contextual information on income, unemployment, ethnicity and
#' religion.
#' 
#' @note
#' Religious identification data are not available for Alaska and Hawaii. The 
#' totals of these columns is generally less than 100, since some people didn't
#' give an answer.
#' The District of Columbia is included, even though it isn't a state.
#' The dataset is not guaranteed to be error free.  Please double check the data 
#' if you want to use it for something serious.
#' 
#' @docType data
#' @name obama_vs_mccain
#' @format A data frame with 52 observations (one for each US state)
#' and the following columns.
#' \describe{
#' \item{State}{The name of the US state.}
#' \item{Region}{The US Federal region.} 
#' \item{Obama}{Percentage of voters who voted for Barack Obama in the 2008 presidential election.}
#' \item{McCain}{Percentage of voters who voted for John McCain in the 2008 presidential election.}
#' \item{Turnout}{Percentage of people who voted in the 2008 presidential election.}
#' \item{Unemployment}{Percentage of people who are unemployed.}
#' \item{Income}{Mean annual income in US dollars.}
#' \item{Population}{Number of people living in the state.}
#' \item{Catholic}{Percentage of people identifying as Catholic.}
#' \item{Protestant}{Percentage of people identifying as Protestant.}
#' \item{Other}{Percentage of people identifying as religious, but not Catholic or Protestant.}
#' \item{Non.religious}{Percentage of people identifying as non-religious.} 
#' \item{Black}{Percentage of people identifying as black.}
#' \item{Latino}{Percentage of people identifying as Latino.}
#' \item{Urbanization}{Percentage of people living in an urban area.}
#' }
#' @references This dataset was kindly compiled and provided by Edwin 
#' Thoen \email{edwinthoen@@hotmail.com}.
#' 
#' The voting information came from 
#' \url{http://www.uselectionatlas.org/}, extracted on 2011-12-09.
#' 
#' The ethnicity, income and urbanisation information came from 
#' \url{http://quickfacts.census.gov}, extracted on 2011-12-09.
#' 
#' The unemployment information came from 
#' \url{http://data.bls.gov/timeseries/LNS14000000}, extracted 2011-12-09.
#' 
#' The religious information came from Table 12 of the American Religious
#' Identification Survey 2008. 
#' \url{http://commons.trincoll.edu/aris/files/2011/08/ARIS_Report_2008.pdf}.
NULL

