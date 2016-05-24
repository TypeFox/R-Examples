#' @name SAT
#' @title This dataset contains variables that address the
#' relationship between public school expenditures and academic
#' performance, as measured by the SAT.
#' @description SAT variables for public school and academic performance
#' @docType data
#' @usage SAT
#' @format A data frame with variables \code{state}, \code{expend}
#' (expenditure per pupil), \code{ratio} (pupil/teacher ratio);
#' \code{salary} (average teacher salary; \code{percentage of SAT
#' takers}; \code{verbal} (verbal score); \code{math} (math score);
#' \code{total} (average total).
#' @source from http://www.amstat.org/publications/jse/datasets/sat.txt
#' @references This data comes from
#' \url{http://www.amstat.org/publications/jse/secure/v7n2/datasets.guber.cfm}. It
#' is also included in the \pkg{mosaic} package and commented on at
#' \url{http://sas-and-r.blogspot.com/2012/02/example-920-visualizing-simpsons.html}. The
#' variables are described at
#' \url{http://www.amstat.org/publications/jse/datasets/sat.txt}.
#'
#' The author references the original source: The variables in this
#' dataset, all aggregated to the state level, were extracted from the
#' 1997 \emph{Digest of Education Statistics}, an annual publication
#' of the U.S. Department of Education.  Data from a number of
#' different tables were downloaded from the National Center for
#' Education Statistics (NCES) website (found at:
#' http://nces01.ed.gov/pubs/digest97/index.html) and merged
#' into a single data file.
NULL



##' @name Medicare
##' @title Sample from "Medicare Provider Charge Data"
##' @description Sample from "Medicare Provider Charge Data"
##' @docType data
##' @usage Medicare
##' @format A data frame with data about billings for procedures at many different hospitals
##' 
##' @source Retrieved from http://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data/index.html
##' 
##' @references This data came from
##' http://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data/index
##' and was referenced in the article
##' \url{http://www.nytimes.com/2013/05/08/business/hospital-billing-varies-wildly-us-data-shows.html},
##' as retrieved on 5/8/2013.
NULL


##' @name wellbeing
##' @title What makes us happy
##' @description What makes us happy
##' @docType data
##' @usage wellbeing
##' @format A data frame with data about what makes people happy (well being) along with several other covariates
##' 
##' @source Found from \url{http://prcweb.co.uk/lab/what-makes-us-happy/}.
##' 
##' @references \url{http://prcweb.co.uk/lab/what-makes-us-happy/} and \url{http://www.nationalaccountsofwellbeing.org/}
NULL


##' @name nisdc
##' @title A data frame measuring daily sea-ice extent from 1978 until 2013.
##' @description A data frame measuring daily sea-ice extent from 1978 until 2013.
##' @docType data
##' @usage nisdc
##' @format A data frame measuring daily sea-ice extent from 1978 until 2013.
##' 
##' @source Original data from
##' \url{ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/NH_seaice_extent_final.csv}
##' and
##' \url{ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/NH_seaice_extent_nrt.csv}. The
##' data can be downloaded with
##' \code{read.table(.,sep=",",skip=2,col.names=col.names)}.
##' Usage is like:
##' \code{ggplot(nisdc, aes(y=Extent,x=factor(Month))) + geom_boxplot() + facet_wrap(~Year)}
##' @references See the blog post \url{http://www.r-bloggers.com/arctic-sea-ice-at-lowest-levels-since-observations-began/} for a description and nice script to play with.
NULL


#' @docType data
#' @keywords datasets
#' @name wchomes
#' @usage data(wchomes)
#' @title A random sample of Wake County, North Carolina residential real estate plots
#' @description This data set comes from a JSE article
#' \url{http://www.amstat.org/publications/jse/v16n3/woodard.pdf} by
#' Roger Woodard. The data is described by: The information for this
#' data set was taken from a Wake County, North Carolina real estate
#' database. Wake County is home to the capital of North Carolina,
#' Raleigh, and to Cary. These cities are the fifteenth and eighth
#' fastest growing cities in the USA respectively, helping Wake County
#' become the ninth fastest growing county in the country. Wake County
#' boasts a 31.18% growth in population since 2000, with a population
#' of approximately 823,345 residents. This data includes 100 randomly
#' selected residential properties in the Wake County registry denoted
#' by their real estate ID number. For each selected property, 11
#' variables are recorded. These variables include year built, square
#' feet, adjusted land value, address, et al.
#' @format A data frame
#' @source \url{http://www.amstat.org/publications/jse/v16n3/woodard.xls}
NULL


#' @docType data
#' @keywords datasets
#' @name ObamaApproval
#' @usage data(ObamaApproval)
#' @title Approval ratings for President Obama
#' @description A collection of approval ratings for President Obama spanning a duration from early 2010 to the summer of 2013.
#' @format A data frame
#' @source Scraped on 7-5-13 from \url{http://www.realclearpolitics.com/epolls/other/president_obama_job_approval-1044.html}
NULL


#' @docType data
#' @keywords datasets
#' @name ceo2013
#' @usage data(ceo2013)
#' @title CEO compensation in 2013
#' @description Data on top 200 CEO compensations in the year 2013
#' @format A data frame
#' @source Scraped from \url{http://www.nytimes.com/interactive/2013/06/30/business/executive-compensation-tables.html?ref=business}
NULL

#' @name movie_data_2011
#' @title movie data for 2011 by weekend
#' @description A data frame with variables \code{Previous} (previous weekend rank), \code{Movie} (title), \code{Distributor}, \code{Genre}, \code{Gross} (per current weekend), \code{Change} (change from previous week), \code{Theaters} (number of theaters), \code{TotalGross} (total gross to date), \code{Days} (days out), \code{weekend} (weekend of report)
#' @docType data
#' @usage movie_data_2011
#' @format A data frame 
#' @source Scraped from pages such as \url{http://www.the-numbers.com/box-office-chart/weekend/2011/04/29}
NULL


#' @docType data
#' @keywords datasets
#' @name snacks
#' @usage data(snacks)
#' @title Snack data from the USDA
#' @description subset of SR26 data on nutrients compiled by the USDA.
#' @format A data frame with some nutrition variables
#' @source This data came from the SR26 data set found at \url{http://www.ars.usda.gov/Services/docs.htm?docid=8964}.
NULL





