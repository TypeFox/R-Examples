##' Monthly Mean Temperature and Total Precipitation for Forstenrieder
##' Park, Munich
##' 
##' This dataset gives the monthly mean temperature and total
##' precipitation at Forstenrieder Park, Munich, Bavaria, Germany.
##' @source Zang C, Pretzsch H, Rothe A (2012) Size-dependent
##' responses to summer drought in Scots pine, Norway spruce and
##' common oak. Trees - Structure and Function, 26, 557-569.
##' @docType data
##' @keywords datasets
##' @name muc_clim
##' @usage data(muc_clim)
##' @format A \code{data.frame} containing four columns with year,
##' month, temperature and precipitation.
NULL

##' Tree-Ring Chronology of a Spruce Population near Munich
##' 
##' This dataset gives the tree-ring indices for \emph{Picea abies} at
##' Forstenrieder Park, Munich, Bavaria, Germany. The chronology
##' represents 20 cores from 10 trees. The series were read in using
##' \code{\link[dplR]{read.rwl}} from package \code{dplR}, detrended
##' using a 60a spline with 50\% frequency cutoff (function
##' \code{\link[dplR]{detrend}}), and averaged to a chronology using a
##' robust mean \code{\link[dplR]{chron}}.
##' @source Zang C, Pretzsch H, Rothe A (2012) Size-dependent
##' responses to summer drought in Scots pine, Norway spruce and
##' common oak. Trees - Structure and Function, 26, 557-569.
##' @docType data
##' @keywords datasets
##' @name muc_spruce
##' @usage data(muc_spruce)
##' @format A \code{data.frame} containing tree-ring indices and
##' replication depth with respective years as \code{rownames}.
NULL

##' Modeled Tree-Ring Chronology of a Spruce Population near Munich
##' 
##' This dataset gives the modelled tree-ring widths for \emph{Picea
##' abies} at Forstenrieder Park, Munich, Bavaria, Germany. Tree
##' growth was modeled as a response of low temperatures in previous
##' and current July and August, high temperatures in current February
##' and March, and high precipitation amounts in current July and
##' August.
##' @source Zang C, Pretzsch H, Rothe A (2012) Size-dependent
##' responses to summer drought in Scots pine, Norway spruce and
##' common oak. Trees - Structure and Function, 26, 557-569.
##' @docType data
##' @keywords datasets
##' @name muc_fake
##' @usage data(muc_fake)
##' @format A \code{data.frame} containing tree-ring indices and
##' replication depth with respective years as \code{rownames}.
NULL

##' Monthly Precipitation Sums for Rothenburg ob der Tauber, Germany
##' 
##' This dataset gives the monthly precipitation sum at Rothenburg ob der Tauber,
##' Bavaria, Germany in a decadal (DENDROCLIM2002-style) format)
##' @source Zang C, Rothe A, Weis W, Pretzsch H (2011) Zur
##' Baumarteneignung bei Klimawandel: Ableitung der
##' Trockenstress-Anfaelligkeit wichtiger Waldbaumarten aus
##' Jahrringbreiten. Allgemeine Forst- und Jagdzeitung, 182, 98-112.
##' @docType data
##' @keywords datasets
##' @name rt_prec
##' @usage data(rt_prec)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of precipitation data in mm rainfall.
NULL

##' Monthly Temperature Means for Rothenburg ob der Tauber, Germany
##' 
##' This dataset gives the monthly temperature means at Rothenburg ob der Tauber,
##' Bavaria, Germany in a decadal (DENDROCLIM2002-style) format)
##' @source Zang C, Rothe A, Weis W, Pretzsch H (2011) Zur
##' Baumarteneignung bei Klimawandel: Ableitung der
##' Trockenstress-Anfaelligkeit wichtiger Waldbaumarten aus
##' Jahrringbreiten. Allgemeine Forst- und Jagdzeitung, 182, 98-112.
##' @docType data
##' @keywords datasets
##' @name rt_temp
##' @usage data(rt_temp)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of temperature data in degree Celsius.
NULL

##' Tree-Ring Chronology of a Spruce Population at Rothenburg ob der Tauber
##' 
##' This dataset gives the tree-ring indices for \emph{Picea abies} at Rothenburg
##' ob der Tauber, Bavaria, Germany. The chronology represents 20 cores from 10
##' trees. The series were read in using \code{\link[dplR]{read.rwl}} from
##' package \code{dplR}, detrended using a 60a spline with 50\% frequency cutoff
##' (function \code{\link[dplR]{detrend}}), and averaged to a chronology using a 
##' robust mean \code{\link[dplR]{chron}}.
##' @source Zang C, Rothe A, Weis W, Pretzsch H (2011) Zur
##' Baumarteneignung bei Klimawandel: Ableitung der
##' Trockenstress-Anfaelligkeit wichtiger Waldbaumarten aus
##' Jahrringbreiten. Allgemeine Forst- und Jagdzeitung, 182, 98-112.
##' @docType data
##' @keywords datasets
##' @name rt_spruce
##' @usage data(rt_spruce)
##' @format A \code{data.frame} containing tree-ring indices and
##' replication depth with respective years as \code{rownames}.
NULL

##' Tree-Ring Chronology of a Pine Population at Penota, Spain
##'
##' This dataset gives the tree-ring indices for \emph{Pinus
##' sylvestris} at a site near Penota, Spain. The chronology
##' represents 18 dated series and spans the period between 1763 and
##' 1991. The chronology was obtained from the ITRDB (see source), and
##' read in using \code{\link[dplR]{read.crn}} from package
##' \code{dplR}.
##' @source Fernandez-Cancio A, Genova Fuster M: Fernandez-Cancio -
##' Penota - PISY - ITRDB SPAI020,
##' http://hurricane.ncdc.noaa.gov/pls/paleox/f?p=519:1:::::P1_STUDY_ID:3255,
##' accessed 2014/07/03.
##' @docType data
##' @keywords datasets
##' @name spai020
##' @usage data(spai020)
##' @format A \code{data.frame} containing tree-ring indices and
##' replication depth with respective years as \code{rownames}.
NULL

##' Monthly Temperature Means for Spain
##' 
##' This data set gives the monthly temperature means at country level
##' for Spain in a decadal (DENDROCLIM2002-style) format from the TYN
##' CY 1.1 data set.
##' @source Mitchell TD, Carter TR, Jones PD, Hulme M, New M, et
##' al. (2004) A comprehensive set of high-resolution grids of monthly
##' climate for Europe and the globe: the observed record (1901-2000)
##' and 16 scenarios (2001-2100). Tyndall Centre for Climate Change
##' Research Working Paper 55, 25.
##' @docType data
##' @keywords datasets
##' @name spain_temp
##' @usage data(spain_temp)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of temperature data in degree Celsius.
NULL

##' Monthly Precipitation Sums for Spain
##' 
##' This data set gives the monthly precipitation sums at country
##' level for Spain in a decadal (DENDROCLIM2002-style) format from
##' the TYN CY 1.1 data set.
##' @source Mitchell TD, Carter TR, Jones PD, Hulme M, New M, et
##' al. (2004) A comprehensive set of high-resolution grids of monthly
##' climate for Europe and the globe: the observed record (1901-2000)
##' and 16 scenarios (2001-2100). Tyndall Centre for Climate Change
##' Research Working Paper 55, 25.
##' @docType data
##' @keywords datasets
##' @name spain_prec
##' @usage data(spain_prec)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of temperature data in degree Celsius.
NULL


##' Monthly Temperature Means for Norway
##' 
##' This data set gives the monthly temperature means at country level
##' for Norway in a decadal (DENDROCLIM2002-style) format from the TYN
##' CY 1.1 data set.
##' @source Mitchell TD, Carter TR, Jones PD, Hulme M, New M, et
##' al. (2004) A comprehensive set of high-resolution grids of monthly
##' climate for Europe and the globe: the observed record (1901-2000)
##' and 16 scenarios (2001-2100). Tyndall Centre for Climate Change
##' Research Working Paper 55, 25.
##' @docType data
##' @keywords datasets
##' @name norway_temp
##' @usage data(norway_temp)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of temperature data in degree Celsius.
NULL

##' Monthly Precipitation Sums for Norway
##' 
##' This data set gives the monthly precipitation sums at country
##' level for Norway in a decadal (DENDROCLIM2002-style) format from
##' the TYN CY 1.1 data set.
##' @source Mitchell TD, Carter TR, Jones PD, Hulme M, New M, et
##' al. (2004) A comprehensive set of high-resolution grids of monthly
##' climate for Europe and the globe: the observed record (1901-2000)
##' and 16 scenarios (2001-2100). Tyndall Centre for Climate Change
##' Research Working Paper 55, 25.
##' @docType data
##' @keywords datasets
##' @name norway_prec
##' @usage data(norway_prec)
##' @format A \code{data.frame} containing thirteen columns with year
##' and twelve months of temperature data in degree Celsius.
NULL

##' Tree-Ring Chronology of a Pine Population at Visdalen, Norway
##'
##' This dataset gives the tree-ring indices for \emph{Pinus
##' sylvestris} at a site near Visdalen, Norway. The chronology
##' represents 34 dated series and spans the period between 1600 and
##' 1983. The raw measurements were obtained from the ITRDB (see
##' source), and read in, detrended with cubic splines (frequency
##' cutoff of 50 percent at two thirds of the series length) and
##' averaged using \code{\link[dplR]{read.rwl}},
##' \code{\link[dplR]{detrend}}, and \code{\link[dplR]{chron}} from
##' package \code{dplR}.
##' @source Briffa K. Briffa - Visdalen - PISY - ITRDB
##' NORW015. http://hurricane.ncdc.noaa.gov/pls/paleox/f?p=519:1:::::P1_STUDY_ID:2861. Accessed
##' 2014/07/03.
##' @docType data
##' @keywords datasets
##' @name norw015
##' @usage data(norw015)
##' @format A \code{data.frame} containing tree-ring indices and
##' replication depth with respective years as \code{rownames}.
NULL
