#' An R package for retrieval, analysis, and anomaly calculation of daily 
#' hydrologic time deries data.
#'
#' This package imports U.S. Geological Survey (USGS) daily hydrologic data 
#' from USGS web services, plots the data, addresses some common data problems, 
#' and calculates and plots anomalies.  For a description of anomalies see 
#' Vecchia (2003), and for examples of the application of streamflow anomalies 
#' in trend analysis of nutrients, pesticides and surface water, see Alexander 
#' and Smith (2006), Ryberg and Vecchia (2006), Ryberg and others (2010), 
#' Sullivan and others (2009), Vecchia (2005), and Vecchia and others (2008).
#'
#' \tabular{ll}{
#' Package: \tab waterData\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.4\cr
#' Date: \tab 2014-11-18\cr
#' License: \tab Unlimited | file LICENSE \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name waterData-package
#' @aliases waterData
#' @docType package
#' @title Hydrologic Data Retrieval and Analysis and Anomaly Calculation
#' @author Karen R. Ryberg \email{kryberg@@usgs.gov} and 
#' Aldo V. Vecchia \email{avecchia@@usgs.gov}
#' @keywords package
#' @references 
#' Alexander, R.B. and Smith, R.A., 2006, Trends in the nutrient enrichment of 
#' U.S. rivers during the late 20th century and their relation to changes in
#' probable stream trophic conditions: Limnology and Oceanography, v. 51, no.
#' 1, Part 2: Eutrophication of Freshwater and Marine Ecosystems, p. 
#' 639--654., accessed August 1, 2012, at \url{http://www.jstor.org/stable/4499617}.
#'
#' Ryberg, K.R. and Vecchia, A.V., 2006, Water-quality trend analysis and 
#' sampling design for the Devils Lake Basin, North Dakota, January 1965 
#' through September 2003: U.S. Geological Survey Scientific Investigations 
#' Report 2006--5238, 64 p., accessed August 1, 2012, at
#'  \url{http://pubs.usgs.gov/sir/2006/5238/}.
#'
#' Ryberg, K.R. and Vecchia, A.V., 2012, waterData---An R package for retrieval, 
#' analysis, and anomaly calculation of daily hydrologic time series data, 
#' version 1.0: U.S. Geological Survey Open-File Report 2012--1168, 8 p. 
#' (Also available at \url{http://pubs.usgs.gov/of/2012/1168/}.)
#' 
#' Ryberg, K.R., Vecchia, A.V., Martin, J.D., Gilliom, R.J., 2010, Trends in 
#' pesticide concentrations in urban streams in the United States, 1992--2008: 
#' U.S. Geological Survey Scientific Investigations Report 2010-5139, 101 p., 
#' accessed August 1, 2012, at \url{http://pubs.usgs.gov/sir/2010/5139/}.
#'
#' Sullivan, D.J., Vecchia, A.V., Lorenz, D.L., Gilliom, R.J., Martin, J.D., 
#' 2009, Trends in pesticide concentrations in corn-belt streams, 1996--2006: 
#' U.S. Geological Survey Scientific Investigations Report 2009-5132, 75 p., 
#' accessed August 1, 2012, at \url{http://pubs.usgs.gov/sir/2009/5132/}.
#'
#' Vecchia, A.V., 2003, Relation between climate variability and stream water 
#' quality in the continental United States, Hydrological Science and 
#' Technology, v. 19 no. 1, 77--98.
#'
#' Vecchia, A.V., 2003, Water-quality trend analysis and sampling design for 
#' streams in North Dakota, 1971--2000: U.S. Geological Survey Scientific 
#' Investigations Report 2003--4094, 73 p., accessed August 1, 2012, at  
#' \url{http://nd.water.usgs.gov/pubs/wri/wri034094/index.html}.
#'
#' Vecchia, A.V., 2005, Water-quality trend analysis and sampling design for 
#' streams in the Red River of the North Basin, Minnesota, North Dakota, and 
#' South Dakota, 1970--2001: U.S. Geological Survey Scientific Investigations 
#' Report 2005--5224, 54 p., accessed August 1, 2012, at 
#' \url{http://pubs.usgs.gov/sir/2005/5224/}.
NULL
