#' Function to calculate short-, medium-, and long-term hydrologic anomalies
#'
#' This function was written with streamflow data in mind because streamflow 
#' is the most commonly used exogenous variable for trend models for water 
#' quality;  however, the function is generic so that users may experiment with 
#' anomalies from other daily hydrologic data.  Examples of the inclusion of 
#' streamflow anomalies in trend analysis of nutrients, pesticides and 
#' surface water can be found in Alexander and Smith (2006), Ryberg and Vecchia 
#' (2006), Ryberg and others (2010), Sullivan and others (2009), Vecchia
#' (2003), Vecchia (2005), and Vecchia and others (2008).
#'
#' @name compAnom
#' @title Calculates  anomalies
#' @param dataset is the daily hydrologic data returned from \link{importDVs} 
#' or data otherwise obtained and in the same format as that produced by 
#' \link{importDVs}.
#' @param which indicates which set of anomalies; 1 calculates the 1-year, 
#' 30-day, and 1-day anomalies; 2 calculates the 100-day, 10-day, and 1-day 
#' anomalies; 3 calculates the 30-day and 1-day anomalies; and 4 calculates 
#' the 10-year, 5-year, 1-year, one-quarter-year (seasonal), and 1-day 
#' anomalies.
#' @keywords ts multivariate
#' @return A list.  In the cases of "which" equal to 1 or 2, the 
#' first element of the list is a data frame containing the station 
#' identification number, dates, streamflow, and long-term, mid-term, and 
#' short-term anomalies.  The next three elements of the list are the length 
#' in days of the long-term, mid-term, and short-term streamflow anomalies.  
#' In the case of "which" equal to 3, the first element of the list is a data 
#' frame containing the station identification number, dates, streamflow, and 
#' mid-term and short-term anomalies.  The next two elements of the list are 
#' the length in days of the mid-term and short-term streamflow anomalies.  
#' In the case of "which" equal to 4, the first element of the list is a data 
#' frame containing the station identification number, dates, streamflow, and 
#' 10-year, 5-year, annual, seasonal, and daily streamflow anomalies.  The 
#' next five elements of the list are the length in days of the 10-year, 
#' 5-year, annual, seasonal, and daily streamflow anomalies.  
#' @export
#' @examples
#' q05054000.85 <- importDVs("05054000", sdate="1985-10-01", edate="2010-09-30")
#' anoms05054000.1 <- compAnom(q05054000.85, which=1)
#' anoms05054000.2 <- compAnom(q05054000.85, which=2)
#' anoms05054000.3 <- compAnom(q05054000.85, which=3)
#' anoms05054000.4 <- compAnom(q05054000.85, which=4)
#' @references
#' Alexander, R.B. and Smith, R.A., 2006, Trends in the nutrient enrichment of 
#' U.S. rivers during the late 20th century and their relation to changes in
#' probable stream trophic conditions: Limnology and Oceanography, v. 51, no.
#' 1, Part 2: Eutrophication of Freshwater and Marine Ecosystems, p. 
#' 639--654, accessed August 1, 2012, at 
#' \url{http://www.jstor.org/stable/4499617}.
#' 
#' Ryberg, K.R. and Vecchia, A.V., 2006, Water-quality trend analysis and 
#' sampling design for the Devils Lake Basin, North Dakota, January 1965 
#' through September 2003: U.S. Geological Survey Scientific Investigations 
#' Report 2006--5238, 64 p., accessed August 1, 2012, at 
#' \url{http://pubs.usgs.gov/sir/2006/5238/}.
#'
#' Ryberg, K.R. and Vecchia, A.V., 2012, waterData---An R packge for retrieval, 
#' analysis, and anomaly calculation of daily hydrologic time series data, 
#' version 1.0: U.S. Geological Survey Open-File Report 2012--1168, 8 p.  
#' (Also available at \url{http://pubs.usgs.gov/of/2012/1168/}.)
#' 
#' Ryberg, K.R., Vecchia, A.V., Martin, J.D., Gilliom, R.J., 2010, Trends in 
#' pesticide concentrations in urban streams in the United States, 1992-2008: 
#' U.S. Geological Survey Scientific Investigations Report 2010--5139, 101 p. 
#' (Also available at \url{http://pubs.usgs.gov/sir/2010/5139/}.)
#'
#' Sullivan, D.J., Vecchia, A.V., Lorenz, D.L., Gilliom, R.J., Martin, J.D., 
#' 2009, Trends in pesticide concentrations in corn-belt streams, 1996--2006: 
#' U.S. Geological Survey Scientific Investigations Report 2009--5132, 75 p., 
#' accessed Ausugst , 2012, at \url{http://pubs.usgs.gov/sir/2009/5132/}.
#'
#' Vecchia, A.V., 2003, Relation between climate variability and stream water 
#' quality in the continental United States, Hydrological Science and 
#' Technology, v. 19, no. 1, 77--98.
#'
#' Vecchia, A.V., 2003, Water-quality trend analysis and sampling design for 
#' streams in North Dakota, 1971--2000: U.S. Geological Survey Scientific 
#' Investigations Report 2003--4094, 73 p., accessed August 1, 2012, at 
#' \url{http://nd.water.usgs.gov/pubs/wri/wri034094/index.html}.
#'
#' Vecchia, A.V., 2005, Water-quality trend analysis and sampling design for 
#' streams in the Red River of the North Basin, Minnesota, North Dakota, and 
#' South Dakota, 1970--2001: U.S. Geological Survey Scientific Investigations 
#' Report 2005--5224, 54 p. accessed August 1, 2012, at 
#' \url{http://pubs.usgs.gov/sir/2005/5224/}.
compAnom <- function(dataset,which=1) {
  if (which==1) {
    lt<-365
    mt<-30
    st<-1
    n1 <- length(dataset[,1])
    if ( n1 >= lt ) {
      qdf <- dataset[,c("staid","dates","val")]
      qlog10 <- log10(dataset$val)
      xmean <- mean(qlog10,na.rm=TRUE)  
      xlt <- rep(NA,n1)
      for (i in (lt:n1)) { 
        xlt[i] <- mean(qlog10[(i - (lt - 1)):i],na.rm=FALSE)
      }
      ltfa <- xlt-xmean
      xmt <- rep(NA,n1)
      for (i in (mt:n1)) { 
        xmt[i] <- mean(qlog10[(i - (mt - 1)):i],na.rm=FALSE)
      }
      mtfa <- xmt-xlt
      # short term anomaly is based on 1 day
      stfa <- qlog10 - xmt
  
      qdf <- cbind(qdf,ltfa,mtfa,stfa)
      list(qdf, lt,mt,st)
    } else {
      stop("Dataset not long enough to calculate 1-year anomaly, 
           try which=2 or which=3")
    }
  }
  else if (which==2) {
    lt<-100
    mt<-10
    st<-1
    n1 <- length(dataset[,1])
    qdf <- dataset[,c("staid","dates","val")]
    qlog10 <- log10(dataset$val)
    xmean <- mean(qlog10,na.rm=TRUE)  
    xlt <- rep(NA,n1)
    for (i in (lt:n1)) { 
      xlt[i] <- mean(qlog10[(i - (lt - 1)):i],na.rm=FALSE)
    }
    ltfa <- xlt-xmean
    xmt <- rep(NA,n1)
    for (i in (mt:n1)) { 
      xmt[i] <- mean(qlog10[(i - (mt - 1)):i],na.rm=FALSE)
    }
    mtfa <- xmt-xlt
    # short term anomaly is based on 1 day
    stfa <- qlog10 - xmt
    
    qdf <- cbind(qdf,ltfa,mtfa,stfa)
    list(qdf,lt,mt,st)
  }
  else if (which==3) {
    st <- 1
    mt <- 30
    n1 <- length(dataset[,1])
    qdf <- dataset[,c("staid","dates","val")]
    qlog10 <- log10(dataset$val)
    xmean <- mean(qlog10,na.rm=TRUE)  
    xmt <- rep(NA,n1)
    for (i in (mt:n1)) { 
      xmt[i] <- mean(qlog10[(i - (mt - 1)):i],na.rm=FALSE)
    }
    mtfa <- xmt-xmean
    # short term anomaly is based on 1 day
    stfa <- qlog10 - xmt
  
    qdf <- cbind(qdf,mtfa,stfa)
    list(qdf, mt, st)
  }
  else if (which==4) {
    ta<-10*365
    fa<-5*365
    lt<-365
    mt<-90
    st<-1
    n1 <- length(dataset[,1])
    if ( n1 >= ta ) {
      qdf <- dataset[,c("staid","dates","val")]
      qlog10 <- log10(dataset$val)
      xmean <- mean(qlog10,na.rm=TRUE)
      xta <- rep(NA,n1)
      for (i in (ta:n1)) { 
        xta[i] <- mean(qlog10[(i - (ta - 1)):i],na.rm=FALSE)
      }
      tfa <- xta-xmean
    
      xfa <- rep(NA,n1)
      for (i in (fa:n1)) { 
        xfa[i] <- mean(qlog10[(i - (fa - 1)):i],na.rm=FALSE)
      }    
      ffa <- xfa-xta
    
      xlt <- rep(NA,n1)
      for (i in (lt:n1)) { 
        xlt[i] <- mean(qlog10[(i - (lt - 1)):i], na.rm=FALSE)
      }    
      ltfa <- xlt-xfa
    
      xmt <- rep(NA,n1)
      for (i in (mt:n1)) { 
        xmt[i] <- mean(qlog10[(i - (mt - 1)):i],na.rm=FALSE)
      }
      mtfa <- xmt-xlt
    
      # short term anomaly is based on 1 day
      stfa <- qlog10 - xmt
  
      qdf <- cbind(qdf,tfa,ffa,ltfa,mtfa,stfa)
      dimnames(qdf)[[2]] <- c("staid", "dates", "val", "tenyranom", "fiveyranom",
                            "annualanom","seasanom","dailyanom")
      list(qdf, ta, fa, lt, mt, st)
    } else { 
    stop("Dataset not long enough to calculate 10-year anomaly, try which=1")
    }
  }
  else {
    stop("which must be a numeric value 1, 2, 3, or 4")
  }
}
