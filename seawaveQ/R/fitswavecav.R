#' Function to prepare data and fit the seawaveQ model.
#' 
#' Fits the seawaveQ model (Vecchia and others, 2008) using a seasonal 
#' wave and continuous ancillary variables (streamflow anomalies and other 
#' continuous variables such as conductivity or sediment) to model water 
#' quality.
#' @name fitswavecav
#' @title Fit seasonal wave and continuous ancillary data for trend 
#' analysis
#' @note The assumed data format is one with columns for water-quality
#' concentration values and a related column for qualification of 
#' those values, such as in the case of left-censored values less 
#' than a particular value.  For example, a water-quality sample
#' was collected and the laboratory analysis indicated that the 
#' concentration was less than 0.01 micrograms per liter.  The 
#' USGS parameter code for simazine is 04035 (U.S. Geological Survey, 
#' 2013b).  When the data are retrieved through the National Water 
#' Information System: Web Interface 
#' (\url{http://waterdata.usgs.gov/nwis}; U.S. Geological Survey, 2013a), 
#' the concentration values are in a column labeled P04035 and the 
#' qualification information, or remark codes, are in a column labeled 
#' R04035.  To use this function, the argument pnames would be the unique 
#' identifier for simazine values and qualifications, 04035, and the 
#' qwcols argument would be c("R", "P") to indicate that the 
#' qualification column starts with an R and the values column starts with 
#' a P. \cr
#' Other users may have data in different format that can be 
#' changed to use with this function.  For example, a user may have
#' concentration values and qualification codes in one column, such
#' as a column labeled simzaine with the values 0.05, 0.10, <0.01, 
#' <0.01, and 0.90.  In this case, the less thans and any other 
#' qualification codes should be placed in a separate column.  The
#' column names for the qualification codes and the concentration values
#' should be the same with the exception of different beginning
#' letters to indicate which column is which.  The columns could be
#' named Rsimazine and Psimazine.  Then the argument pnames = "simazine" 
#' and the argument qwcols = c("R", "P"). \cr
#' Users should exercise caution when their water-quality data have 
#' multiple censoring limits and may want to recensor the data to a 
#' single censoring level.  Censoring and recensoring issues are discussed
#' in the text and Appendix 1 of Ryberg and others (2010).
#' @param cdat is the concentration data
#' @param cavdat is the continuous (daily) ancillary data
#' @param tanm is a character identifier that names the trend 
#' analysis run.  It is used to label output files.
#' @param pnames are the parameters (water-quality constituents) to 
#' analyze (omit the the starting character, for example for sulfate data 
#' indicated by P00945, enter "00945").  
#' @param yrstart is the starting year of the analysis (treated as January
#' 1 of that year).  Zero means the start date will be determined by the 
#' start date of cavdat, the continuous ancillary data.
#' @param yrend is the ending year of the analysis (treated as December 31
#' of that year).  Zero means the end date will be determined by the end 
#' date of cavdat, the continuous ancillary data.
#' @param tndbeg is the beginning (in whole or decimal years) of the 
#' trend period. Zero means the begin date will be the beginning of the
#' concentration data, cdat.
#' @param tndend is the end of the trend (treated as December 31
#' of that year). Zero means the end date will be the end of the 
#' concentration data, cdat.
#' @param iwcav is a character vector indicating which continuous
#' ancillary variables to include, if none are used for analysis,
#' use iwcav=c("none").
#' @param dcol is the column name for the dates, should be the same for 
#' both cdat and cavdat
#' @param qwcols is a character vector with the beginning of the
#' column headers for remarks code (default is R), and beginning of 
#' column headers for concentration data (default is P for parameter).
#' @param mclass has not been implemented yet but will provide
#' additional model options.
#' @keywords models regression ts survival
#' @return a pdf file containing plots of the data and modeled 
#' concentration, a text file containing a summary of the survival 
#' regression call for each model selected, and a list.  The first element
#' of the list is a data frame described under format.  The second element
#' of the list is the summary of the survival regression call.  The third 
#' element is the observed concentration data (censored and uncensored). 
#' The fourth element is the concentration data predicted by the model.  
#' The fifth element provides summary statistics for the predicted 
#' concentrations.
#' @format The data frame returned has one row for each parameter analyzed 
#' and the number of columns depend on the number of continuous ancillary
#' variables used. The general format is as follows: \cr
#' \tabular{lll}{
#'  pname \tab character \tab Parameter analyzed \cr
#'  mclass \tab numeric \tab Currently a value of 1 \cr
#'  jmod \tab numeric \tab The choice of pulse input function, an 
#'  integer 1--14. \cr
#'  hlife \tab numeric \tab the model half-life in months, an integer, 1 to 
#'  4 months \cr
#'  cmaxt \tab numeric \tab the decimal season of maximum concentration \cr
#'  scl \tab numeric \tab the scale factor from the 
#'  \code{survreg.object} \cr
#'  loglik \tab numeric \tab the log-likelihood for the model \cr
#'  cint \tab numeric \tab coefficient for model intercept \cr
#'  cwave \tab numeric \tab coefficient for the seasonal wave \cr
#'  ctnd \tab numeric \tab coefficient for the trend component of model \cr
#'  c[alphanumeric] \tab numeric \tab 0 or more coefficients for the 
#'  continuous ancillary variables\cr
#'  seint \tab numeric \tab standard error for the intercept \cr
#'  sewave \tab numeric \tab standard error for the seasonal wave \cr
#'  setnd \tab numeric \tab standard error for the trend \cr
#'  se[alphanumeric] \tab numeric \tab  0 or more standard errors for the 
#'  continuous ancillary variables\cr
#'  pvaltnd \tab numeric \tab the p-value for the trend line \cr
#' }
#' @seealso The functions that \code{fitswavecav} calls internally: \cr
#' \code{\link{prepData}} and \code{\link{fitMod}}.
#' @export
#' @author Aldo V. Vecchia and Karen R. Ryberg
#' @examples
#' data(swData)
#' modMoRivOmaha<-combineData(qwdat=qwMoRivOmaha, cqwdat=cqwMoRivOmaha)
#' myfit1 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha, 
#' tanm="myfit1", pnames=c("04035", "04037", "04041"), yrstart=1995, 
#' yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("flowa30","flowa1"), 
#' dcol="dates", qwcols=c("R","P"))
#' myfit2 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha, 
#' tanm="myfit2", pnames=c("04035", "04037", "04041"), yrstart=1995, 
#' yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("seda30","seda1"), 
#' dcol="dates", qwcols=c("R","P"))
#' myfit3 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha, 
#' tanm="myfit3", pnames=c("04035", "04037", "04041"), yrstart=1995, 
#' yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("flowa30","flowa1", 
#' "seda30", "seda1"), dcol="dates", qwcols=c("R","P"))
#' # trend model results
#' myfit3[[1]]
#' # example regression call
#' myfit3[[2]][[1]]
#' # first few lines of observed concentrations
#' head(myfit3[[3]])
#' # first few lines of predicted concentrations
#' head(myfit3[[4]])
#' # summary statistics for predicted concentrations
#' head(myfit3[[5]])
#' @references
#' Ryberg, K.R., Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2010, 
#' Trends in pesticide concentrations in urban streams in the United 
#' States, 1992--2008: U.S. Geological Survey Scientific Investigations 
#' Report 2010-5139, 101 p. (Also available at 
#' \url{http://pubs.usgs.gov/sir/2010/5139/}.)
#'
#' U.S. Geological Survey, 2013a, National Water Information System: 
#' Web Interface, accessed Febaruary 26, 2013, at
#' \url{http://waterdata.usgs.gov}.
#' 
#' U.S. Geological Survey, 2013b, Parameter code definition: National 
#' Water Information System: Web Interface, accessed Febaruary 26, 2013,
#' at \url{http://nwis.waterdata.usgs.gov/usa/nwis/pmcodes}.
#' 
#' Vecchia, A.V., Martin, J.D., and Gilliiom, R.J., 2008, Modeling 
#' variability and  trends in pesticide concentrations in streams: 
#' Journal of the American Water Resources Association, v. 44, no. 5, p. 
#' 1308-1324, \url{http://dx.doi.org/10.1111/j.1752-1688.2008.00225.x}.
fitswavecav <- function(cdat, cavdat, tanm="trend1", pnames, yrstart=0, 
                        yrend=0, tndbeg=0, tndend=0, iwcav=c("none"), 
                        dcol="dates", qwcols=c("R", "P"), mclass=1) {
  # perform data checks
							
  dtmes <- c("yrstart, yrend, tndbeg, tndend should all be numeric, 
            greater than or equal to 0.")
  if (!is.numeric(c(yrstart, yrend, tndbeg, tndend))) stop(dtmes)
  if ( yrstart < 0 | yrend < 0 | tndbeg < 0 | tndend < 0 ) stop(dtmes)
  if(yrstart > yrend) {yrstart <- 0; yrend <- 0}
  if(tndbeg > tndend) { tndbeg <- 0; tndend <- 0}
  
  # year function is from the lubridate package
  if(yrstart != 0) { 
    yrstart <- max(yrstart, year(min(cavdat[,dcol])))
  }
  if(yrstart == 0) { yrstart <- min(year(cavdat[, dcol])) }
  if(yrend != 0) { 
    yrend <- min(yrend, year(max(cavdat[, dcol]))) + 1
  }
  if(yrend == 0) { yrend <- year(max(cavdat[, dcol])) + 1 }
  if(tndbeg == 0) {tndbeg <- yrstart}
  if(tndend == 0) {tndend <- yrend}
  
  dtmsg <- paste("Trend begin year is ", tndbeg, "; trend end year is ", 
                 tndend, ".", sep="")
  message(dtmsg)
							
  npars <- length(pnames)
  nparsmes <- c("There are no parameters to analyze. User must pass
				at least one parameter name to the function using the pnames 
				argument.")
  if (npars < 1) stop(nparsmes)
  
  # prepare concentration data
  # prepare continous ancillary data
  myfun <- function(x) deparse(substitute(x)) 
  prepmsg <- paste("Preparing the data.")
  message(prepmsg)
  myData <- prepData(cdat, cavdat, yrstart, yrend, dcol, pnames, 
                     iwcav, qwcols) 
  cdat <- myData[[1]]
  cavdat <- myData[[2]]
  pnamesf<-vector(length=0)
  for  (iipar in (1:npars)) {
    if ( exists("stpars") ) { rm(stpars)}
    if ( exists("aovout") ) { rm(aovout)}
    if ( exists("obsDat") ) { rm(obsDat)}
    if ( exists("predDat") ) { rm(predDat)}
    if ( exists("predSummary") ) { rm(predSummary)}
    # for individual parameters, need to remove missing data
    # cdat columns for trend analysis of single parameter
    matches <- unique (grep(paste(iwcav, collapse="|"), names(cdat)))
    spcols <- c(1:4, grep(pnames[iipar], names(cdat)), matches )
    cdatiipar <- cdat[, spcols]
    cdatsub <- subset(cdatiipar, !is.na(cdatiipar[, paste(qwcols[2], pnames[iipar],
                                                sep="")]))
    cencol<-paste(qwcols[1], pnames[iipar], sep="")
    centmp <- cdatsub[, cencol]=='<'
    # check to see if at least 10 noncensored values
    if(sum(!centmp) > 9) {
      # fit model to data
      fitmsg <- paste("Fitting model for ", pnames[iipar], ".", sep="")
      message(fitmsg)
      myRes <- fitMod(cdatsub, cavdat, yrstart, yrend, tndbeg, tndend, 
                      tanm, pnames=pnames[iipar], qwcols, mclass=1)
      stpars <- myRes[[1]]
      aovout <- myRes[[2]]
      obsDat <- myRes[[3]][[1]]
      predDat <- myRes[[3]][[2]]
      mycol<-paste("P", pnames[iipar], sep="")
      predSummary <- data.frame(analysis=tanm, pname=pnames[iipar], 
                                predMeanConc=round(mean(predDat[,mycol]), 
                                               digits=5), 
                                matrix(round(quantile(predDat[, mycol], 
                                                probs=c(0.10, 0.25, 0.5, 0.75, 0.9)), 
                                             digits=5), 
                                       nrow=1, 
                                       dimnames=list(NULL,c("predQ10", "predQ25", 
                                                         "predQ50", "predQ75", "predQ90"))),
                                stringsAsFactors=FALSE)
      if ( !exists("aovoutall") ) {
        aovoutall <- myRes[[2]]
      } else { 
        aovoutall <- c(aovoutall, myRes[[2]])
      }
      pnamesf <- c(pnamesf,pnames[iipar])
    } else {
      notenoughmes <- paste("Less than 10 uncensored values for P",  
                            pnames[iipar], ".\n", 
                            "Analysis not performed.", sep="")
      message(notenoughmes)   
    }
    #  prepare output
    if ( exists("stpars") ) {
      stpars <- round(stpars, 5)
      row.names(stpars) <- NULL
      stparsout <- matrix(stpars[1,], nrow=1)
      if(iipar==1) {
        stparsoutall <- stparsout
      } else if (iipar > 1 & exists("stparsoutall") )  {
        stparsoutall <- rbind(stparsoutall, stparsout)
      } else {
        stparsoutall <- stparsout
      }
    }
    if ( exists("obsDat") ) {
      if(iipar==1) {
        obsdatall <- obsDat
      } else if (iipar > 1 & exists("obsdatall") )  {
        obsdatall <- merge(obsdatall, obsDat, all=TRUE)
      } else {
        obsdatall <- obsDat
      }
    }
    if ( exists("predDat") ) {
      if(iipar==1) {
        preddatall <- predDat
      } else if (iipar > 1 & exists("preddatall") )  {
        preddatall <- merge(preddatall, predDat, all=TRUE)
      } else {
        preddatall <- predDat
      }
    }
    if ( exists("predSummary") ) {
      if(iipar==1) {
        predSumAll <- predSummary
      } else if (iipar > 1 & exists("predSumAll") )  {
        predSumAll <- rbind(predSumAll, predSummary)
      } else {
        predSumAll <- predSummary
      }
    }
  }
	
  if ( exists("stparsoutall") ) {
    mod1 <- floor((stparsoutall[,2] - 1) / 4) + 1
    hlife1 <- stparsoutall[,2] - (mod1 - 1) * 4
    
    nxtmp <- length(stparsoutall[1,])
    stparsoutall <- cbind(mclass, mod1, hlife1, stparsoutall[,nxtmp], 
                        matrix(stparsoutall[, -c(1, 2, nxtmp)], 
                               nrow=dim(stparsoutall)[1]))
    stparsoutall <- data.frame(pnamesf, stparsoutall)
    if(iwcav[1] != 'none') { 
      names(stparsoutall) <- c('pname', 'mclass', 'jmod', 'hlife', 
                               'cmaxt', 'scl', 'loglik', 
                             paste('c', c('int', 'wave', 'tnd', iwcav), 
                                   sep=''), 
                             paste('se', c('int', 'wave', 'tnd', iwcav),
                                   sep=''), 'pvaltnd')
    } else if (iwcav[1] == 'none') {
      names(stparsoutall) <- c('pname', 'mclass', 'jmod', 'hlife', 
                               'cmaxt', 'scl', 'loglik',
                             paste('c', c('int', 'wave', 'tnd'), 
                                   sep=''),
                             paste('se', c('int', 'wave', 'tnd'), 
                                   sep=''), 'pvaltnd')
    }
    obsdatall$dectime<-round(obsdatall$dectime, digits=3)
    preddatall$dectime<-round(preddatall$dectime, digits=3)
    fitRes <- list(stparsoutall, aovoutall, obsdatall, preddatall, predSumAll)
    fitRes
  } else { message("No constituent had 10 or more uncensored values.")}
}
