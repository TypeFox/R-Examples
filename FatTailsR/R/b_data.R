

#' @include a_FatTailsR-package.R



#' @title Get DS Dataset
#'
#' @description 
#' A function to extract the log-returns 
#' of 16 financial series and time series provided by the packages \code{datasets} 
#' (EuStockMarkets, sunspot.year) and \code{timeSeries} (USDCHF, MSFT, LPP2005REC).
#' The 16 datasets are converted to a list of numeric without any reference 
#' to the original dates. This list is usually called \code{DS}, hence the name.
#' 
#' @details
#' The dataset is usually created by the instruction \code{DS <- getDSdata()}.
#' Then, it is used with a call to DS[[j]] with j in 1:16. 
#' \enumerate{
#'   \item{ "USDCHF" (USDCHF, timeSeries) }
#'   \item{ "MSFT" (MSFT, timeSeries) }
#'   \item{ "DAX" (EuStockMarkets, datasets) }
#'   \item{ "SMI" (EuStockMarkets, datasets) }
#'   \item{ "CAC" (EuStockMarkets, datasets) }
#'   \item{ "FTSE" (EuStockMarkets, datasets) }
#'   \item{ "SBI" (LPP2005REC, timeSeries) }
#'   \item{ "SPI" (LPP2005REC, timeSeries) }
#'   \item{ "SII" (LPP2005REC, timeSeries) }
#'   \item{ "LMI" (LPP2005REC, timeSeries) }
#'   \item{ "MPI" (LPP2005REC, timeSeries) }
#'   \item{ "ALT" (LPP2005REC, timeSeries) }
#'   \item{ "LPP25" (LPP2005REC, timeSeries) }
#'   \item{ "LPP40" (LPP2005REC, timeSeries) }
#'   \item{ "LPP60" (LPP2005REC, timeSeries) }
#'   \item{ "sunspot" (sunspot.year, datasets) }
#' }
#' Note that \code{sunspot.year} is regularly updated with each new version of  
#' \code{R}. The generated dataset is \code{logreturn(sunspot.year + 1000)}.
#' 
#' @seealso    
#' \code{\link[datasets]{EuStockMarkets}}, \code{\link[datasets]{sunspot.year}}, 
#' \code{\link[timeSeries]{TimeSeriesData}}, \code{\link{regkienerLX}}, 
#' \code{\link{fitkienerX}}
#' 
#' @examples    
#' 
#' require(timeSeries) 
#' 
#' getDSdata
#' DS  <- getDSdata()
#' attributes(DS)
#' sapply(DS, length)
#' sapply(DS, head)
#' 
#' @export
#' @name getDSdata
getDSdata <- function() {
DSenv <- new.env(parent = baseenv())
utils::data("USDCHF", "MSFT", "LPP2005REC", package = "timeSeries", envir = DSenv)
utils::data("EuStockMarkets", "sunspot.year", package = "datasets", envir = DSenv)
prices2returns <- function(x) { 100*diff(log(as.numeric(x))) }
DS <- list(
	"USDCHF"	= prices2returns(DSenv$USDCHF),
	"MSFT"		= prices2returns(DSenv$MSFT[,4]),
	"DAX"		= prices2returns(DSenv$EuStockMarkets[,1]),
	"SMI"		= prices2returns(DSenv$EuStockMarkets[,2]),
	"CAC"		= prices2returns(DSenv$EuStockMarkets[,3]),
	"FTSE"		= prices2returns(DSenv$EuStockMarkets[,4]),
	"SBI"		= 100*as.numeric(DSenv$LPP2005REC[,1]),
	"SPI"		= 100*as.numeric(DSenv$LPP2005REC[,2]),
	"SII"		= 100*as.numeric(DSenv$LPP2005REC[,3]),
	"LMI"		= 100*as.numeric(DSenv$LPP2005REC[,4]),
	"MPI"		= 100*as.numeric(DSenv$LPP2005REC[,5]),
	"ALT"		= 100*as.numeric(DSenv$LPP2005REC[,6]),
	"LPP25"		= 100*as.numeric(DSenv$LPP2005REC[,7]),
	"LPP40"		= 100*as.numeric(DSenv$LPP2005REC[,8]),
	"LPP60"		= 100*as.numeric(DSenv$LPP2005REC[,9]),
	"sunspot"	= prices2returns(DSenv$sunspot.year+1000) )
return(DS)
}



#' @title Datasets dfData, tData, xData, zData, extractData : dfData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the data.frame format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name dfData
NULL

#' @title Datasets dfData, tData, xData, zData, extractData : tData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the timeSeries format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name tData
NULL

#' @title Datasets dfData, tData, xData, zData, extractData : xData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the xts format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name xData
NULL

#' @title Datasets dfData, tData, xData, zData, extractData : zData
#'
#' @description
#' A list of datasets in data.frame, timeSeries, xts and zoo formats. 
#' This is the zoo format. 
#' Visit \code{\link{extractData}} for more information. 
#' 
#' @keywords datasets
#' @docType data
#' @name zData
NULL



#' @title Datasets dfData, tData, xData, zData, extractData : extractData
#'
#' @description 
#' dfData, tData, xData, zData are datasets made of lists of data.frame, timeSeries, 
#' xts and zoo components. They describe dates and prices of 10 financial series
#' used in the documents and demos presented at 8th and 9th R/Rmetrics conferences  
#' (2014, 2015). See the references. 
#' The last serie (CHF, interest rates in Switzerland) exhibits negative prices. 
#' The (log)returns must be calculated. All distributions of returns exhibit fat tails.
#' Function \code{extractData} converts subsets of tData, xData, zData from a list 
#' format to a 2 dimensions format.
#' 
#' @param    pr	   	 character. Extract prices or returns: \code{c("p","r","prices","returns")}. 
#' @param    ft    	 character. Output format among \code{c("tss","xts","zoo","dfr","bfr","mat")}. 
#' @param    start   character. Start date.
#' @param    end	 character. End date.
## ## @param    ttData	 dataset. The dataset used by extractData. Do not change it.
#'
#' @details
#' 10 financial series presented in four different formats for convenient use: 
#' dfData is the complete dataset in data.frame format with dates as row.names 
#' and one or several columns. 
#' tData, xData and zData are respectively of timeSeries, xts and zoo formats  
#' and display one column only per stock. 
#' \enumerate{
#'   \item{ "GOLD" from 1999-01-04 to 2013-12-31, dim 3694x1 (df, t, x, z). }
#'   \item{ "DEXIA" from 2008-10-27 to 20009-10-26, dim 255x1 (t, x, z), 255x5 (df). }
#'   \item{ "SOCGEN" from 1992-07-20 to 2013-12-31, dim 5445x1 (t, x, z), 5445x5 (df). }
#'   \item{ "VIVENDI" from 1992-07-20 to 2013-12-31, dim 5444x1 (t, x, z) 5444x5 (df). }
#'   \item{ "EURUSD" from 1999-01-03 to 2013-12-31, dim 3843x1 (t, x, z), 3843x4 (df). }
#'   \item{ "VIX" from 2004-01-02 to 2013-12-31, dim 2517x1 (t, x, z), 2517x4 (df). }
#'   \item{ "CAC40" from 1988-01-04 to 2013-12-31, dim 6574x1 (t, x, z), 6574x4 (df). }
#'   \item{ "DJIA" from 1896-05-26 to 2013-12-31, dim 32064x1 (df, t, x, z). }
#'   \item{ "SP500" from 1957-01-02 to 2013-12-31, dim 14350x1 (df, t, x, z). }
#'   \item{ "CHF" from 1995-01-02 to 2013-09-13, dim 4880x1 (t, x, z), 4880x8 (df). 
#'             Interst rates in Switzerland. Include negative prices at the end 
#'             of the dataset. Care is required to calculate the returns! }
#' }
#' 
#' Function \code{extractData} extracts 8 financiel series (DEXIA and CHF are not converted) 
#' and converts them into a 2 dimensions object with any of the following class:
#' \itemize{
#'   \item{ "tss" is for timeSeries. }
#'   \item{ "xts" is for xts. }
#'   \item{ "zoo" is for zoo. }
#'   \item{ "dfr" is the usual R data.frame with the Date as index. }
#'   \item{ "bfr" is a data.frame with the Date in the first column. }
#'   \item{ "mat" is a matrix with the Date in the rownames. }
#' }
#' The \code{start} date must be posterior to "2007-01-01" (default) and the 
#' \code{end} date must be anterior to "2013-12-31" (default).
#'  
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from: 
#' \url{http://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#' 
#' @seealso 
#' \code{\link{tData}}, \code{\link{xData}}, \code{\link{zData}}, \code{\link{dfData}},
#' \code{\link{getDSdata}}.
#' 
#' @examples    
#' 
#' require(zoo) 
#' require(xts) 
#' require(timeSeries) 
#' 
#' ### dfData, tData, xData, zData : prices only
#' attributes(dfData); attributes(tData); attributes(xData); attributes(zData) 
#' lapply(dfData, head, 3)
#' lapply( tData, head, 3)
#' lapply( xData, head, 3)
#' lapply( zData, head, 3)
#' 
#' ### extractData : prices and logreturns
#' head(ptD <- extractData("p", "tss", "2009-01-01", "2012-12-31")) ; tail(ptD)
#' head(rtD <- extractData("r", "tss")) 
#' head(pxD <- extractData("p", "xts")) 
#' head(rxD <- extractData("r", "xts")) 
#' head(pzD <- extractData("p", "zoo")) 
#' head(rzD <- extractData("r", "zoo")) 
#' head(pbD <- extractData("p", "bfr")) 
#' head(rbD <- extractData("r", "bfr")) 
#' head(pmD <- extractData("p", "mat")) 
#' head(rmD <- extractData("r", "mat")) 
#' 
#' @export
#' @name extractData
extractData <- function(pr = "p", ft = "tss",
						start = "2007-01-01", end = "2013-12-31") {
# tDataEnv <- new.env(parent = baseenv())
# utils::data("tData", package = character(0)) # v1.3-10 no visible binding
if (is.na(charmatch(strtrim(pr, 1)[1], 
          c("p", "r")))) { 
		  stop("pr must be p or r, prices or returns.") 
		  }
if (is.na(charmatch(ft, 
          c("tss","xts","zoo","dfr","bfr","mat")))) { 
		  stop("ft must be tss, xts, zoo, dfr, bfr, mat.") 
		  }
	start3   <- as.Date("2006-12-31")
	end3     <- as.Date("2014-01-01")
	start2   <- as.Date(start)
	end2     <- as.Date(end)
	startend <- c(start3, start2, end2, end3)
if (checkquantiles(startend)) {
		start1   <- start2
		end1     <- end2
	} else {
		start1   <- start3
		end1     <- end3
		warning("start or end dates are not in range 2007-01-01 - 2013-12-31. Ignored")
	}
	# 24/10/2015 10 heures. il faut remplacer tData par 
	# get("tData") ou get("tData", pos = parent.frame(n = 1)) 
	lData  <- lapply(get("tData"), timeSeries::window, start = start1, end = end1) 
	ptD    <- cbind(lData[[1]], lData[[3]], lData[[4]], lData[[5]], 
				    lData[[6]], lData[[7]], lData[[8]], lData[[9]])
	colnames(ptD) <- names(lData)[-c(2, 10)]
	rtD <- 100 * diff(log(ptD))
	z   <- if (strtrim(pr, 1)[1] == "p") {
				switch(ft,
				"tss" = ptD,
				"xts" = xts::as.xts(ptD),
				"zoo" = zoo::as.zoo(ptD),
				"dfr" = as.data.frame(ptD, stringsAsFactors = FALSE),
				"bfr" = data.frame(Date = as.Date(rownames(ptD)), ptD, 
						     row.names = NULL, stringsAsFactors = FALSE),
				"mat" = as.matrix(ptD))
			} else {
				switch(ft,
				"tss" = rtD,
				"xts" = xts::as.xts(rtD),
				"zoo" = zoo::as.zoo(rtD),
				"dfr" = as.data.frame(rtD, stringsAsFactors = FALSE),
				"bfr" = data.frame(Date = as.Date(rownames(rtD)), rtD, 
						     row.names = NULL, stringsAsFactors = FALSE),
				"mat" = as.matrix(rtD))
			}
	# dimnames(z) <- list( "DATES" = dimnames(z)[[1]], 
                        # "STOCKS" = dimnames(z)[[2]])
return(z)
}






#' @title Several Vectors of Probabilities
#' 
#' @description
#' Several vectors of probabilities used in FatTailsR. 
#' Remark: pprobs5 <- sort(c(pprobs2, pprobs3, pprobs4)).
#' 
#' pprobs0 <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
#' 
#' pprobs1 <- c(0.01, 0.05, 0.95, 0.99)
#' 
#' pprobs2 <- c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99)
#' 
#' pprobs3 <- c(0.001, 0.0025, 0.005, 0.995, 0.9975, 0.999)
#' 
#' pprobs4 <- c(0.0001, 0.00025, 0.0005, 0.9995, 0.99975, 0.9999)
#'
#' pprobs5 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999)
#'
#' pprobs6 <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.50, 
#' 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999) 
#'
#' pprobs7 <- c(0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99) 
#' 
#' pprobs8 <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999) 
#' 
#' pprobs9 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#' 0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
#' 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 
#' 
#' 
#' @seealso 
#' The conversion function \code{\link{getnamesk}}
#'  
#' @export
#' @name pprobs0
pprobs0 <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99) 
#' @export
#' @rdname pprobs0
pprobs1 <- c(0.01, 0.05, 0.95, 0.99) 
#' @export
#' @rdname pprobs0
pprobs2 <- c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99) 
#' @export
#' @rdname pprobs0
pprobs3 <- c(0.001, 0.0025, 0.005, 0.995, 0.9975, 0.999)
#' @export
#' @rdname pprobs0
pprobs4 <- c(0.0001, 0.00025, 0.0005, 0.9995, 0.99975, 0.9999)
#' @export
#' @rdname pprobs0
pprobs5 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 
#' @export
#' @rdname pprobs0
pprobs6 <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.50, 
             0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999) 
#' @export
#' @rdname pprobs0
pprobs7 <- c(0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
             0.95, 0.975, 0.99)
#' @export
#' @rdname pprobs0
pprobs8 <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999) 
#' @export
#' @rdname pprobs0
pprobs9 <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
             0.10, 0.17, 0.25, 0.33, 0.41, 0.50, 0.59, 0.67, 0.75, 0.83, 0.90, 
			 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999) 



#' @title Parameter Subsets
#' 
#' @description
#' Some vectors of parameter names to be used with parameter \code{exfitk} in 
#' functions regkienerLX(.., exfitk = ...) and \code{fitkienerX(.., exfitk = ...)}
#' or to subset the vector (or matrix) \code{fitk } obtained after regression 
#' \code{fitk <- regkienerLX(..)$fitk} or estimation \code{fitk <- fitkienerX(..)}. 
#' Visit \code{\link{fitkienerX}} for details on each parameter.
#' 
#' \code{exfit0 <- c("lh", "ret")}
#' 
#' \code{exfit1 <- c("m", "g", "a", "k", "w", "d", "e")}
#' 
#' \code{exfit2 <- c("m1", "sd", "sk", "ke", "m1x", "sdx", "skx", "kex")}
#' 
#' \code{exfit3 <- c("q.01", "q.05", "q.95", "q.99", "ltm.025", "rtm.975")}
#' 
#' \code{exfit4 <- c("VaR.01", "VaR.05", "VaR.95", "VaR.99", "ES.025", "ES.975")}
#' 
#' \code{exfit5 <- c("c.01", "c.05", "c.95", "c.99", "h.025", "h.975")}
#' 
#' \code{exfit6 <- c(exfit1, exfit2, exfit3, exfit4, exfit5)}
#' 
#' \code{exfit7 <- c(exfit0, exfit1, exfit2, exfit3, exfit4, exfit5)}
#' 
#' 
#' @examples     
#' 
#' require(minpack.lm)
#' require(timeSeries)
#' 
#' ### Load the datasets and select one number j in 1:16
#' j      <- 5
#' DS     <- getDSdata()
#' (fitk  <- regkienerLX(DS[[j]])$fitk)
#' fitk[exfit3]
#' fitkienerX(DS[[j]], exfitk = exfit3)
#' 
#' 
#' @export
#' @name exfit0
exfit0 <- c("lh", "ret") 
#' @export
#' @rdname exfit0
exfit1 <- c("m", "g", "a", "k", "w", "d", "e")
#' @export
#' @rdname exfit0
exfit2 <- c("m1", "sd", "sk", "ke", "m1x", "sdx", "skx", "kex")
#' @export
#' @rdname exfit0
exfit3 <- c("q.01", "q.05", "q.95", "q.99", "ltm.025", "rtm.975")
#' @export
#' @rdname exfit0
exfit4 <- c("VaR.01", "VaR.05", "VaR.95", "VaR.99", "ES.025", "ES.975")
#' @export
#' @rdname exfit0
exfit5 <- c("c.01", "c.05", "c.95", "c.99", "h.025", "h.975")
#' @export
#' @rdname exfit0
exfit6 <- c(exfit1, exfit2, exfit3, exfit4, exfit5)
#' @export
#' @rdname exfit0
exfit7 <- c(exfit0, exfit1, exfit2, exfit3, exfit4, exfit5)


