#' Function to create boxplots of water-quality data that include
#' censored values.
#'
#' This function determines the columns within the data set that have
#' concentration data, based on them having column headings that start
#' with P (or a user-specificed indicater in the second element of the 
#' qwcols argument) followed by a number.  The function determines the 
#' associated remark, or qualification columns, based on them having 
#' column headings that start with R (or a user-specificed indicater in 
#' the first element of the qwcols argument) followed by numbers that 
#' match the associated concentration data. Then it determines which 
#' values are censored, indicated by a less than symbol in the R columns, 
#' performs regression on order statistics, \link{ros}, using the NADA 
#' package, and estimates values for the censored concentrations for 
#' constituents with less than 90-percent censoring.  The water-quality 
#' concentrations are then depicted by boxplots.
#' @name rosBoxPlot
#' @title Boxplot of water-quality data
#' @param data is the dataset with columns that begin with P
#' followed by a number indicating concentration data and 
#' columns that begin with R followed by numbers that match those
#' of the concentration data indicating qualification codes.  See 
#' example data sets for more information about the data format, 
#' \link{IllRivValleyCty} and \link{qwMoRivOmaha}.
#' @param qwcols is a character vector with the beginning of the
#' column headers for remarks code (default is R), and beginning of 
#' column headers for concentration data (default is P for parameter).
#' @param site is a label for the plot title indicating the site where
#' the water-quality samples were collected.
#' @param ... arguments to be passed to \link{plot} method.
#' @return a boxplot
#' @note The regression on order statistics function in R package NADA 
#' (Lee, 2012), \code{\link{ros}}, is an implementation of a regression on 
#' order statistics designed for multiply-censored analytical-chemistry 
#' data (Helsel, 2005).  The method assumes data contains zero to many 
#' left-censored (less-than) values.  For highly censored data, \link{ros} 
#' may produce a warning message. Such as, 
#' \preformatted{Warning messages:
#' 1: In ros(my.list$obs, my.list$cen) :
#'  Input > 80\% censored -- Results are tenuous.}
#' The boxplot will still be generated, but the user should consider the
#' warning message when interpreting the plots.
#' See Oblinger Childress and others (1999) for information on the 
#' remark codes used by the U.S. Geological Survey.
#' @keywords hplot
#' @author Karen R. Ryberg
#' @export
#' @references
#' Helsel, D.R., 2005, Nondetects and data analysis: New York, John 
#' Wiley and Sons.
#' 
#' Lee, Lopaka, 2012, Nondetects and data analysis for environmental 
#' data: R package version 1.5-4, 
#' \url{http://CRAN.R-project.org/package=NADA}.
#' 
#' Oblinger Childress, C.J., Foreman, W.T., Connor, B.F., and Maloney, 
#' T.J., 1999, New reporting procedures based on long-term method 
#' detection levels and some considerations for interpretations of 
#' water-quality data provided by the U.S. Geological Survey: U.S. 
#' Geogolical Survey Open-File Report 99--193, 19 p. (Also available at 
#' \url{http://water.usgs.gov/owq/OFR_99-193/index.html}.)
#' @examples
#' data(swData)
#' # summary of water-quality concentrations
#' apply(IllRivValleyCty[, grep("P[[:digit:]]", 
#'       dimnames(IllRivValleyCty)[[2]])], 2, summary)
#' # simple boxplot of water-quality concentrations
#' rosBoxPlot(IllRivValleyCty)
#' # same boxplot function with many additional plotting arguments
#' rosBoxPlot(IllRivValleyCty, 
#'            site="05586100 Illinois River at Valley City, IL", 
#'            log="y", yaxt="n", ylim=c(0.0000001, 1), qwcols=c("R", "P"),
#'            ylab=c("Concentration, micrograms per liter"), 
#'            col="skyblue1", cex.axis=0.7, cex.sub=0.8, 
#'            par(tcl=0.5, las=1, yaxs="i", mgp=c(3, 0.5, 0), 
#'                mar=c(5, 5, 2, 2),cex.main=0.9))
#' axis(2, 
#'      at=c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
#'      labels=c("0.0000001", "0.000001","0.00001", "0.0001", "0.001",
#'               "0.01", "0.1", "1"), cex.axis=0.7)
rosBoxPlot <- function(data, site="", qwcols=c("R", "P"), ...) {
  # parameter columns
  pmatch <- paste(qwcols[2], "[[:digit:]]", sep="")
  parms <- grep(pmatch, dimnames(data)[[2]])
  # qualifier columns
  rmatch <- paste(qwcols[1], "[[:digit:]]", sep="")
  quals <- grep(rmatch, dimnames(data)[[2]])
  obs <- as.matrix(data[,parms])
  mynames <- dimnames(data)[[2]][parms]
#   if (dim(obs)[[2]]==1) {
#     colnames(obs)<-mynames
#   }
  censored <- as.matrix(data[,quals])
  censored <- censored=="<"
  my.res<-list()
  for (i in 1:length(parms) ) {
    my.list <- list(obs=obs[, i][!is.na(obs[, i])], 
                    cen=censored[, i][!is.na(obs[, i])])
    n<-length(my.list$obs)
    m <- length(my.list$cen[my.list$cen==TRUE])
    if (n > 0 & m/n < .90) {
      my.ros<-ros(my.list$obs, my.list$cen)
      my.res[[paste(mynames[i], sep="")]] <- my.ros$modeled
    }
  }
  boxplot(my.res, ...)
  title(sub="Censored values estimated using regression on order 
        statistcs", line=3)
  if (is.character(site)) {
    title(main=site)
  }
}
