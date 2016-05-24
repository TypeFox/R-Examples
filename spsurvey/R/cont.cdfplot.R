cont.cdfplot <- function(pdffile="cdf2x2.pdf", cdfest, units.cdf="Percent",
   ind.type=rep("Continuous", nind), logx=rep("", nind), xlbl=NULL,
   ylbl="Percent", ylbl.r=NULL, legloc="BR", cdf.page=4, width=10, height=8,
   confcut=5, cex.main=1.2, ...) {

################################################################################
# Function: cont.cdfplot
# Purpose: Create CDF plots from the cont.analysis output object
# Programmers: Tony Olsen
#              Tom Kincaid
# Date: March 26, 2007
# Last Revised: June 23, 2014
# Description:
#   This function creates CDF plots.  Input data for the plots is provided by a
#   data frame utilizing the same structure as the data frame named "CDF" that
#   is included in the output object produced by function cont.analysis.  Plots
#   are produced for every combination of Type of population, subpopulation
#   within Type, and indicator.  Output from the function is placed in a PDF
#   file.
# Arguments:
#   pdffile = name of the PDF file.  The default is "cdf2x2.pdf".
#   cdfest = data frame utilizing the same structure as the data frame named
#     "CDF" that is included in the output object produced by function
#     cont.analysis.
#   units.cdf = indicator for the type of units in which the CDF is plotted,
#     where "Percent" means the plot is in terms of percent of the population,
#     and "Units" means the plot is in terms of units of the population.  The
#     default is "Percent".
#   ind.type = character vector consisting of the values "Continuous" or
#     "Ordinal" that controls the type of CDF plot for each indicator.  The
#     default is "Continuous" for every indicator.
#   logx = character vector consisting of the values "" or "x" that controls
#     whether the x axis uses the original scale ("") or the base 10 logarithmic
#     scale ("x") for each indicator.  The default is "" for every indicator.
#   xlbl = character vector consisting of the x-axis label for each indicator.
#     If this argument equals NULL, then indicator names are used as the labels.
#     The default is NULL.
#   ylbl =  character string providing the y-axis label.  The default is
#     "Percent".
#   ylbl.r = character string providing the label for the right side y-axis,
#     where NULL means a label is not created, and "Same" means the label is the
#     same as the left side label (i.e., argument ylbl).  The default is NULL.
#   legloc = indicator for location of the plot legend, where "BR" means bottom
#     right, "BL" means bottom left, "TR" means top right, and "TL" means top
#     left.  The default is "BR". 
#   cdf.page = number of CDF plots on each page, which must be chosen from the
#     values: 1, 2, 4, or 6.  The default is 4.
#   width = width of the graphic region in inches.  The default is 10.
#   height = height of the graphic region in inches.  The default is 8.
#   confcut = numeric value that controls plotting confidence limits at the CDF
#     extremes.  Confidence limits for CDF values (percent scale) less than
#     confcut or greater than 100 minus confcut are not plotted.  A value of
#     zero means confidence limits are plotted for the complete range of the
#     CDF.  The default is 5.
#   cex.main = expansion factor for the plot title.  The default is 1.2.
#   ... = additional arguments passed to the cdf.plot function.
# Output:
#   A PDF file containing the CDF plots.
# Other Functions Required:
#   cdf.plot - plot the CDF and associated confidence limits
# Example:
#   mysiteID <- paste("Site", 1:100, sep="")
#   mysites <- data.frame(siteID=mysiteID, Active=rep(TRUE, 100))
#   mysubpop <- data.frame(siteID=mysiteID, All.Sites=rep("All Sites",100),
#      Resource.Class=rep(c("Good","Poor"), c(55,45)))
#   mydesign <- data.frame(siteID=mysiteID, wgt=runif(100, 10, 100),
#      xcoord=runif(100), ycoord=runif(100), stratum=rep(c("Stratum1",
#      "Stratum2"), 50))
#   ContVar <- rnorm(100, 10, 1)
#   mydata.cont <- data.frame(siteID=mysiteID, ContVar=ContVar)
#   mypopsize <- list(All.Sites=c(Stratum1=3500, Stratum2=2000),
#      Resource.Class=list(Good=c(Stratum1=2500, Stratum2=1500),
#      Poor=c(Stratum1=1000, Stratum2=500)))
#   myanalysis <- cont.analysis(sites=mysites, subpop=mysubpop, design=mydesign,
#      data.cont=mydata.cont, popsize=mypopsize)
#   cont.cdfplot("myanalysis.pdf", myanalysis$CDF, ylbl.r="Stream Length (km)")
################################################################################

# Open the PDF file

pdf(file=pdffile, width=width, height=height)

# Set up the number of plots per page

if(cdf.page == 1) mf <- c(1,1)
else if(cdf.page == 2) mf <- c(2,1)
else if(cdf.page == 4) mf <- c(2,2)
else if(cdf.page == 6) mf <- c(3,2)
else mf <- c(2,2)

# Set graphical parameter values

op <- par(mfrow=mf, mgp=c(1.5,0.7,0))

# Assign the vectors of population Type names and indicator names

typenames <- unique(cdfest$Type)
indnames <- unique(cdfest$Indicator)
nind <- length(indnames)

# If not supplied, set up the x-axis labels

if(is.null(xlbl)) {
   xlbl <- as.character(indnames)
   names(xlbl) <- as.character(indnames)
}

# Obtain the confidence level

conflev <- as.numeric(substr(names(cdfest)[8], 4, 5))

# Create the plots

for(itype in 1:length(typenames)) {
   tsttype <- cdfest$Type == typenames[itype]
   subnames <- unique(cdfest$Subpopulation[tsttype])
   for(jsub in 1:length(subnames)) {
      tstsub <- tsttype & cdfest$Subpopulation == subnames[jsub]
      temp <- match(unique(cdfest$Indicator[tstsub]), indnames)
      for(kin in temp) {
	    tstind <- tstsub & cdfest$Indicator == indnames[kin]
         cdf.plot(cdfest[tstind,], units.cdf=units.cdf, type.cdf=ind.type[kin],
            logx=logx[kin], xlbl=xlbl[indnames[kin]], ylbl=ylbl, ylbl.r=ylbl.r, 
            figlab=paste(typenames[itype], " - ", subnames[jsub], ": ",
            indnames[kin], sep=""), legloc=legloc, confcut=confcut,
            conflev=conflev, cex.main=cex.main, ...)
      }
   }
}


# Reset graphical parameter values

par(op)

# Close the PDF file

graphics.off()

invisible(NULL)
}
