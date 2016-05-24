.dist.Rmd <-
function(Y, dname, fun.call, d) {

  explain <- TRUE

  fncl <- .fun.call.deparse(fun.call) 
  fc <- .rm.arg("Rmd", fncl) 

  if (explain) show <- "" else show <- ", echo=FALSE"

  tx <- character(length = 0)

  tx[length(tx)+1] <- "---"
  tx[length(tx)+1] <- "output:"
  tx[length(tx)+1] <- "  html_document:"
  tx[length(tx)+1] <- "    fig_height: 4.5"
  tx[length(tx)+1] <- "    fig_width: 5.5"

  tx[length(tx)+1] <- "---"


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "***"

  v <- packageVersion("lessR")
  tx[length(tx)+1] <- paste(
"_", format(Sys.time(), "%a %b %d, %Y at %H:%M"), " &nbsp; with ",
"lessR version ", v, "_",
sep="")

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "***"



  tx[length(tx)+1] <- paste("# Distribution of ", Y, sep="")


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <- "suppressPackageStartupMessages(library(lessR))  # load lessR"
  tx[length(tx)+1] <- "```"


  tx[length(tx)+1] <- paste(
"The purpose of this analysis is to examine the distribution of ",
"values of the variable ",
 Y, ".",
 sep="")

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Data"


  rdcall <- getOption("read.call")
  if (is.null(rdcall)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
       "To generate a knitr output file, first read the data for this\n",
       "regression analysis with the lessR function Read.\n\n")
  }

  #ref <- .get.arg("ref", rdcall)
  #if (nchar(ref) == 0) {
      #cat("\n"); stop(call.=FALSE, "\n","------\n",
       #"To generate a knitr output file, need to specify a file name\n",
       #"to Read the data for this regression analysis.\n\n")
  #}


  tx[length(tx)+1] <- "First read the data."
  ref <- .get.arg("ref", rdcall)
  if (ref %in% c("Employee", "Reading", "Cars93", "Jackets", "Learn", "Mach4")) {
    ref  <- paste(ref, "\"", ", format=\"lessR", sep="")
    tx[length(tx)+1] <- paste(
"Here read from a data file included with the `lessR` package.",
sep="")
  }


  tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
  tx[length(tx)+1] <- paste(dname, " <- ", rdcall, sep="")
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- paste(
"Data from the following variables are available for analysis: ",
"`r xAnd(names(", dname, "))`. ",
sep="")






  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Histogram"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"Obtain the histogram and related analyses with the `lessR` function ",
"`Histogram`. Can also use the abbreviation `hs`.",
sep="")
  tx[length(tx)+1] <- "```"
  tx[length(tx)+1] <- paste("h <- Histogram(", Y, ")", sep="")  # not run here (no {r})
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"The first part of the output provides the descriptive statistics for ",
"this sample of data of ", Y, ".",
sep="")
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <- paste("s <- ss.brief(", Y, ")", sep="") 
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <-"s$out_stats" 
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- paste("\n",
"Of the `r s$n+s$n.miss` cases presented for analysis, `r s$n` are ",
"retained, so the number of deleted data values due to missing data ",
"is `r s$n.miss`. The sample mean of ", Y, " is `r xP(s$mean", d, ")` with ",
"a standard deviation of `r xP(s$sd", d, ")` ranging from ",
"`r xP(s$min", d, ")` to `r xP(s$max", d, ")`. ",
sep="")
  
  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <- paste("h <- hs(", Y, ")", sep="")
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- paste("\n",
"The histogram presents a graphical display of the frequency of the ",
"different values of ", Y, " grouped into bins. The width of each bin is ",
"`r h$bin_width`.",
sep="")

  tx[length(tx)+1] <- paste("The tabular presentation of the distribution ", 
"lists the bins, the count of the values in each bin, the corresponding ",
"proportions and the cumulative counts and proportions.",
sep="")
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <-"h$out_freq" 
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- paste("The analysis of any variable should consider if ",
"outliers are present. Here the definition of an outlier is Tukey's ",
"from the boxplot, any value more than 1.5 interquartile ranges from ",
"the median.",
sep="")
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <-"h$out_outliers" 
  tx[length(tx)+1] <- "```"





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Density Plot"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"Obtain the density plot and related analyses of ", Y,
" with the `lessR` function ",
"`Density`. Can also use the abbreviation `dn`.",
sep="")
  tx[length(tx)+1] <- "```{r}"
  tx[length(tx)+1] <- paste("Density(", Y, ", quiet=TRUE)", sep="")  # no stats
  tx[length(tx)+1] <- "```"





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Box Plot"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"Obtain the box plot and related analyses of ", Y,
" with the `lessR` function ",
"`BoxPlot`. Can also use the abbreviation `bx`.",
sep="")
  tx[length(tx)+1] <- "```{r}"
  tx[length(tx)+1] <- paste("BoxPlot(", Y, ", quiet=TRUE)", sep="")  # no stats
  tx[length(tx)+1] <- "```"





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Scatter Plot"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"Obtain the 1-dimensional scatter plot and related analyses of ", Y,
" with the ",
"`lessR` function `ScatterPlot`. Can also use the abbreviation `sp`.",
sep="")
  tx[length(tx)+1] <- "```{r}"
  tx[length(tx)+1] <- paste("ScatterPlot(", Y, ", quiet=TRUE)", sep="")  # no stats
  tx[length(tx)+1] <- "```"

  return(tx)

}
