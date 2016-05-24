.corfa.Rmd <-
function(mimm, nm.mimm=NULL, dname, fun.call, n.inds, n.factors, iter, item.cor, 
         explain, interpret, results) {

  fncl <- .fun.call.deparse(fun.call) 
  fc <- .rm.arg("Rmd", fncl) 

  # set parameters
  d <- 3

  if (explain) show <-  "" else show <- ", echo=FALSE"
  

  tx <- character(length = 0)

  tx[length(tx)+1] <- "---"
  #tx[length(tx)+1] <- "title: \"Factor Analysis\""
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


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "# Confirmatory Factor Analysis "

  tx[length(tx)+1] <- paste(
    "The purpose of this analysis is to identify the factors that underlie ", 
    "the observed variables, the items.",
    sep="")

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
  tx[length(tx)+1] <- "suppressPackageStartupMessages(library(lessR))  # load lessR"
  tx[length(tx)+1] <- "```"




  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Data"

  rdcall <- getOption("read.call")
  if (is.null(rdcall)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
       "To generate a knitr output file, first read the data\n",
       "with the lessR function Read.\n\n")
  }

  if (explain && !is.null(dname)) {
    tx[length(tx)+1] <- paste("Read the data with the `lessR` ",
"function `Read`. ")
    ref <- .get.arg("ref", rdcall)  # only works for Read, not rd or rd.brief
    if (ref %in% c("Employee", "Reading", "Cars93", "Jackets", "Learn", "Mach4")) {
      ref  <- paste(ref, "\"", ", format=\"lessR", sep="")
      tx[length(tx)+1] <- paste(
"Here read from a data file included with the `lessR` package.",
sep="")
    }
  }

  if (results && !is.null(dname)) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- xW(paste(dname, " <- ", rdcall, sep=""))
    tx[length(tx)+1] <- "```"
  }

  if (!is.null(dname)) tx[length(tx)+1] <- paste(
"Data from the following variables are available for analysis: ",
"`r xAnd(names(", dname, "))`. ",
sep="")



  if (explain) tx[length(tx)+1] <- paste("\n",
"The analysis proceeds from the correlations. ",
sep="")


  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- paste("mycor <- cr(", dname, ")", sep="")
    tx[length(tx)+1] <- "```"
  }





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Specified Model"

  loc <- regexec("heat.map = FALSE", fc)
  if (loc == -1) fc <- sub(")$", ", heat.map=FALSE)", fc)

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    if (!is.null(mimm))
      tx[length(tx)+1] <- paste(nm.mimm, " <- \"", format(mimm), "\"")
    tx[length(tx)+1] <- xW(paste("c <-", fc), 85)
    tx[length(tx)+1] <- "```"
  }

  if (explain) tx[length(tx)+1] <- paste(
"The output begins with a specification of the variables, the items, ",
"in the model and the content of each item if available. ",
"Each set of items forms a scale, scored as an unweighted composite. ",
"Corresponding to each observed scale score is an underlying factor. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_labels" 
    tx[length(tx)+1] <- "```"
  }







  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Scale Reliability"


  if (explain) tx[length(tx)+1] <- paste(
"Reliability is of the composite, the unweighted total score, ",
"for each scale. Alpha assumes equal item reliabilities.",
sep="")

  if (explain  &&  iter>0) tx[length(tx)+1] <- paste(
"The more generally preferred ", 
"Omega uses each indicator\'s communality in the computation of reliability. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_reliability" 
    tx[length(tx)+1] <- "```"
  }







  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Solution"

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_indicators" 
    tx[length(tx)+1] <- "```"
  }


    if (explain && item.cor) tx[length(tx)+1] <- paste(
"* _Item Correlation_: Correlation of two items with each other\n",
"* _Communality_, in the diagonal of the item correlations: Proportion",
"of the correlation of an item with itself that is due only to ",
"its underlying factor ",
sep="")

    tx[length(tx)+1] <- paste(
"* _Factor Loading_: Correlation of an item with a factor\n",
"* _Pattern Coefficient_: Regression coefficient of an item on its, ",
"  underlying factor, a special case of a factor loading\n",
"* _Factor Correlation_: Correlation of two factors with each other ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_solution" 
    tx[length(tx)+1] <- "```"
  }






  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Residuals"

    tx[length(tx)+1] <- paste(
"Difference between an item correlation and its value imposed ",
"by the estimated multiple indicator measurement model. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_residuals" 
    tx[length(tx)+1] <- "```"
  }


  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_res_stats" 
    tx[length(tx)+1] <- "```"
  }




  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## lavaan Code"

    tx[length(tx)+1] <- paste(
"lavaan code for confirmatory factor analysis of the model, ",
"fully standardized, maximum likelihood solution. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"c$out_lavaan" 
    tx[length(tx)+1] <- "```"
  }

  return(tx)

}
