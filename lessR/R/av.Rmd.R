.av.Rmd <-
function(nm, dname, fun.call, d) {

  explain <- TRUE
  results <- TRUE

  fncl <- .fun.call.deparse(fun.call) 
  fc <- .rm.arg("Rmd", fncl) 

  if (explain) show <- "" else show <- ", echo=FALSE"

  n.vars <- length(nm)
  n.pred <- n.vars - 1
  Y <- nm[1]
  pred <- character(length=0)
  for (i in 1:n.pred) pred[i] <- nm[i+1]
  X <- xAnd(pred)

  if (n.pred == 2) {
    if (grepl("*", fc, fixed=TRUE)) bet.grp  <- TRUE else bet.grp <- FALSE
    if (grepl("+", fc, fixed=TRUE)) wth.grp  <- TRUE else wth.grp <- FALSE
  }




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



  tx[length(tx)+1] <- paste("# ANOVA of ", Y, sep="")


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "```{r, echo=FALSE}"
  tx[length(tx)+1] <- "suppressPackageStartupMessages(library(lessR))  # load lessR"
  tx[length(tx)+1] <- "```"


  tx[length(tx)+1] <- paste(
"The purpose of this analysis is the analysis of variance of ",
"values of the variable ", Y, " for different levels of the values ",
"of ", X, ".", 
 sep="")



  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Data"

  rdcall <- getOption("read.call")
  if (is.null(rdcall)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
       "To generate a knitr output file, first read the data for this\n",
       "analysis with the lessR function Read.\n\n")
  }

  tx[length(tx)+1] <- "Read the data with the `lessR` function `Read`."
  ref <- .get.arg("ref", rdcall)
  if (ref %in% c("Employee", "Reading", "Cars93", "Jackets", "Learn", "Mach4")) {
    ref  <- paste(ref, "\"", ", format=\"lessR", sep="")
    tx[length(tx)+1] <- paste(
"Here read from a data file included with the `lessR` package.",
sep="")
  }

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
  tx[length(tx)+1] <- paste(dname, " <- ", rdcall, sep="")
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- paste(
"Data from the following variables are available for analysis: ",
"`r xAnd(names(", dname, "))`. ",
sep="")





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Analysis of Variance"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste(
"Obtain the analysis with the `lessR` function `ANOVA`.",
sep="")


  if (n.pred == 1) {
    tx[length(tx)+1] <- paste(
"This is a one-way ANOVA of ", Y, " with treatment factor ", X, ".",
 sep="")
}

  if (n.pred == 2) {
    if (bet.grp)
      tx[length(tx)+1] <- paste(
"This is a two-way between groups ANOVA of ", Y, " with ", 
"treatment factors ", X, ".",
sep="")
    else
      tx[length(tx)+1] <- paste(
"This is a randomized blocks ANOVA of ", Y, " with one treatment ",
"factor, ", pred[1], ", ",
"and one blocking factor, ", pred[2], ".",
 sep="")
}

  locTRUE <- regexec("graphics = TRUE", fc)
  if (locTRUE == -1) {
    loc <- regexec("graphics = FALSE", fc)
    if (loc == -1) fc <- sub(")$", ", graphics = FALSE)", fc)
  }

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
  tx[length(tx)+1] <- paste("a <- ", fc, sep="")
  tx[length(tx)+1] <- "```"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Background"

  if (explain) tx[length(tx)+1] <- paste(
"The output begins with a specification of the variables ",
"in the model and a brief description of the data. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"a$out_background" 
    tx[length(tx)+1] <- "```"
  }




  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Descriptive Statistics"

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste(
"For this sample of data, ",
"the first part of the output provides the descriptive statistics ",
"of ", Y, " for ",
"the different levels of ", X,
".",
sep="")

  if (n.pred == 1) {
    tx[length(tx)+1] <- "```{r, echo=FALSE}"
    tx[length(tx)+1] <-"a$out_descriptive" 
    tx[length(tx)+1] <- "```"
  }

  if (n.pred == 2) {

    if (bet.grp) {

      tx[length(tx)+1] <- paste(
"First, the number of data values in each cell. ",
sep="")

      tx[length(tx)+1] <- "```{r, echo=FALSE}"
      tx[length(tx)+1] <-"a$out_cell.n" 
      tx[length(tx)+1] <- "```"

    tx[length(tx)+1] <- paste(
"The mean of cell in the design follows. ",
sep="")

      tx[length(tx)+1] <- "```{r, echo=FALSE}"
      tx[length(tx)+1] <-"a$out_cell.means" 
      tx[length(tx)+1] <- "```"

    }

    tx[length(tx)+1] <- paste(
"The marginal means of each row and column assist in interpreting ",
"any main effects. ",
sep="")

      tx[length(tx)+1] <- "```{r, echo=FALSE}"
      tx[length(tx)+1] <-"a$out_marginals" 
      tx[length(tx)+1] <- "```"

    tx[length(tx)+1] <- paste(
"The mean of all the data, the grand mean, follows. ",
sep="")

      tx[length(tx)+1] <- "```{r, echo=FALSE}"
      tx[length(tx)+1] <-"a$out_gm" 
      tx[length(tx)+1] <- "```"

    tx[length(tx)+1] <- paste(
"Consider also the variation in each cell, assessed by the ",
"standard deviation. ",
sep="")

      tx[length(tx)+1] <- "```{r, echo=FALSE}"
      tx[length(tx)+1] <-"a$out_cell.sd" 
      tx[length(tx)+1] <- "```"

  }





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Summary Table"

  if (explain) tx[length(tx)+1] <- paste("\n",
"The analysis of variance (ANOVA) partitions this total ",
"sum of squares for ", Y, " into the residual variability, ",
"$\\sum e^2_i$, and the ",
"sum of squares for ", X, ". The ANOVA table ",
"displays these sources of variation. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"a$out_anova" 
    tx[length(tx)+1] <- "```"
  }


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Effects"

  if (explain) tx[length(tx)+1] <- paste("\n",
"The ANOVA summary table provides the significance test ",
"for ", X, ", in terms of its relationship with ", Y, ". This test ",
"shows if the effect likely exists or not. ",
"Also of interest is the extent of the effect. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"a$out_effects" 
    tx[length(tx)+1] <- "```"
  }





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Pairwise Differences"

  if (explain) tx[length(tx)+1] <- paste("\n",
"A significant effect shows if ", X, " impacts ", Y, ". ",
"But do all levels of ", X, " impact ", Y, " or just some? ",
"Answer the question by examining the pairwise comparisons, ",
"adjusted for the overall level of significance. ",
"These tests are the search for Honestly Significant Differences, ",
"HSD.",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"a$out_hsd" 
    tx[length(tx)+1] <- "```"
  }



  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Residuals"

  if (explain) tx[length(tx)+1] <- paste("\n",
"ANOVA is a kind of regression analysis. ",
"As such, individual data values can be compared to their forecasted ", 
"values, which for ANOVA are the corresponding cell means. ",
"The difference between obtained and forecasted is the residual. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"a$out_res" 
    tx[length(tx)+1] <- "```"
  }

  return(tx)

}
