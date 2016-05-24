.reg.Rmd <-
function(nm, dname, fun.call, res.rows, pred.rows, res.sort,
         digits.d, results, explain, interpret, document, code,
         pvalues, tolerances, resid.max, numeric.all, X1.new,
         new.val=matrix(nrow=n.vars-1, ncol=2, byrow=TRUE)) {

  fncl <- .fun.call.deparse(fun.call) 
  if (regexec("Rmd", fncl)[1] > 0) fc <- .rm.arg("Rmd", fncl) 
  if (regexec("explain", fc)[1] > 0) fc <- .rm.arg.ns("explain", fc) 
  if (regexec("interpret", fc)[1] > 0) fc <- .rm.arg.ns("interpret", fc) 
  if (regexec("results", fc)[1] > 0) fc <- .rm.arg.ns("results", fc) 

  show <- ifelse (code, "", ", echo=FALSE")

  # set parameters
  n.vars <- length(nm)
  n.pred <- n.vars - 1
  d <- digits.d
  Y <- nm[1]
  pred <- character(length=0)
  for (i in 1:n.pred) pred[i] <- nm[i+1]
  X <- xAnd(pred)

  # get variable labels and units if exist
  mylabels <- attr(get(dname, pos=.GlobalEnv), which="variable.labels")
  myunits <- attr(get(dname, pos=.GlobalEnv), which="variable.units")
  var.lbl <- character(length=0)
  var.unit <- character(length=0)
  for (i in 1:n.vars) {
    var.lbl[i] <- mylabels[which(names(mylabels) == nm[i])]
    var.unit[i] <- myunits[which(names(myunits) == nm[i])] 
  }

  if (n.pred > 1) {
    pl <- "s" 
    et <- "each "
    cnst <- ", with the values of all remaining predictor variables held constant"
  }
  else {
    pl <- ""
    et <- "the "
    cnst <- ""
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

  tx[length(tx)+1] <- paste("\n",
"_Output Options: explain=", explain, ", interpret=", interpret,
", results=", results, ", document=", document, ", code=", code, "_",
sep="")

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "***"



  tx[length(tx)+1] <- ""
  if (n.pred > 1)
    tx[length(tx)+1] <- paste("# Multiple Regression of ", Y, sep="")
  else
    tx[length(tx)+1] <- paste("# Regression of ", Y, " on ", X, sep="")

  tx[length(tx)+1] <- "```{r echo=FALSE}"
  tx[length(tx)+1] <- "suppressPackageStartupMessages(library(lessR))  # load lessR"
  tx[length(tx)+1] <- "```"

  if (!is.na(var.lbl[1])) {
    Ylbl <- var.lbl[1]
    txtY <- paste(", ", Ylbl, sep="")
  }
  else
    txtY <- ""

  txtX <- ""
  if (n.pred == 1) {
    if (!is.na(var.lbl[2])) {
      Xlbl <- var.lbl[2]
      txtX <- paste(", ", Xlbl, sep="")
    }
  }

  tx[length(tx)+1] <- paste(
"The variable of primary interest is the _response variable_, ",
Y, txtY, ". ",
"The purpose of this analysis is to account for the values of ",
Y, " in ",
" terms ", "of the values of the _predictor variable", pl, "_ ",
X, txtX, ".", 
sep="")







  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Data"

  rdcall <- getOption("read.call")
  if (is.null(rdcall)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
       "To generate an R markdown file, first read the data for this\n",
       "regression analysis with the lessR function Read.\n\n")
  }


  if (document) {
    tx[length(tx)+1] <- "Read the data with the `lessR` function `Read`."
    ref <- .get.arg("ref", rdcall)  # only works for Read, not rd or rd.brief
    if (ref %in% c("Employee", "Reading", "Cars93", "Jackets", "Learn", "Mach4")) {
      ref  <- paste(ref, "\"", ", format=\"lessR", sep="")
      tx[length(tx)+1] <- paste(
"Here read from a data file included with the `lessR` package.",
sep="")
    }
  }
  else
    tx[length(tx)+1] <- "Read the data."
  

  if (explain) tx[length(tx)+1] <- paste(
"The corresponding data values for the ",
"variables in the model comprise the _training ",
"data_, from which the model is estimated. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- paste(dname, " <- ", rdcall, sep="")
    tx[length(tx)+1] <- "```"
  }

  tx[length(tx)+1] <- paste(
"Data from the following variables are available for analysis: ",
"`r xAnd(names(", dname, "))`. ",
sep="")







  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## The Model"

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Specified Model"

  tx[length(tx)+1] <- paste(
"Express ", Y, " as a linear function of ", xNum(n.pred), " ", 
"predictor variable", pl, ": ", X, ". ",
sep="")

  if (explain) tx[length(tx)+1] <- paste(
"Within the context of the model, indicate the response variable with a Y, ",
"subscripted by the variable name, $Y_{", Y, "}$. Write each predictor variable as ",
"a subscript to an X. ",
"From the training data compute $\\hat Y_{", Y, "}$, the _fitted value_ of the ",
"response variable from the model for ",
"a specific set of values for ", X, ". ",
sep="")

  cv <- paste("$$\\hat Y_{", Y, "} = b_0 + b_1 X_{", nm[2], "}", sep="")
  if (n.vars > 2)
    for (i in 3:n.vars) cv <- paste(cv, " + b_", i-1, " X_{", nm[i], "}", sep="")
  cv <- paste(cv, "$$", sep="")
  tx[length(tx)+1] <- cv

    if (explain) tx[length(tx)+1] <- paste(
"The _intercept_, $b_0$, indicates ",
"the fitted value of ", Y, ", for values of ", X, " all equal to zero.", sep="")

    if (n.pred > 1) {
        txt2 <- paste("through $b_", n.pred, "$", sep="")
      }
    else {
      txt2 <- ""
    }

      if (explain) tx[length(tx)+1] <- paste(
xU(et), "_slope coefficient_, ",
"$b_1$ ", txt2, ", is the ",
"average change in the value of ",
"response variable, ", Y, ", for a one-unit increase in the value of ",
"the corresponding predictor variable",
cnst, ". The values of these estimated coefficients only apply to ",
"the interpretation of the training data from which they were estimated. ",
sep="")

    if (explain) tx[length(tx)+1] <- paste("\n",
"To compute $\\hat Y_{", Y, "}$ from a specific set of values ",
"for ", X, " requires the estimated values of ", 
"the coefficients of the model, the values of ",
"each regression coefficient, $b_j$. ",
"This estimation procedure ",
"depends on the _residual_, the difference between the ",
"actual value of ", Y, " for each row of data ",
"and the corresponding value fitted by the model. ", 
"Define the residual as a variable across all rows of data. ",
"Use the subscript _i_ ",
"to indicate the $i^{th}$ row of data ",
"to emphasize that the expression applies to _each_ row of training data. ",
"The name of the response variable in this notation is understood ",
"and so is omitted for simplicity. ",
sep="")

    if (explain) tx[length(tx)+1] <- paste(
"$$e_i = Y_i - \\hat Y_i$$",
sep="")

    if (explain) tx[length(tx)+1] <- paste(
"Estimate the coefficients with ordinary least squares (OLS), ",
"which provides the one set of estimates that minimize the ",
"sum of the squared ",
"residuals, $\\sum e^2_i$, across all the rows of the training data. ",
sep="")

    if (document) tx[length(tx)+1] <- paste("\n",
"Accomplish the estimation and related computations with the `lessR` ",
"function `Regression`, abbreviated as `reg`. ",
"Keep graphics separate, so generate these later. ",
sep="")

  locTRUE <- regexec("graphics = TRUE", fc)
  if (locTRUE == -1) {
    loc <- regexec("graphics = FALSE", fc)
    if (loc == -1) fc <- sub(")$", ", graphics = FALSE)", fc)
  }

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- paste("r <-", fc)
    tx[length(tx)+1] <- "```"
  }

  if (explain) tx[length(tx)+1] <- paste(
"The output begins with a specification of the variables ",
"in the model and a brief description of the data. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"r$out_background" 
    tx[length(tx)+1] <- "```"
  }

  if (interpret) {
    tx[length(tx)+1] <- paste("\n",
"Of the `r r$n.obs` cases presented for analysis, `r r$n.keep` are retained, ",
"so the number of deleted cases due to missing data is `r r$n.obs - r$n.keep`.",
sep="")
  }




  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Estimated Model"

  if (explain) tx[length(tx)+1] <- paste(
"The analysis of the model begins with the estimation of each sample ",
"regression coefficient, $b_j$, from the training data. ",
"Of greater interest is each corresponding population value, $\\beta_j$, ",
"in the _population model_. ",
sep="")

  cv <- paste("$$\\hat Y_{", Y, "} = \\beta_0 + \\beta_1 X_{", nm[2], "}", sep="")
  if (n.vars > 2)
    for (i in 3:n.vars) cv <- paste(cv, " + \\beta_", i-1, " X_{", nm[i], "}", sep="")
  cv <- paste(cv, "$$", sep="")
  if (explain)tx[length(tx)+1] <- cv

  if (explain) tx[length(tx)+1] <- paste(
"The associated inferential analyses for each estimate are the ",
"hypothesis test and confidence interval.",
sep="")

  if (explain) tx[length(tx)+1] <- paste(
"Each _t_-test evaluates the _null hypothesis_ that ",
"the corresponding _individual_ population regression ",
"coefficient is 0, here for the $j^{th}$ coefficient. ",
sep="")

  if (explain) tx[length(tx)+1] <- paste(
"$$H_0: \\beta_j=0$$\n",
"$$H_1: \\beta_j \\ne 0$$",
sep="")

  if (explain) tx[length(tx)+1] <- paste(
"The confidence interval provides ",
"the range of likely values of the corresponding $\\beta_j$. ",
"Each 95% confidence interval is the margin of error on ",
"either side of the corresponding ",
"estimated intercept or slope coefficient, $b_j$. ",
sep="")
  
  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- "r$out_estimates"
    tx[length(tx)+1] <- "```"
  }

  tx[length(tx)+1] <- paste(
"This estimated model is the specific linear function that yields a ",
"fitted value of ", Y, " from the provided value", pl, " of ", X, ".",
sep="")

  cv <- paste(
"$$\\hat Y_{", Y, "} = `r xP(r$coefficients[1],", d, ")` ", 
"`r ifelse(sign(r$coefficients)==1, \"+\", \"-\")[2]` ",
"`r xP(abs(r$coefficients[2]),", d, ")` X_{", nm[2], "}", 
sep="")
  if (n.vars > 2)
    for (i in 3:n.vars)
      cv <- paste(
cv, "`r ifelse(sign(r$coefficients)==1, \"+\", \"-\")[", i, "]` ",
" `r xP(abs(r$coefficients[", i, "]),", d, ")`",
" X_{", nm[i], "}",
sep="")
  cv <- paste(cv, "$$", sep="")
  tx[length(tx)+1] <- cv

  gt05 <- length(which(pvalues[2:length(pvalues)] > 0.05)) 
  if (gt05 > 0) {
    if (gt05  > 1) {
      txt1 <- "these"
      txt3 <- "have _p_-values"
      txt5 <- "each "
      pl2 <- "s"
    }
    else {
      txt1 <- "this"
      txt3 <- "has a _p_-value"
      txt5 <- "the "
      pl2 <- ""
    }

    if (interpret) tx[length(tx)+1] <- paste(
xU(xNum(gt05)), " predictor variable", pl2, " ", txt3, " larger than ",
"$\\alpha$ = 0.05: ", 
"`r xAnd(names(which(r$pvalues[2:length(r$pvalues)] > 0.05)))`. ", 
sep="")

    if (interpret) tx[length(tx)+1] <- paste(
xU(txt5), "null hypothesis of no ",
"relationship could not be rejected, so there is a reasonable possibility ",
"that ", txt5, " predictor variable may not contribute to ",
"explaining the values of ", Y, cnst, ". ",
sep="")
  }  #  p > .05


  n.sig <- length(which(pvalues[2:length(pvalues)] <= 0.05)) 
  if (n.sig > 0) {

    if (length(which(pvalues[2:length(pvalues)] <= 0.05)) > 1) {
      txt1 <- "These predictor variables each have "
      txt2 <- "their"
      txt3 <- "these corresponding coefficients"
      txt4 <- "these"
      pl3 <- "s"
    }
    else {
      txt1 <- "This predictor variable has "
      txt2 <- "its"
      txt3 <- "this corresponding coefficient"
      txt4 <- "this"
      pl3 <- ""
    }

    if (interpret) tx[length(tx)+1] <- paste("\n",
txt1, "a _p_-value less than or equal to $\\alpha$ = 0.05: ", 
"`r xAnd(names(which(r$pvalues[2:length(r$pvalues)] <= 0.05)))`. ",
sep="")

    if (interpret && n.pred > 1  && n.sig < n.pred)
      tx[length(tx)+1] <- paste(
"The possibility should be further explored in the remainder of this ",
"analysis that ", txt4," ", xNum(n.sig), " variable", pl3, " ", 
"may form an equally effective ",
"but more parsimonious model in terms of ", txt2, " ",
"cumulative contribution to explaining the values of ", Y, ", ", 
"compared to the current model with ", xNum(n.pred), " predictor variables. ", 
sep="")

    if (interpret) tx[length(tx)+1] <- paste("\n",
"To extend the ",
"results beyond this sample to the population from which the sample ",
"was obtained, interpret the meaning of ", txt3, " ",
"in terms of ", txt2, " confidence interval", pl3, ". ",
sep="")

    if (interpret) {

        # response variable labels/units
        l.nm <- which(names(mylabels) == Y)
        if (!is.null(mylabels[l.nm])) {
          if (nzchar(mylabels[l.nm])) {
            lY <- mylabels[which(names(mylabels) == Y)]
            txtY <- lY
          }
        }
        else
          txtY <- Y

        # if process units, then labels must already exist
        if (!is.null(myunits[which(names(myunits) == Y)])) {
          uY <- myunits[which(names(myunits) == Y)]
          uYq <- paste("\"", uY, "\"", sep="")  # unit with quotes
          if (uY != "dollar")
            txtY <- paste("the value of", uY, "of", lY)
          else
            txtY <- lY
        }
        else
          uYq <- ""

      # sig predictor variable labels
      for (i in 1:n.sig) { 
        j <- which(pvalues[2:length(pvalues)] <= 0.05)[i] 
        if (i == 1 && n.pred > 1) tx[length(tx)+1] <- ""
        if (n.pred > 1) tx[length(tx)+1] <- paste( "* _", pred[j], "_: ", sep="")

        l.nm <- which(names(mylabels) == pred[j])
        if (!is.null(mylabels[l.nm])) {
          if (nzchar(mylabels[l.nm])) {
            l <- mylabels[which(names(mylabels) == pred[j])]
            txt <- paste("unit of the value of", l)
          }
        }
        else
          txt <- paste("unit of the value of", pred[j])

        # sig predictor variable units
        u.nm <- which(names(myunits) == pred[j])
        if (!is.null(myunits[u.nm])) {
          if (nzchar(myunits[u.nm])) {
            u <- myunits[u.nm]
            txt <- paste(u, "of", l)
          }
          else
            txt <- paste("unit of the value of", pred[j])
        }

        tx[length(tx)+1] <- paste(
"With 95% confidence, for each additional ",
txt,
", on average, ", txtY, " changes somewhere between ",
"`r xP(r$cilb[", j+1, "],", d, ",", uYq, ")`", " to ",
"`r xP(r$ciub[", j+1, "],", d, "," ,uYq, ")`",
sep="") 
        remn <- sub(pred[j], "", pred)  # leaves an empty vector value
        remain <- character(length=0)
        for (k in 1:n.pred)  # collate into a single string
          if (nzchar(remn[k])) remain[length(remain)+1] <- remn[k]
        remain <- xAnd(remain)
        if (n.pred > 1)
          tx[length(tx)] <- paste(tx[length(tx)], ", with the values of ",
            remain, " held constant",
  sep="")
        tx[length(tx)] <- paste(tx[length(tx)], ".", sep="")
      }  #  i in n.sig
    }  # end interpret

  }  #  n.sig > 0, i.e., at least one predictor significant

  else
    uYq <- ""
  




  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Model Fit"

  if (explain) tx[length(tx)+1] <- paste(
"An estimated model is not necessarily a useful model. ",
sep="")

  if (n.pred == 1)
    geom <- "line"
  else if (n.pred == 2)
    geom <- "plane"
  else
    geom <- "surface"

  if (explain) tx[length(tx)+1] <- paste(
"A preliminary consideration is the fit of the model, ",
"based on the extent that the ",
"values of ", Y, " fitted by the model to the training data match the ",
"corresponding training data values of ", Y, ". Are the residuals typically ",
"close to their mean of zero, or are they scattered about ",
"the regression ", geom, " with ",
"relatively large positive and negative values? ",
sep="")



  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Partitioning of Variance"


  if (explain) tx[length(tx)+1] <- paste("\n",
"The analysis of fit depends on the adequacy of the model to account for ",
"the variability of the data values of ", Y, ", expressed in model notation ",
"as  $Y_{", Y, "}$. ",
"The core component of variability is the _sum of squares_, short for ",
"the sum of some type of squared deviations. ",
"The _total variability_ of $Y_{", Y, "}$ depends on ",
"the deviations of its data values from its mean, ",
"$Y_{", Y, "} - \\bar Y_{", Y, "}$, and then the resulting sums of squares, $SS_{", Y, "}$. ",
sep="")

  if (explain) tx[length(tx)+1] <- paste("\n",
"The analysis of the residuals, ",
"$e = Y_{", Y, "} - \\hat Y_{", Y, "}$, follows from the ",
"corresponding sum of squares, ", 
"the value minimized by the least squares estimation procedure, ",
"$\\sum e^2_i$ = $SS_{Residual}$. ",
"This residual sum of squares represents variation of $Y_{", Y, "}$ _not_ ",
"accounted for by $\\hat Y_{", Y, "}$. ",
"The complement to the residual sum of squares is the ",
"Model (or Regression) sum of squares, ",
"the deviations of the fitted values about the mean, ",
"$\\hat Y_{", Y, "} - \\bar Y_{", Y, "}$. ", 
sep="")

  if (explain) tx[length(tx)+1] <- paste("\n",
"The analysis of variance (ANOVA) partitions this total ",
"sum of squares into the residual variability, $\\sum e^2_i$, and the Model ",
"sum of squares, $SS_{Model}$. ",
"The ANOVA table ",
"displays these various sources of variation. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"r$out_anova" 
    tx[length(tx)+1] <- "```"
  }

  if (results) tx[length(tx)+1] <- paste(
"$$SS_{",Y,"} = SS_{Model} + SS_{Residual} = ",
"`r xP(r$anova_model[\"ss\"],", d, ")` + ",
"`r xP(r$anova_residual[\"ss\"],", d, ")` = ",
"`r xP(r$anova_total[\"ss\"],", d, ",", uYq, ", semi=TRUE)` $$",
sep="")

  if (explain && n.pred == 1) tx[length(tx)+1] <- paste("\n",
"This decomposition of the sums of squares ",
"of ", Y, " into what is explained by the model and what is ",
"not explained is fundamental to assessing the fit of the model. ",
sep="")

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste("\n",
"The ANOVA further partitions ",
"the overall model sums of squares by predictor variable. ",
sep="")

  cv <- paste("$$SS_{Model} = SS_{", nm[2], "}", sep="")
  if (n.vars > 2)
    for (i in 3:n.vars) cv <- paste(cv, " + SS_{", nm[i], "}", sep="")

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste(
cv, " = `r xP(r$anova_model[\"ss\"],", d, ")`$$",
sep="")

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste("\n",
"The sum of squares for a predictor variable ",
"is called a _sequential sums of squares_. ",
"It represents the effect of a predictor variable ", 
"after the effects of all previously entered variables ",
"in the model have already been accounted for. ",
"Unless the predictor variables are uncorrelated, its value depends ",
"on the sequence of the variables as specified in the model. ",
"Usually the interpretation of a sequential effect is more useful ",
"when the variables ",
"are entered in order of perceived importance. ",
sep="")

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste("\n",
"Progressing through the table of the sequential sums of squares for each ",
"predictor variable from the first entry, ", pred[1], ", through ",
"the last last entry, ", pred[n.pred],  ", forms ",
"a sequence of increasingly larger _nested models_ that ",
"successively contain more variables. ",
"For example, the _p_-value of ",
"`r r$pvalues[r$n.vars]` ",
"for the last variable entered into the model, ",
"`r all.vars(r$formula)[r$n.vars]`, is the same for both the ANOVA ",
"table and its regression slope coefficient because in both ",
"situations the effects of all other predictor variables are ",
"partialled out. ",
sep="")

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste("\n",
"This fundamental relationship of these various sums of squares ",
"components provides the basis for assessing fit.",
sep="")





  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Fit Indices"

  if (explain) tx[length(tx)+1] <- paste("\n",
"From the ANOVA two types of primary indicators of fit are derived: ",
"standard deviation of the residuals and several ",
"$R^2$ statistics. ",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"r$out_fit" 
    tx[length(tx)+1] <- "```"
  }

  if (explain) tx[length(tx)+1] <- paste("\n",
"The _standard deviation of the residuals_, $s_e$, directly assesses ",
"the variability of the data values of ", Y, " about the ",
"corresponding fitted values for the training data, the particular ",
"data set from which the model was estimated. ",
"Each mean square in the ANOVA table is a variance, a ",
"sum of squares divided by the corresponding degrees of freedom, _df_. ",
"By definition, the standard deviation, ",
"$s_e$ is the square root of the mean square of the residuals. ",
sep="")

  if (explain) tx[length(tx)+1] <- "$$s_e = "
  if (explain) tx[length(tx)] <- paste(tx[length(tx)], 
"\\sqrt{MS_{Residual}} = ",
"\\sqrt{`r xP(r$anova_residual[\"ms\"],", d, ")`} = ",
sep="")
  if (explain) tx[length(tx)] <- paste(tx[length(tx)],
"`r xP(r$se,", d, ",", uYq, ", semi=TRUE)`$$",
sep="")

  if (interpret) tx[length(tx)+1] <- paste(
"To interpret $s_e$ = `r xP(r$se,", d, ",", uYq, ")`, consider the estimated range ",
"of 95% of the values of a normally distributed variable, which ",
"depends on the corresponding 2.5% cutoff from the $t$-distribution ",
"for df=`r r$anova_residual[\"df\"]`: ",
"`r xP(-qt(0.025, df=r$anova_residual[\"df\"]),3)`. ",
sep="")

    if (interpret) tx[length(tx)+1] <- paste(
"$$95\\% \\;  Range: 2 * t_{cutoff} * s_e = ",
"2 * `r xP(-qt(0.025, df=r$anova_residual[\"df\"]),3)` * `r xP(r$se,", d, ")` = 
`r xP(r$resid_range,", d, ",", uYq, ", semi=TRUE)`$$",
sep="")

    if (interpret) tx[length(tx)+1] <- paste("\n",
"This range of the residuals for the fitted values is the lower limit ",
"of the range of forecasting error presented later. ",
sep="")

  if (explain) tx[length(tx)+1] <- paste("\n",
"A second type of fit index is ",
"$R^2$, the proportion of the overall variability of ",
"response variable ", Y, " that is accounted for by the model, ",
"applied to the training data, expressed either in terms of ",
"$SS_{Residual}$ or $SS_{Model}$. ",
sep="")

  if (explain) tx[length(tx)+1] <- "$$R^2 = "
  if (explain) tx[length(tx)] <- paste(tx[length(tx)], 
"1 - \\frac{SS_{Residual}}{SS_{", Y, "}} = ",
"\\frac{SS_{Model}}{SS_{", Y, "}} = ",
"\\frac{`r xP(r$anova_model[\"ss\"],", d, ")`} ",
"{`r xP(r$anova_total[\"ss\"],", d, ")`} = ",
sep="")
  if (explain) tx[length(tx)] <- paste(tx[length(tx)],
"`r xP(r$Rsq,3)` $$ ",
sep="")

  if (explain) tx[length(tx)+1] <- paste(
"Unfortunately when any new predictor variable is added ",
"to a model, useful or not, $R^2$ necessarily increases. ",
"Use the adjusted version, $R^2_{adj}$, to more appropriately ",
"compare ",
"models estimated from the same training data with different ",
"numbers of predictors. ",
"$R^2_{adj}$ helps to avoid overfitting a model because it only ",
"increases if a new predictor variable added to the model ",
"improves the fit more than would be expected by ",
"chance. The adjustment considers the ",
"number of predictor variables relative to the number of rows of data ",
"(cases). Accomplish this adjustment with the degrees of freedom, ",
"to shift from the Sum of Squares ",
"to the corresponding Mean Squares. ",
sep="")

  if (explain) tx[length(tx)+1] <- "$$R^2_{adj} = "
  if (explain) tx[length(tx)] <- paste(tx[length(tx)], 
"1 - \\frac{SS_{Residual} \\; / \\; `r r$anova_residual[\"df\"]`}{SS_{", Y, "} \\; / \\; `r r$anova_total[\"df\"]`} = ",
"1 - \\frac{MS_{Residual}}{MS_{", Y, "}} = ",
"1 - \\frac{`r xP(r$anova_residual[\"ms\"],", d, ")`} ",
"{`r xP(r$anova_total[\"ms\"],", d, ")`} = ", 
sep="")
  if (explain) tx[length(tx)] <- paste(tx[length(tx)],
"`r xP(r$Rsqadj,3)`$$",
sep="")

  if (interpret) tx[length(tx)+1] <- paste(
"From this analysis compare $R^2$ = `r xP(r$Rsq,3)` to the ",
"adjusted value of $R^2_{adj}$ = `r xP(r$Rsqadj,3)`, a difference of ", 
"`r xP((r$Rsq-r$Rsqadj), 3)`. A large difference indicates that too many ",
"predictor variables in the model for the available data yielded an overfitted ",
"model. ",
sep="")

  if (explain)  tx[length(tx)+1] <- paste("\n",
"Both $R^2$ and $R^2_{adj}$ describe the fit of the model to the training ",
"data. Base the fit of the model to forecasts from new data ",
"on the _predictive residual_ (PRE). ",
"To calculate this residual for a row of data (case), first ",
"estimate a model with the case deleted, that is, from only all the ",
"remaining cases in the training data, ",
"what is called a _case-deletion_ statistic. Repeat for all rows of data. ",
"$SS_{PRE}$, or PRESS, ",
"is the sum of squares of all the predictive residuals in a data set. ",
"From $SS_{PRE}$ define the predictive $R^2$. ",
sep="")

  if (explain) tx[length(tx)+1] <- "$$R^2_{PRESS} = "
  if (explain) tx[length(tx)] <- paste(tx[length(tx)], 
"1 - \\frac{SS_{PRE}}{SS_{", Y, "}} = ",
"1 - \\frac{`r xP(r$PRESS,", d, ")`} ",
"{`r xP(r$anova_total[\"ss\"],", d, ")`} = ",
sep="")
  if (explain) tx[length(tx)] <- paste(tx[length(tx)],
"`r xP(r$RsqPRESS,3)` $$ ",
sep="")

  if (interpret) tx[length(tx)+1] <- paste("\n",
"Because an estimated model at least to some extent overfits the training data, ",
"the more useful ",
"$R^2_{PRESS}$ = `r xP(r$RsqPRESS,3)` is lower than both $R^2$ and $R^2_{adj}$. ",
sep="")




  tx[length(tx)+1] <- ""
  if (explain && n.pred > 1) 
    tx[length(tx)+1] <- "### Hypothesis Tests of Multiple Coefficients"

  if (explain && n.pred > 1) tx[length(tx)+1] <- paste(
"The ANOVA table presents the overall ",
"hypothesis test that evaluates if _all_ the predictor variables ",
"as a set -- ", X, " -- ", 
"are related to ", Y, " as specified by the model. ",
sep="")

  cv <- paste("\\beta_{", nm[2], "} = \\beta_{", nm[3], "}", sep="")
  if (n.vars > 2)
    for (i in 3:n.vars)
      cv <- paste(cv, " = \\beta_{", nm[i], "}", sep="")
  cv <- paste(cv, "= 0", sep="")

  if (explain  &&  n.pred > 1)  tx[length(tx)+1] <- paste(
"$$\n",
"\\begin{aligned}\n",
"H_0&: \\;", cv, " \\\\\\\\ \n",
"H_1&: \\; at \\; least \\;  one \\;  \\beta_j \\ne 0\n",
"\\end{aligned}\n",
"$$",
sep="")

  if (explain  &&  n.pred > 1) tx[length(tx)+1] <- paste(
"From the sums of squares for the Model and Residuals, ",
"with degrees of freedom of ",
"`r as.integer(r$anova_model[\"df\"])` ",
"and `r as.integer(r$anova_residual[\"df\"])`, ",
"the test statistic is _F_ = `r xP(r$anova_model[\"fvalue\"],", d, ")`, ", 
"with a _p_-value of `r xP(r$anova_model[\"pvalue\"],", d, ")`." ,  
sep="")

  if (explain && n.pred > 2) tx[length(tx)+1] <- paste("\n",
"To help identify predictors that contribute little beyond ",
"those variables previously included in the model, generally ",
"list the more important variables first in the model specification. ",
"Add together the sequential sum of squares from the ANOVA table ", 
"for variables listed last in the table to form a nested model. ",
"Then test if the designated ",
"_subset_ of the regression coefficients are all equal to zero. ",
"To illustrate, consider the hypothesis test that ",
"the slope coefficients for the last two variables, ",
nm[n.vars-1], " and ", nm[n.vars], ", are both equal to zero. ",
sep="")

  cv <- paste("\\beta_{", nm[n.vars-1], "} = \\beta_{", nm[n.vars], "}", sep="")
  cv <- paste(cv, "= 0", sep="")

  if (explain && n.pred > 2) tx[length(tx)+1] <- paste(
"$$\n",
"\\begin{aligned}\n",
"H_0&: \\;", cv, " \\\\\\\\ \n",
"H_1&: \\; at \\; least \\;  one \\;  \\beta_{", nm[n.vars-1], "}, \\beta_{", nm[n.vars], "} \\ne 0\n",
"\\end{aligned}\n",
"$$",
sep="")

  if (document) tx[length(tx)+1] <- paste("\n",
"Compare two nested models with the `lessR` function `Nest`. Specify the ",
"response variable ", Y, ", the variables in ",
"the reduced model, and then the additional variables in the full model. ",
"`Nest` also ensures that the same data values are compared when there ",
"is missing data that might otherwise leave more data in the analysis of ", 
"the reduced model. ",
sep="")

  cv <- "n <- Nest("
  cv <- paste(cv, nm[1], ", c(", sep="")
  for (i in 2:(n.vars-2)) {
    txt <- ","
    if (i == n.vars-2) txt <- ")"
    cv <- paste(cv, nm[i], txt, sep="")
  }
  cv <- paste(cv,", c(", nm[n.vars-1], ", ", nm[n.vars],"))", sep="")

  if (results && n.pred > 2) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- cv 
    tx[length(tx)+1] <- "```"
  }

  if (explain && n.pred > 2) tx[length(tx)+1] <- paste("\n",
"First verify that the reduced and full models are properly specified. ",
sep="")

  if (explain && n.pred > 2) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"n$out_models" 
    tx[length(tx)+1] <- "```"
  }


  if (explain && n.pred > 2) tx[length(tx)+1] <- paste("\n",
"Compare nested models to evaluate ",
"if the additional variables in the full model provide ",
"a detectable increase in fit beyond that of the reduced model. ",
"Evidence for accepting the reduced model is to have the test ",
"_not_ be significant, which means that the evaluated ",
"coefficients from the additional variables in the full model ",
"are not detected to be different from 0, and so ",
"perhaps need not be included in the model. ", 
sep="")

  if (explain && n.pred > 2) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <-"n$out_anova" 
    tx[length(tx)+1] <- "```"
  }

  tx[length(tx)+1] <- "`r reject <- \"Reject the null hypothesis of the tested regression coefficients equal to 0 because of the small _p_-value of\"`"

  tx[length(tx)+1] <- "`r accept <- \"No difference of the coefficients from zero detected because of the relatively large _p_-value of\"`"

  if (explain && n.pred > 2) tx[length(tx)+1] <- paste("\n",
"`r if ((n$anova_tested[5]< 0.05)) reject else accept` ",
"`r xP(n$anova_tested[5],3)`. ", 
sep="")

  if (explain && n.pred > 2)tx[length(tx)+1] <- paste(
"Realize that if the reduced model was constructed solely from ", 
"the regression output of the initial training data, and then analyzed ",
"from the same data, ",
"the analysis is post-hoc. The _p_-value in this situation ",
"is an interesting and perhaps useful heuristic, but cannot be literally ",
" interpreted. ",
sep="")









  tx[length(tx)+1] <- ""
  if (n.pred > 1)
    tx[length(tx)+1] <- "## Relations Among the Variables"
  else
    tx[length(tx)+1] <- paste("## Relation Between ", Y, " and ", X, sep="")

  if (numeric.all) {      

    tx[length(tx)+1] <- ""
    if (n.pred > 1) tx[length(tx)+1] <- "### Scatter Plot"

    if (explain) tx[length(tx)+1] <- paste("\n",
"How do the variables in the model relate to each other? ",
"The correlation", pl, " of response variable ", Y, " with predictor ",
"variable", pl, " ", X, " should be high. ",
sep="")

    if (explain  &&  n.pred > 1) tx[length(tx)+1] <- paste(
"The correlations of the predictor variables ",
"with each other should be relatively small. ",
sep="")

      if (n.pred == 1) tx[length(tx)+1] <- paste(
"The correlation of ", Y, " with ", X, " in the training data is $r$ = ",
"`r xP(r$cor[2,1],3)`. ",
sep="")


      if (n.pred > 1)
        txt <- "all the relationships among the variables"
      else
        txt <- paste("the relationship of ", Y, " and ", X, sep="")
    
      if (explain) tx[length(tx)+1] <- paste("\n",
"Visually summarize ", txt, " in the model ", 
"with the scatterplot. ", 
sep="") 

      if (explain  &&  n.pred > 1) tx[length(tx)+1] <- paste(
"The model has multiple predictor variables, so the ",
"different scatter plots are presented in a scatter plot matrix. ",
"Each scatter plot in the matrix also contains a non-linear best-fit ",
"curve. ",
sep="")

    if (results && n.pred > 1) tx[length(tx)+1] <- paste(
"Express the linear numeric variable relationships among the variables ",
"in the model with their correlations. ",
sep="")

    if (document) tx[length(tx)+1] <- paste(
"Plot the scatter plot separately with the ",
"`lessR` function `regPlot`. Specify option 1 to indicate this specific plot.",
sep="")

    if (results) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      if (n.pred > 1)
        tx[length(tx)+1] <- "regPlot(r, 1)  # 1: scatter plot matrix"
      else
        tx[length(tx)+1] <- "regPlot(r, 1, pred.intervals=FALSE)  # 1: scatter plot "
      tx[length(tx)+1] <- "```"
    }


    if (n.pred >1) {

      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- "### Collinearity Analysis"

      if (explain) tx[length(tx)+1] <- paste(
"The collinearity analysis assesses the extent that the ",
"predictor variables -- ", X, " -- linearly depend upon each other, ",
"which in the simplest case is a high pairwise correlation. ",
"Although collinearity diminishes neither the fit of the model, ",
"nor forecasting efficacy, it ",
"typically indicates an overly complex model. ",
"Collinear variables have relatively large ",
"standard errors of their estimated slope coefficients, ",
"which yield unstable estimates. The result is that any ",
"unique effects of collinear variables ",
"cannot be statistically disentangled without a very ",
"large sample size. ",
sep="")

      if (explain) tx[length(tx)+1] <- paste("\n",
"To assess collinearity for predictor variable $X_j$, regress ",
"that predictor onto all of the remaining predictor variables. A high ",
"resulting $R^2_j$ ",
"indicates collinearity for that predictor. ",
"Usually express this result in terms of ",
"the _tolerance_ of the predictor, $1 - R^2_j$, ",
"the proportion of variance for $X_j$ _not_ due do the ",
"remaining predictor variables. ",
"Because each $R^2_j$ should be low, presumably at least less than 0.8, ",
"tolerance should be high, at least larger ",
"than approximately 0.20. ",
"An equivalent statistic is the _variance inflation factor_ (VIF), ",
"which indicates the extent collinearity inflates the variance of ",
"the estimated coefficient. VIF is the reciprocal of tolerance, so ",
"usually should be at least less than approximately 5.",
sep="")

      if (results) {
        tx[length(tx)+1] <- ""
        tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
        tx[length(tx)+1] <-"r$out_collinear" 
        tx[length(tx)+1] <- "```"
      }

      l20 <- length(which(tolerances < 0.20))
      if (interpret  &&  l20 > 1)
        hh <- "s have tolerances"
      else
        hh <- " has a tolerance"

      if (interpret && l20 > 0)
        tx[length(tx)+1] <- paste(
"Collinearity is indicated. ", 
xU(xNum(l20)), " variable", hh, " less than the ", 
"cutoff of 0.20: ",  # r$tolerances has no names, so get them at r$model
"`r xAnd(names(r$model[which(r$tolerances < 0.20) + 1]))`. ",
sep="")

      l2030 <- length(which(tolerances >= 0.20 & tolerances < 0.30))
      if (l2030 > 1)
        hh <- "s have tolerances"
      else
        hh <- " has a tolerance"

      if (interpret  &&  l2030 > 0) 
        tx[length(tx)+1] <- paste(
xU(xNum(l2030)), " variable", hh, " greater than the ", 
"cutoff of 0.20, but still somewhat low, less than 0.30: ",
"`r xAnd(names(r$model[which(r$tolerances >= 0.20 & r$tolerances < 0.30) + 1]))`. ",
sep="")

      if (interpret && length(which(tolerances < 0.30)) == 0)
        tx[length(tx)+1] <- paste(
"No collinearity exists according to the tolerance cutoff of 0.30. ",
sep="")

      pl2 <- ifelse (length(which(tolerances == min(tolerances))) > 1, "s", "")
      vrb <- ifelse (length(which(tolerances == min(tolerances))) > 1,
                     "are", "is")

      if (results  &&  n.pred > 1) tx[length(tx)+1] <- paste(
"The predictor variable", pl2, " with the lowest tolerance ", vrb, " ",
"`r xAnd(names(r$model[which(r$tolerances == min(r$tolerances)) + 1]))` at ",
"`r xP(min(r$tolerances),3)`.",
sep="")



      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- "### All Possible Subset Models"

      if (explain) tx[length(tx)+1] <- paste(
"Especially when collinearity is present, can a simpler model be ",
"more, or at least almost, as effective as the current model? ",
sep="")

      if (explain) tx[length(tx)+1] <- paste("\n",
"To investigate, assess the fit for the models that ",
"correspond to all possible combinations of the predictor variables. ",
"Each row of the analysis defines a different model. ",
"A 1 means the predictor variable is in the model, a 0 ",
"means it is excluded from the model.",
sep="")

      if (results) {
        tx[length(tx)+1] <- ""
        tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
        tx[length(tx)+1] <-"r$out_subsets" 
        tx[length(tx)+1] <- "```"
      }

      if (explain) tx[length(tx)+1] <- paste("\n",
"The goal of this analysis is _parsimony_, to obtain the ",
"most explanatory power, here assessed with $R^2_{adj}$, with ",
"the least number of predictor variables, presumably guided also ",
"by the content and meaningfulness of the variables in the model. ",
sep="")

      if (explain) tx[length(tx)+1] <- paste("\n",
"Note that this analysis only describes the available data. ",
"This subset analysis is a ",
"descriptive heuristic that can effectively help eliminate unnecessary ",
"predictor variables from your model, but all resulting inferential ",
"statistics such as _p_-values are no longer valid. ",
"Ultimately a model revised from the training data ",
"requires cross-validation on a new data set. ",
sep="")

    }  # end n.pred > 1

  }  # end numeric.all

  else
     tx[length(tx)+1] <- 
       "No analysis due to non-numeric variables"






  
  if (res.rows >0) {

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- "## Analysis of Residuals and Influence"

    if (explain) tx[length(tx)+1] <- paste(
"Values of ", Y, " fitted by the estimated model do not generally ",
"equal the corresponding data values. Which cases (rows of data) ",
"contribute the most to this lack of fit? ",
sep="")

    if (explain) tx[length(tx)+1] <- paste(
"The identification of cases that have a large residual ",
"and/or undue influence on the estimation of the model helps ",
"detect potential outliers. For each case, ",
"in addition to the data values, fitted value and corresponding residual, ",
"the analysis provides the following values . ",
sep="")

    if (explain) tx[length(tx)+1] <- paste("\n",
"* _residual_: Value of the response variable ", Y, " minus its ",
"fitted value, $e = Y_{", Y, "} - \\hat Y_{", Y, "}$ \n",
"* _rstudent_: Externally Studentized residual, standardized value of the ",
"residual from a model estimated without the case present \n",
"* _dffits_: Standardized difference between a fitted value with and without ",
"the case present \n",
"* _cooks_: Cook's Distance, the aggregate influence of the case ",
"on all the fitted values with each fitted value calculated with the case deleted \n",
sep="")

    if (results) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      tx[length(tx)+1] <-"r$out_residuals" 
      tx[length(tx)+1] <- "```"
    }

    if (res.sort != "off") {
      cv[1] <- paste("`r xP(r$resid.max[1],2)`", sep="")
      for (i in 2:4)
        cv <- paste(cv, ", `r xP(r$resid.max[", i, "],2)`", sep="")
      cv <- paste(cv, " and `r xP(r$resid.max[5],2)`.", sep="")
      

      if (res.sort == "cooks") txt <- "Cook's distances"
      else if (res.sort == "rstudent") txt <- "Studentized residuals" 
      else if (res.sort == "dffits") txt <- "dffits values" 
 
      if (results) tx[length(tx)+1] <- paste(
"From this analysis the five largest ", txt, ": ", cv,
sep="")

      if (interpret  && res.sort == "cooks") tx[length(tx)+1] <- paste("\n",
"An informal criterion ",
"for evaluating the size of Cook\'s distance is a cutoff value of 1 ",
"to indicate too large of a large size. ",
sep="")

    if (length(which(resid.max > 1)) > 1)
      hh <- "s have values"
    else
      hh <- " has a value"
  
    lbl <- xRow(resid.max)
    if (interpret  && res.sort == "cooks"  && (length(which(resid.max > 1)) > 0))
      tx[length(tx)+1] <- paste(
"The following case", hh, " more than the ", 
"cutoff of 1: ", xAnd(lbl[which(resid.max > 1)]), ". ",
"For larger sample sizes this guideline should be reduced as the ",
"influence of any case tends to diminish as the sample size increases. ",
"The best basis for understanding a high Cook's distance value is to ",
"also consider the substantive nature of the underlying data values ",
"of the case and to verify if they were sampled from the same population ",
"as the remaining cases. ",
sep="")


    if (interpret  && res.sort == "cooks"  && (length(which(resid.max > 1)) == 0))
      tx[length(tx)+1] <- paste(
"No cases have a Cook's distance larger than 1 in this analysis. ", 
sep="")


    }  # res.sort is on

  }  # res.rows > 0



  if (pred.rows > 0) {

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- "## Prediction Intervals"

    if (explain) tx[length(tx)+1] <- paste(
"Ultimately the analysis moves beyond the training sample. ",
"Prediction is from _new_ data values of ", X, ", what may be ",
"called the _prediction sample_. Applying these new data values ",
"to the estimated model yields the predicted values. ",
"For data values from the training sample, the fitted value and ",
"predicted value are the same, ",
"calculated from the same estimated model, but are ",
"different concepts with different interpretations. ",
sep="")
    
    if (explain) tx[length(tx)+1] <- paste("\n",
"Unfortunately, prediction is not perfect. ",
"The range of values likely to contain the ",
"actual data value for ", Y, " predicted from specific values of ", X, " ",
"quantifies the _forecasting error_. ",
"The standard deviation of the residuals, $s_e$, assumed to be the same ",
"for all sets of values of the predictor variables, specifies the ",
"_modeling error_ of the fitted values from the training data, error ",
"due to imperfections in the model. ",
"However, for predictions of future values of ", Y, ", new data ",
"are collected. So sampling ",
"error of a value on the regression line, ",
"indicated with $s_{\\hat Y}$, must also be considered in the ",
"assessment of forecasting error. Consideration of both sources of ",
"error results in the _standard error of forecast_. ",
sep="")

    if (explain) tx[length(tx)+1] <- paste("\n",
"$$s_f = \\sqrt{s^2_e + s^2_{\\hat Y}}$$",
sep="")

    if (explain) tx[length(tx)+1] <- paste("\n",
"Unlike modeling error, the amount of ",
"sampling error varies ",
"depending on the values of the predictor variables, so each ",
"row of data has its own value, $s_f$. ",
"Each prediction interval ",
"is the margin of error, the _t_-cutoff multiplied by the corresponding ",
"$s_f$, added and subtracted on either side ",
"of $\\hat Y_{", Y, "}$.",
sep="")




  if (n.pred <= 6  &&  numeric.all  &&  is.null(X1.new)) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- "### Prediction Intervals from the Training Data"
  }

    if (results) tx[length(tx)+1] <- paste("\n",
"The analysis provides each row of data values, _as if_ they were  ",
"new data, with a predicted value based ",
"on the model estimated ",
"from the training data, as well as the standard error of forecast. From ",
"these values obtain the lower and upper bounds of the corresponding ",
"95% prediction interval. By default, only the first three, middle three ",
"and last three rows of data are presented, sufficient to indicate the ",
"ranges of prediction error encountered throughout the ranges of data values",
sep="")

    if (results) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      tx[length(tx)+1] <-"r$out_predict" 
      tx[length(tx)+1] <- "```"
    }

    if (interpret) tx[length(tx)+1] <- paste("\n",
"The size of the prediction intervals for the range of data found ",
"in the input data table vary from a minimum of ",
"`r xP(r$pred_min_max[1], ", d, ",", uYq, ")` for `r xAnd(xRow(r$pred_min_max[1]))` ",
"to a maximum of ",
"`r xP(r$pred_min_max[2], ", d, ",", uYq, ")` for `r xAnd(xRow(r$pred_min_max[2]))`. ",
sep="")

    if (explain  &&  n.pred == 1)
      tx[length(tx)+1] <- paste("\n",
"The confidence intervals for the points on the regression line, ",
"and the much larger prediction intervals for the individual data points, ",
"are illustrated with an enhancement of the original scatter plot.",
sep="")

    if (document) tx[length(tx)+1] <- paste(
"Plot the scatter plot with prediction intervals separately with the ",
"`lessR` function `regPlot`. Specify option 1 to indicate this specific plot.",
sep="")

    if (results  &&  n.pred == 1) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      tx[length(tx)+1] <- paste(
"regPlot(r, 1)  # 1: scatter plot with prediction intervals",
sep="")
      tx[length(tx)+1] <- "```"
    }



  if (n.pred <= 6  &&  numeric.all  &&  is.null(X1.new)) {

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- "### Prediction Intervals from New Data"

    tx[length(tx)+1] <- paste("\n",
"New data values from which to obtain ",
"a forecast, different from the training data, can be entered with ",
"the options X1.new, X2.new, up to X6.new, where each option name refers to ",
"the position of the corresponding predictor variable in the specification ",
"of the regression model. ", 
"Any number of values can be specified for each predictor variable. ",
"Suppose, for example, that there are two values of interest ",
" for ", et, "predictor variable ",
"from which to make a forecast, listed below. ",
sep="")
    
    tx[length(tx)+1] <- ""
    for (i in 1:n.pred) {
      tx[length(tx)+1] <- paste(
pred[i], ": ", new.val[i,1], ", ", new.val[i,2], "  ",  
sep="")
    }

    tx[length(tx)+1] <- paste("\n",
"Re-run the analysis to obtain the prediction intervals with these ",
"specified values. ",
sep="")

    tx[length(tx)+1] <- ""
    cv <- ",\n        "
    for (i in 1:n.pred)
      cv <- paste(cv,
" X", i, ".new=c(", new.val[i,1], ",", new.val[i,2], ")",
ifelse(i == n.pred, "", ","), 
sep="")
    cv <- paste(cv, ",\n         graphics = FALSE", sep="")
    fc <- sub(", graphics = FALSE", cv, fc, fixed=TRUE)

    if (results) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      tx[length(tx)+1] <- paste("r <-", fc)
      tx[length(tx)+1] <- "```"
    }

    if (n.pred > 1)
      tx[length(tx)+1] <- paste("\n",
"The new data values are specified for each variable separately, but ",
"a row of data consists of data values for all the predictor values. ",
"Accordingly, calculate a prediction interval for each combination ",
"of the specified new values for each predictor variable. ",
sep="")
    else
      tx[length(tx)+1] <- paste("\n",
"Calculate the prediction intervals only for the new data values. ",
sep="")

    if (results) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
      tx[length(tx)+1] <-"r$out_predict" 
      tx[length(tx)+1] <- "```"
    }

    if (interpret) tx[length(tx)+1] <- paste("\n",
"The size of the prediction intervals for the range of data found ",
"in the newly specified values vary from a minimum of ",
"`r xP(r$pred_min_max[1], ", d, ",", uYq, ")` for `r xAnd(xRow(r$pred_min_max[1]))` ",
"to a maximum of ",
"`r xP(r$pred_min_max[2], ", d, ",", uYq, ")` for `r xAnd(xRow(r$pred_min_max[2]))`. ",
"The rows in the output display, however, are re-ordered according to the ",
"combinations of the ordered values of the predictor variables. ",
sep="")

    }

  }







  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "## Model Validity"

  if (explain) tx[length(tx)+1] <- paste(
"The residuals should be independent, normal random variables with a ",
"mean of zero and constant variance. ",
sep="")

  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Distribution of Residuals"

  if (explain) tx[length(tx)+1] <- paste(
"For the inferential analyses to be valid, ",
"the residuals should be normally distributed. ",
sep="")

  if (explain) {
    tx[length(tx)+1] <- paste(
"Violation of normality does not bias the estimates, but ",
"it does render the inferential tests invalid. ",
sep="")
  }

    if (document) tx[length(tx)+1] <- paste(
"Plot the distribution of residuals separately with the `lessR` ",
"function `regPlot`. Specify option 2 to indicate this specific plot.",
sep="")

  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- "regPlot(r, 2)  # 2: distribution of residuals"
    tx[length(tx)+1] <- "```"
    tx[length(tx)+1] <- ""
  }


  tx[length(tx)+1] <- ""
  tx[length(tx)+1] <- "### Fitted Values vs Residuals"

  if (explain) tx[length(tx)+1] <- paste(
"The residuals should represent random variation, free of ",
"any pattern or structure. ",
"They should satisfy the _homoscedasticity_ assumption, ",
"randomly scattered about 0, with approximately ",
"the same level of variability across the range of the fitted values ",
"within a horizontal band around the zero-line. Otherwise they ",
"demonstrate _heteroskedasticity_. ",
sep="")

    if (document) tx[length(tx)+1] <- paste(
"Plot the scatter plot of fitted values with residuals separately with the ",
"`lessR` function `regPlot`. Specify option 3 to indicate this specific plot.",
sep="")

  
  if (results) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("```{r", show, "}", sep="")
    tx[length(tx)+1] <- "regPlot(r, 3)  # 3: scatter plot of fitted with residuals"
    tx[length(tx)+1] <- "```"
  }

  return(tx)

}
