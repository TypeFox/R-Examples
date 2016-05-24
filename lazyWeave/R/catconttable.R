#' @name ComparisonTable
#' @export catconttable
#' @importFrom Hmisc label.default
#' @importFrom Hmisc label.data.frame
#' @importFrom Hmisc 'label<-.default'
#' @importFrom Hmisc 'label<-.data.frame'
#' @importFrom Hmisc print.labelled
#' @importFrom Hmisc '[.labelled'
#' @importFrom Hmisc latexTranslate
#' @importFrom Hmisc smean.cl.boot
#' 
#' @title Comparison Tables
#' @description Produce a table of comparisons for reports and manuscripts
#' 
#' @param data A \code{ccf.df} or \code{data.frame} with the variables to be
#'    compared and the grouping variable.
#' @param vars A character vector naming the variables to be compared
#' @param vars.cat A character vector that can be used to specify which, if
#'    any, numeric variables in \code{vars} should be treated as categorical.
#' @param byVar A character(1) giving the grouping variable.  This allows more
#'    than one level.  Numeric variables are coerced to factors.
#' @param fisher A character vector giving the names of variables that should
#'    be compared with Fisher's Exact test.  Currently, there is no
#'    implementation to determine this automatically.
#' @param fisher.arg A list of additional arguments to pass to \code{fisher.test}
#' @param cmh A character vector giving the names of variables that should
#'    be compared with Manthel-Haenszel's Test for Linear Trend.  This is 
#'    not yet written and will be ignored.
#' @param row.score Currently ignored
#' @param col.score Currently ignored
#' @param mcnemar a character vector giving the names of variables that should
#'    be compared using McNemar's test.
#' @param correct Character vector giving the variables for which a continuity
#'    correction should be applied in McNemar's test.
#' @param odds A character vector giving the names of variables for which
#'    odds ratios should be calculated.  For categorical measures, this is
#'    the primary test of comparison.  For numeric measures, this is calculated
#'    in addition to another test.
#' @param odds.scale For numeric variables only.  A list with named elements
#'    that gives the scale on which the odds ratio should be presented.  For
#'    example, if the odds for variable \code{x} should be presented in 5 year
#'    increments, we would use \code{odds.scale=list(x = 5)}.
#' @param odds.unit For numeric variables only.  A list with named elements 
#'    that gives the units on which the odds ratio should be presented.  For
#'    example, if the odds of variable \code{x} should be presented in 5 year
#'    increments, we would use \code{odds.unit=list(x="years")}.
#' @param none A character vector naming variables in \code{vars} for which no 
#'    comparison should be made.
#' @param row.p Toggles if row or column proportions are calculated.
#' @param normal A character vector that assigns variables in \code{vars} as
#'    normally distributed.
#' @param var.equal A character vector that assigns variables in \code{vars} as
#'    having equal variance.  This is used to determine the proper form of
#'    a t-test.
#' @param median A character vector that assigns variables that shoudl be 
#'    summarized with a median, quartiles, or min and max.  
#' @param alpha Significance levels for tests.
#' @param B The number of Bootstrap samples for bootstrapped confidence 
#'    intervals.
#' @param seed The seed to use in starting the Bootstrapping.
#' @param minl Minimum length for levels abbreviations.  The function
#'    \code{abbreviate} is used to create unique rownames for each level of 
#'    a variable in the output data frame.  If the abbreviations are short, 
#'    they may not be readable.  This allows the user to make the length longer.
#'    
#' @details   \code{catconttable} is a wrapper that determines the type of 
#' variable and calls either cattable or conttable as appropriate.  For this 
#' to work properly, all factor variables must be defined before the function 
#' call.
#' 
#' In contrast, if cattable is called directly, variables are coerced to 
#' factors, which could lead to peculiar results if a numeric value is given.
#' 
#' @author Benjamin Nutter
#' 
#' @seealso \code{\link{write.ctable}}
#' 
#' @examples
#'
#' #Read in the delivery dataset from the lazyWeave package
#' data(Delivery)
#' 
#' #Use conttable to summarize maternal age, ga weeks, weight (grams) 
#' #and grava by delivery type.  The dataset name is specified under the "data="
#' #option, the variables of interest are listed under "vars=", and the K-level by variable 
#' #is specified under "byVar=".
#' 
#' #Default is to report mean and bootstrapped 95% CI for mean.  Tests of location are by 
#' #default either Wilcoxon rank sum (K=2) or Kruskal-Wallis (K>2) rank sum.  The "seed="
#' #option allows for reproducibility by setting the seed for getting bootstrapped samples.
#' 
#' d_type.contable <- conttable(data=Delivery,
#'                              vars=c("maternal.age", "ga.weeks", "wt.gram", "grava"),
#'                                     byVar="delivery.type")
#' 
#' #Specifying weights by delivery type as a normally distributed variables, reports means, 
#' #standard deviations and a t-test of equality of the means for delivery type.  Variables listed 
#' #under "var.equal=" are assumed to have equal variances in all levels of byVar.  Otherwise, 
#' #variances are allowed to be unequal.
#' 
#' d_type.conttable <- conttable(data=Delivery,
#'                               vars=c("maternal.age", "ga.weeks", "wt.gram", "grava", "apgar1"),
#'                               byVar="delivery.type",
#'                               normal=c("wt.gram", "maternal.age"),
#'                               var.equal="ga.weeks")
#'                               
#' #List variables under "median=" to report median, 25th and 75th percentiles.
#' d_type.conttable <- conttable(data=Delivery,
#'                               vars=c("maternal.age", "ga.weeks", "wt.gram", "grava", "apgar1"),
#'                               byVar="delivery.type",
#'                               normal=c("wt.gram", "maternal.age"),
#'                               var.equal="ga.weeks",
#'                               median=c("grava","apgar1"))
#' 
#' #Use cattable to summarize child sex, laceration, and laceration degree by delivery type.
#' #Row percent, overall counts, and counts by delivery type are reported.  Column percents can 
#' #be specified by the row.p=FALSE option.
#' #By default chi-square tests of independence are performed.
#' 
#' d_type.cattable <- cattable(data=Delivery,
#'                             vars=c("child.sex", "laceration"),
#'                             byVar="delivery.type")
#' 
#' #For variables listed under "fisher=" Fisher's exact test of independence is performed.
#' #The reported test statistic is the odds ratio.
#' 
#' d_type.cattable <- cattable(data=Delivery,
#'                             vars=c("child.sex", "laceration"),
#'                             fisher=c("child.sex"),
#'                             byVar="delivery.type")
#' 
#' 
#' #All variables listed in a single table
#' 
#' d_type.catconttable <- catconttable(data=Delivery,
#'                                     vars=c("maternal.age", "ga.weeks", "child.sex", "wt.gram",
#'                                            "grava", "apgar1", "laceration"),
#'                                     median=c("grava", "apgar1"),
#'                                     normal="maternal.age",
#'                                     fisher="child.sex",
#'                                     byVar="delivery.type")
#' 
#' \dontrun{
#'   #Code for writing ctable objects to a file.  See write.ctable() for more information
#'   
#'   #Write to PDF
#'   options(lazyReportFormat='latex')
#'   lazy.write(
#'     lazy.file.start(),
#'     write.ctable(d_type.catconttable),
#'     lazy.file.end(),
#'     OutFile="SampleOutput.tex")
#'     
#' #Generate a pdf in the working directory
#'   lazy.build("SampleOutput.tex")
#'   
#'   unlink("SampleOutput.tex")
#'   unlink("SampleOutput.pdf")
#' } 
#' 
catconttable <- function(data, vars, byVar, vars.cat=NULL, fisher=NULL, fisher.arg=NULL,
                     cmh=NULL, row.score=NULL, col.score=NULL,
                     normal = NULL, var.equal = NULL, 
                     median=NULL, odds=NULL, odds.scale=NULL, odds.unit=NULL,
                     none=NULL,
                     row.p=TRUE, alpha=0.05, B=1000, seed=NULL){
  
  if (missing(byVar)){
    byVar <- "PlAcE_hOlDeR_fOr_CaTcOnTtAbLe"
    data[, byVar] <- factor("")
  }
  
  if (!all(vars %in% names(data))){
    bad.vars <- vars[!vars %in% names(data)]
    bad.vars.msg <- paste("The following variables are not found in 'data':", paste(bad.vars, collapse=", "))
    stop(bad.vars.msg)
  }

  all.missing <- sapply(data[, c(vars, byVar)], function(x) all(is.na(x)))
  if (any(all.missing)){
    miss.vars <- c(vars, byVar)[all.missing]
    miss.vars.msg <- paste("The following variables contain only missing values:", paste(miss.vars, collapse=", "))
    stop(miss.vars.msg)
  }
  
  if ("tbl_df" %in% class(data)) data <- as.data.frame(data)
  
  var.info <- function(v, ...){
    if (!is.numeric(data[, v]) | v %in% vars.cat)
      cattable(data=data, vars=v, byVar=byVar, fisher=fisher, fisher.arg=fisher.arg,
                    cmh=cmh, row.score=row.score, col.score=col.score,
                    odds=odds,
                    none=none, row.p=row.p, alpha=0.05)
    else conttable(data=data, vars=v, byVar=byVar,
                 normal = normal, var.equal = var.equal, median=median,
                 odds = odds, odds.scale=odds.scale, odds.unit=odds.unit,
                 alpha = alpha, B=B, seed=seed)
  }

  ctable <- do.call("rbind", lapply(vars, var.info))
  ctable$type <- factor(ctable$type)
  attributes(ctable)$byVar <- data[, byVar]
  Hmisc::label(attributes(ctable)$byVar) <- Hmisc::label(data[, byVar])
  attributes(ctable)$vars <- vars  
  return(ctable)
}

             
             
