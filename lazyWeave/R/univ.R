#' @name univ
#' @export univ
#' @importFrom Hmisc label.default
#' @importFrom Hmisc label.data.frame
#' @importFrom Hmisc 'label<-.default'
#' @importFrom Hmisc 'label<-.data.frame'
#' @importFrom Hmisc print.labelled
#' @importFrom Hmisc '[.labelled'
#' @importFrom Hmisc latexTranslate
#' 
#' @title Univariable Table
#' @description Similar to the QHS SAS macro, provides a simple
#'   summary of numerical variables
#'   
#' @param data A \code{ccf.df} or \code{data.frame} containing the variables in
#'   \code{vars} and \code{byVar}.
#' @param vars Character vector naming the variables in \code{data} to be
#'   summarized.
#' @param byVar A categorical variables giving groups by which statistics 
#'   should be computed.
#' @param alpha significance level for determining the confidence limits.
#' @param test The test to use for comparing between groups.  Currently 
#'   limited to t-test and ANOVA (parametric tests).
#' @param test.args a list of additional arguments to pass to the function in 
#' \code{test}.
#' 
#' @details
#' Statistics available in univ, and the calls to get them are:
#' \enumerate{
#' \item{\code{n}}{number of non-missing values.}
#' \item{\code{missing}}{number of missing values}
#' \item{\code{mean}}{arithmetic mean}
#' \item{\code{median}}{median value}
#' \item{\code{sd}}{standard deviation}
#' \item{\code{lcl}}{lower confidence limit}
#' \item{\code{ucl}}{upper confidence limit}
#' \item{\code{min}}{minimum value}
#' \item{\code{max}}{maximum value}
#' \item{\code{p25}}{25th percentile}
#' \item{\code{p75}}{75th percentile}
#' \item{\code{cv}}{coefficient of variation}
#' }
#' 
#' \code{univ} does not perform any kind of inference on the varaibles 
#' supplied in the argument.  It returns only summary statistics.  If 
#' comparisons are desired, try using \code{conttable}.  Confidence limits 
#' are determined using the t-distribution.
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' data(Delivery)
#' #Read in the delivery dataset from the CCFmisc library
#' #use Hmisc library to labeling variables in univariate tables
#' Hmisc::label(Delivery$maternal.age) <- "Maternal Age"
#' Hmisc::label(Delivery$ga.weeks) <- "Gestation weeks"
#' Hmisc::label(Delivery$wt.gram) <- "Weight (g)"
#' 
#' 
#' #a univariate table of the variables maternal age,
#' #ga.weeks and wt.grams.  The object resulting
#' #from univ() can be used in other functions to create html or
#' #LaTeX tables.
#' 
#' uni <- univ(Delivery,
#'             vars=c("maternal.age", "ga.weeks", "wt.gram"))
#' 
#' #a univariate table of the variables maternal age,
#' #ga.weeks and wt.grams by delivery.type.  The object resulting
#' #from univ() can be used in other functions to create html or
#' #LaTeX tables.
#' 
#' deliv.uni <- univ(Delivery,
#'                   vars=c("maternal.age", "ga.weeks", "wt.gram"),
#'                   byVar="delivery.type")
#' 
#' #if you want to take advantage of the confidence interval
#' #output from univ() different alpha levels can be set
#' #by the alpha= argument.
#' 
#' deliv_99.uni <- univ(Delivery,
#'                      vars=c("maternal.age", "ga.weeks", "wt.gram"),
#'                      byVar="delivery.type",
#'                      alpha=0.01)
#'

'univ' <- function(data, vars, byVar, alpha=0.05,
                   test=c("t.test", "aov", "wilcox.test", "kruskal.test"), test.args=NULL){

#********************************************************************
#* 3. Provide dummy byVar if necessary
#* 4. Coerce byVar to a factor variable
#********************************************************************

#*** 3. Provide dummy byVar if necessary
  if(missing(byVar)){
    byVar <- "SomePlaceHolderVariable"
    data[,byVar] <- factor(1)  
    byVar.miss <- TRUE }
  else byVar.miss <- FALSE

#*** 4. Coerce byVar to a factor variable
  if(!is.factor(data[,byVar])){
    data[,byVar] <- factor(data[,byVar])
    warning(paste(expression(byVar),"was coerced to a factor variable"))  }
  
  test <- match.arg(test, c("t.test", "aov", "wilcox.test", "kruskal.test"))
 
#********************************************************************
#* 1. Generic function to obtain statistics
#* 2. Generic function to obtain counts
#********************************************************************

#*** 1. Generic function to obtain statistics
  stat.func <- function(v,func,...){
    st <- suppressWarnings(tapply(data[,v],data[,byVar],func,na.rm=TRUE,...))
    st <- ifelse(is.finite(st), st, NA)
    return(st)
  }

#*** 2. Generic function to obtain counts
  count.func <- function(v,missing=FALSE){
    if(missing) tapply(is.na(data[,v]),data[,byVar],sum)
    else        tapply(!is.na(data[,v]),data[,byVar],sum)  }

#********************************************************************
#* Get Factor and Group names and corresponding statistics
#********************************************************************

  if(byVar.miss) Factor <- Hmisc::label(data[vars], default=names(data[, vars]))
  else Factor <- unlist(lapply(vars,function(v) 
                    c(Hmisc::label(data[, v], default=v),
                      rep(NA,nlevels(data[,byVar])-1))))
  Group <- rep(levels(data[,byVar]),length(vars))
  N <- as.vector(sapply(vars,count.func))
  MISSING <- as.vector(sapply(vars,count.func,missing=TRUE))
  MEAN <- as.vector(sapply(vars,stat.func,func="mean"))
  MEAN <- ifelse(is.nan(MEAN), NA, MEAN)
  SD <- as.vector(sapply(vars,stat.func,func="sd"))
  LCL <- suppressWarnings(ifelse(N > 1, MEAN - stats::qt(1-alpha/2, N-1) * SD, NA))
  UCL <- suppressWarnings(ifelse(N > 1, MEAN + stats::qt(1-alpha/2, N-1) * SD, NA))
  MIN <- as.vector(sapply(vars,stat.func,func="min"))
  P25 <- as.vector(sapply(vars,stat.func,func="quantile",probs=0.25))
  MEDIAN <- as.vector(sapply(vars,stat.func,func="median"))
  P75 <- as.vector(sapply(vars,stat.func,func="quantile",probs=0.75))
  MAX <- as.vector(sapply(vars,stat.func,func="max"))
  CV <- SD/MEAN
 
  if (nlevels(data[, byVar]) > 1){  
    PVAL <- lapply(vars, function(v) do.call(test, c(list(data[, v] ~ data[, byVar]), test.args))) 
    PVAL <- if (test != "aov") unlist(lapply(PVAL, function(x) c(x$p.value, rep(NA, nlevels(data[, byVar]) - 1))))
            else if (test == "aov") unlist(lapply(PVAL, function(x) c(stats::anova(x)[1, 5], rep(NA, nlevels(data[, byVar]) - 1))))
  }
  else PVAL <- NA

#********************************************************************
#* 1. Prepare Output Data Frame
#* 2. Prepare Report Data Frame
#* 3. Format Non-Count Values in Report
#* 4. Remove Group Column if byVar is not provided
#********************************************************************

#*** 1. Prepare Output Data Frame
  output <- data.frame(Factor=Factor, Group=Group, N=N, MISSING=MISSING,
                       MEAN=MEAN, SD=SD, LCL=LCL, UCL=UCL,
                       MIN=MIN, P25=P25, MEDIAN=MEDIAN, P75=P75,
                       MAX=MAX, CV=CV, PVAL=PVAL, 
                       stringsAsFactors=FALSE)
  names(output) <- c("Factor", "Group", "N", "Missing", "Mean", "SD", "LCL",
                     "UCL", "Min", "P25", "MEDIAN", "P75", "MAX", "CV", "PVAL")

  rownames(output) <- NULL
  class(output) <- c("univ", "data.frame")
  return(output)
}
