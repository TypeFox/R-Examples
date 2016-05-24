#' @title Calculate SIR or SMR
#' @author Matti Rantanen, Joonas Miettinen
#' @description Poisson modelled standardised incidence or mortality ratios (SIRs / SMRs) i.e. 
#' indirect method for calculating standardised rates. SIR is a ratio of observed and expected cases.
#' Expected cases are derived by multiplying the strata-specific population rate with the
#' corresponding person-years of the cohort.
#' 
#' @details \code{sir} is a comprehensive tool for modelling SIRs/SMRs with flexible 
#' options to adjust and print SIR's, test homogeneity and utilize 
#' multistate data. The cohort data and the variable names for observation 
#' counts and person-years are required.
#' The reference data is optional, since the cohort data 
#' can be stratified (with \code{print}) and compared to total.
#' 
#' 
#' \strong{Adjust and print}
#' 
#' A SIR can be adjusted by the covariates found in both \code{coh.data} and \code{ref.data}.
#' Variable to adjust are given in \code{adjust}.
#' Variable names needs to match in both \code{coh.data} and \code{ref.data}. 
#' Typical variables to adjust by are gender, age group and calendar period.
#' 
#' \code{print} is used to stratify the SIR output. In other words, the variables 
#' assigned to \code{print} are the covariates of the Poisson model.
#' Variable levels are treaded as categorical.
#' Variables can be assigned in both \code{print} and \code{adjust}. 
#' This means the output it adjusted and printed by these variables.
#' 
#' \code{print} can also be a list of expressions. This allows changing variable 
#' names or transforming variables with functions such as \code{cut} and\code{round}.
#' For example, the existing variables \code{agegroup} and \code{year} could be
#' transformed to new levels using \code{cut} by
#' 
#' \code{print = list( age.category = cut(agegroup, breaks = c(seq(0,85,5), 120)), year.cat = cut(year, seq(1950,2015,10)))}
#' 
#' 
#' \strong{ref.rate or ref.obs & ref.pyrs}
#' 
#' The population rate variable can be given to the \code{ref.rate} parameter. 
#' That is, when using e.g. the \code{popmort} or a comparable data file, one may
#' supply \code{ref.rate} instead of \code{ref.obs} and \code{ref.pyrs}, which
#' will be ignored if \code{ref.rate} is supplied. 
#' 
#' 
#' Note that if all the stratifying variables in 
#' \code{ref.data} aren't listed in \code{adjust}, 
#' or when the categories are otherwise combined,
#' the (unweighted) mean of rates is used for computing expected cases.
#' This might incur a small bias in comparison to when exact numbers of observations
#' and person-years are available. 
#' 
#' 
#' 
#' \strong{mstate}
#' 
#' E.g. with \code{lexpand} it's possible to compute counts for several outcomes
#' so that the population at risk is same for each 
#' outcome such as a certain kind of cancer. 
#' The transition counts are in wide data format, 
#' and the relevant columns can be supplied to \code{sir}
#' in a vector via the \code{coh.obs} argument. 
#' The name of the corresponding new column in \code{ref.data} is given in
#' \code{mstate}. It's recommended to include the \code{mstate} variable in \code{adjust},
#' so the corresponding information should also be available in \code{ref.data}.
#' 
#' This approach is analogous to where SIRs are calculated separately their 
#' own function calls.
#' 
#' 
#' \strong{Other parameters}
#' 
#' The univariate multiple-comparison-adjusted p-value uses \code{\link[stats]{p.adjust}}.
#' Univariate confidence intervals are calculated using exact 
#' Poisson intervals (poisson.ci). The multivariate result
#' is based on a poisson regression model with profile-likelihood confidence intervals 
#' when possible. Otherwise Wald's normal-approximation is used.
#'
#' The p-value is a test for the levels of \code{print}. The test can be either 
#' \code{"homogeneity"}, a likelihood ratio test where the model with variable(s) in
#' \code{print} (categorical factor) is compared to the constant model. 
#' Option \code{"trend"} is the same likelihood ratio test except the 
#' variable(s) in \code{print} are/is continous.
#' 
#' 
#' \strong{EAR: Excess Absolute Risk}
#' 
#' A simple way to quantify the absolute difference between cohort risk and 
#' population risk.
#' Make sure that the person-years are calculated accordingly before using EAR.
#' 
#' Formula for EAR:
#' \deqn{EAR = \frac{observed - expected}{person years} \times 1000.}{EAR = (obs - exp)/pyrs * 1000.}
#' 
#' \strong{Data format}
#' 
#' The data should be given in aggregated format, i.e the number of observations 
#' and person-years are represented for each stratum.
#' The extra variables and levels are reduced automatically before estimating SIRs. 
#' Example of data format:
#' 
#' \tabular{rrrrr}{
#'   sex \tab age \tab period \tab obs \tab pyrs \cr
#'   0 \tab 1 \tab 2010 \tab 0 \tab 390 \cr
#'   0 \tab 2 \tab 2010 \tab 5 \tab 385 \cr
#'   1 \tab 1 \tab 2010 \tab 3 \tab 308 \cr
#'   1 \tab 2 \tab 2010 \tab 12 \tab 315
#' }
#' 
#' 
#' @param coh.data aggregated cohort data, see e.g. \code{\link{lexpand}}
#' @param coh.pyrs variable name for person years in cohort data; quoted or unquoted
#' @param coh.obs variable name for observed cases; quoted or unquoted
#' @param ref.data population data. Can be left NULL if \code{coh.data} is stratified in \code{print}.
#' @param ref.rate population rate variable (cases/person-years). Overwrites arguments
#' \code{ref.pyrs} and \code{ref.obs}; quoted or unquoted
#' @param ref.pyrs variable name for person-years in population data; quoted or unquoted
#' @param ref.obs variable name for observed cases; quoted or unquoted
#' @param subset logical condition to select data from \code{coh.data} before any computations
#' @param adjust variable names for adjusting without stratifying output; quoted vector or unquoted list
#' @param print variable names to stratify results; quoted vector or unquoted named list with functions
#' @param mstate set column names for cause specific observations; quoted or unquoted. Relevant only
#' when \code{coh.obs} length is two or more. See details.
#' @param test.type Test for equal SIRs. Test available are 'homogeneity' and 'trend'.
#' @param EAR logical; TRUE calculates Excess Absolute Risks for univarite SIRs.
#' (see details)
#' @param alpha level of type-I error in confidence intervals, default 0.05 is 95\% CI.
#' @param p.adj add multiple comparison p-value adjust type for univariate model, 
#' check \code{help(p.adjust)} for options. Default NULL doesn't add adjusted p-values.
#' @param round.by set number of digits in results
#' @param round.by.pvalue set number of digits in p-values

#' 
#' @examples 
#' data(popmort)
#' data(sire)
#' c <- lexpand( sire, status = status, birth = bi_date, exit = ex_date, entry = dg_date,
#'               breaks = list(per = 1950:2013, age = 1:100, fot = c(0,10,20,Inf)), 
#'               aggre = list(fot, agegroup = age, year = per, sex) )
#' ## SMR due other causes: status = 2
#' se <- sir( coh.data = c, coh.obs = 'from0to2', coh.pyrs = 'pyrs', 
#'            ref.data = popmort, ref.rate = 'haz', 
#'            adjust = c('agegroup', 'year', 'sex'), print = 'fot')
#' se
#' ## for examples see: vignette('sir')
#' 
#' 
#' @seealso \code{\link{plot.sir}}, \code{\link{lexpand}}
#' \href{../doc/sir.html}{A SIR calculation vignette}
#' 
#' @return A list of 5: 3 \code{data.table} objects, vector of starta variables and global p-value.
#' 
#' @export sir
#' 
#' @import data.table
#' @import stats




sir <- function( coh.data, 
                 coh.obs,
                 coh.pyrs,
                 ref.data = NULL,
                 ref.obs = NULL,
                 ref.pyrs = NULL, ref.rate = NULL,
                 subset = NULL,
                 print = NULL,
                 adjust = NULL,
                 mstate = NULL,
                 test.type = 'homogeneity',
                 alpha = 0.95,
                 p.adj = NULL,
                 EAR = FALSE,
                 round.by = 2,
                 round.by.pvalue = 4){

  coh.data <- data.table(coh.data)
  
  ## subsetting---------------------------------------------------------------
  ## no copy taken of data!
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data = coh.data, substiset = subset)
  coh.data <- coh.data[subset,]
  
  
  # print list --------------------------------------------------------------
  
  # env1 <- environment() # set environment where to assign new print
  # coh.data <- data_list(data = coh.data, arg.list = substitute(print), env = env1)
  
  mstate <- as.character(substitute(mstate))
  if(length(mstate) == 0) {
    mstate <- NULL
  }
  if(!is.null(mstate)) {
    coh.data[,(mstate) := 0L] 
  }
  
  # evalPopArg
  coh.obs <- substitute(coh.obs)
  c.obs <- evalPopArg(data = coh.data, arg = coh.obs)
  coh.obs <- names(c.obs)
  
  coh.pyrs <- substitute(coh.pyrs)
  c.pyr <- evalPopArg(data = coh.data, arg = coh.pyrs)
  coh.pyrs <- names(c.pyr)
  
  print <- substitute(print)
  c.pri <- evalPopArg(data = coh.data, arg = print)
  print <- names(c.pri)

  adjust <- substitute(adjust)
  c.adj <- evalPopArg(data = coh.data, arg = adjust)
  adjust <- names(c.adj)
  
  # collect data
  coh.data <- cbind(c.obs, c.pyr)
  if(!is.null(print))  coh.data <- cbind(coh.data, c.pri) 
  if(!is.null(adjust)) coh.data <- cbind(coh.data, c.adj)
  
  if( !is.null(ref.data) ){
    ref.obs <- as.character(substitute(ref.obs))
    ref.pyrs <- as.character(substitute(ref.pyrs))
    ref.rate <- as.character(substitute(ref.rate))
    
    if (length(ref.obs) == 0) ref.obs <- NULL
    if (length(ref.pyrs) == 0) ref.pyrs <- NULL
    if (length(ref.rate) == 0) ref.rate <- NULL
  }


  # print(coh.data)
  
  st <- sir_table( coh.data = coh.data, 
                   coh.obs = coh.obs,
                   coh.pyrs = coh.pyrs,
                   ref.data = ref.data,
                   ref.obs = ref.obs,
                   ref.pyrs = ref.pyrs, 
                   ref.rate = ref.rate,
                   print = print,
                   adjust = adjust,
                   mstate = mstate)
  
 
  results <- sir_est( table = st,
                      print = print,
                      adjust = adjust,                 
                      test.type = test.type,
                      alpha = alpha,
                      p.adj = p.adj,
                      EAR = EAR,
                      round.by = round.by,
                      round.by.pvalue = round.by.pvalue)
  
  # Output as data.frame if wanted  

  if (!getOption("popEpi.datatable")) {
    for (i in 1:3) {
      if (!is.null(results[[i]])) {
        setDFpe(results[[i]])
      }
    }  
  }
  
  setclass(results, c('sir', 'pe', class(results)))
  return(results)
}


#' @title Estimate splines for SIR or SMR
#' @author Matti Rantanen, Joonas Miettinen
#' 
#' @description Splines for standardised incidence or mortality ratio. A useful 
#' tool to e.g. check whether a constant SIR can be assumed for all calendar periods, a
#' gegroups or follow-up intervals. Splines can be fitted for these time dimensions
#' separately or in the same model.
#' 
#' @param coh.data cohort data with observations and at risk time variables
#' @param coh.pyrs variable name for person-years in cohort data
#' @param coh.obs variable name for observed cases
#' @param ref.data aggregated population data
#' @param ref.rate population rate observed/expected. This overwrites the parameters
#' \code{ref.pyrs} and \code{ref.obs}.
#' @param ref.pyrs variable name for person-years in population data
#' @param ref.obs variable name for observed cases
#' @param subset logical condition to subset \code{coh.data} before any computations
#' @param adjust variable names for adjusting the expected cases
#' @param print variable names for which to estimate SIRs/SMRs and 
#' associated splines separately 
#' @param mstate set column names for cause spesific observations. Relevant only
#' when coh.obs length is two or more. See help for \code{sir}.
#' @param spline variable name(s) for the splines
#' @param knots number knots (vector),  pre-defined knots (list of vectors) or for optimal number of knots left NULL
#' @param dependent.splines logical; if TRUE, all splines are fitted in same model.
#' @param reference.points fixed reference values for rate ratios. If left \code{NULL}
#' the smallest value is the reference point (where SIR = 1). 
#' Ignored if \code{dependent.splines = FALSE}
#' 
#' 
#' @details 
#' 
#' See \code{\link{sir}} for help on SIR/SMR estimation in general; usage of splines
#' is discussed below.
#' 
#' \strong{The spline variables}
#' 
#' The model can include one, two or three splines variables.
#' Variables can be included in the same model selecting \code{dependent.splines = TRUE}
#' and SIR ratios are calculated (first one is the SIR, others SIR ratios). 
#' Reference points vector can be set via \code{reference.points}
#' where first element of the vector is the reference point for first ratio.
#' 
#' Variable(s) to fit splines are given as a vector in argument \code{spline}.
#' Order will affect the results.
#' 
#' 
#' \strong{dependent.splines} 
#' 
#' By default dependent.splines is FALSE and all splines are fitted in separate models. 
#' If TRUE, the first variable in \code{spline} is a function of a SIR and other(s) are ratios.
#' 
#' \strong{knots}
#' 
#' There are three options to set knots to splines:
#' 
#' Set the number of knots for each spline variable with a \strong{vector}. 
#' The knots are automatically placed to the quantiles of observed cases in cohort data. 
#' The first and last knots are always the maximum and minimum values, so knot
#' value needs to be at least two.
#' 
#' Predefined knot places can be set with a \strong{list} of vectors.
#' The vector for each spline in the list specifies the knot places. The lowest 
#' and the largest values are the boundary knots and these should be checked beforehand.
#' 
#' If \code{knots} is left \strong{NULL}, the model searches the optimal number 
#' of knots by model AIC by fitting models iteratively from 2 to 15 knots and 
#' the one with smallest AIC is selected.
#' If \code{dependent.splines = TRUE}, the number of knots is searched by fitting each spline
#' variable separately.
#' 
#' 
#' \strong{print}
#' 
#' Splines can be stratified by the levels of variable given in \code{print}. If 
#' \code{print} is a vector, only the first variable is accounted for. The knots 
#' are placed globally for all levels of \code{print}. This also ensures that the likelihood 
#' ratio test is valid.
#' Splines are also fitted independently for each level of \code{print}.
#' This allows for searching interactions, e.g. by fitting spline for period 
#' (\code{splines='period'}) for each agegroup (\code{print = 'agegroup'}).
#' 
#' 
#' \strong{p-values}
#' 
#' The outputted p-value is a test of whether the splines are equal (homogenous)
#' at different levels of \code{print}. 
#' The test is based on the likelihood ratio test, where the full model 
#' includes \code{print} and is 
#' compared to a null model without it.
#' When \code{(dependent.splines = TRUE)} the p-value returned is a global p-value.
#' Otherwise the p-value is spline-specific.
#' 
#' 
#' @return A list of date.frames and vectors.
#' Three spline estimates are named as \code{spline.est.A/B/C} and the corresponding values
#' in \code{spline.seq.A/B/C} for manual plotting
#' 
#' 
#' @seealso \code{\link{plot.sirspline}}, \code{\link{sir}}, \code{\link{splitMulti}} 
#' \href{../doc/sir.html}{A SIR calculation vignette}
#' 
#' @export sirspline
#' @import data.table 
#' @import splines
#' @import stats
#' 
#' @examples \dontrun{
#' ## for examples see: vignette('sir')
#' }

sirspline <- function( coh.data, 
                       coh.obs,
                       coh.pyrs,
                       ref.data = NULL,
                       ref.obs = NULL,
                       ref.pyrs = NULL, 
                       ref.rate = NULL,
                       subset = NULL,
                       print = NULL,
                       adjust = NULL,
                       mstate = NULL,
                       spline,
                       knots = NULL,
                       reference.points = NULL,
                       dependent.splines = TRUE){
  
  coh.data <- data.table(coh.data)
  
  ## subsetting-----------------------------------------------------------------
  ## no copy taken of data!
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data = coh.data, substiset = subset)
  coh.data <- coh.data[subset,]
  
  # print list --------------------------------------------------------------

  env1 <- environment()
  coh.data <- data_list(data = coh.data, arg.list = substitute(print), env = env1)
  
  mstate <- as.character(substitute(mstate))
  if(length(mstate) == 0) {
    mstate <- NULL
  }
  if(!is.null(mstate)) {
    coh.data[,(mstate) := 0L] 
  }
  
  # evalPopArg
  
  spline <- substitute(spline)
  c.spl <- evalPopArg(data = coh.data, arg = spline)
  spline <- names(c.spl)
  
  coh.obs <- substitute(coh.obs)
  c.obs <- evalPopArg(data = coh.data, arg = coh.obs)
  coh.obs <- names(c.obs)
  
  coh.pyrs <- substitute(coh.pyrs)
  c.pyr <- evalPopArg(data = coh.data, arg = coh.pyrs)
  coh.pyrs <- names(c.pyr)
  
  print <- substitute(print)
  c.pri <- evalPopArg(data = coh.data, arg = print)
  print <- names(c.pri)
  
  adjust <- substitute(adjust)
  c.adj <- evalPopArg(data = coh.data, arg = adjust)
  adjust <- names(c.adj)

  # collect data
  coh.data <- cbind(c.obs, c.pyr, c.spl)
  if(!is.null(print))  {
    coh.data <- cbind(coh.data, c.pri[, print[!print %in% spline], with=FALSE])
  }
  if(!is.null(adjust)) {
    coh.data <- cbind(coh.data, c.adj[, adjust[!adjust %in% spline], with=FALSE])
  }
  
  if( !is.null(ref.data) ){
    ref.obs <- as.character(substitute(ref.obs))
    ref.pyrs <- as.character(substitute(ref.pyrs))
    ref.rate <- as.character(substitute(ref.rate))
    
    if (length(ref.obs) == 0) ref.obs <- NULL
    if (length(ref.pyrs) == 0) ref.pyrs <- NULL
    if (length(ref.rate) == 0) ref.rate <- NULL
  }
  
  st <- sir_table( coh.data = coh.data, 
                   coh.obs = coh.obs,
                   coh.pyrs = coh.pyrs,
                   ref.data = ref.data,
                   ref.obs = ref.obs,
                   ref.pyrs = ref.pyrs, ref.rate = ref.rate,
                   print = print,
                   adjust = adjust,
                   mstate = mstate,
                   spline = spline)

  results <- sir_spline( table = st, 
                         print = print,
                         adjust = adjust,
                         spline = spline,
                         knots = knots,
                         reference.points = reference.points,
                         dependent.splines = dependent.splines)
  
  setclass(results, c('sirspline', 'pe', class(results)))
  return(results)  
}






# Input: two data.table:s
# output: one data.table including rates
#' @import stats
#' @import data.table
sir_table <- function( coh.data, 
                       coh.obs,
                       coh.pyrs,
                       ref.data = NULL,
                       ref.obs = NULL,
                       ref.pyrs = NULL, 
                       ref.rate = NULL,
                       print = NULL,
                       adjust = NULL,
                       spline = NULL,
                       mstate = NULL) {


  # initial checks -------------------------------------------------
  
  if(is.null(ref.data)) {
    if(is.null(print)){
      stop('Both ref.data and print cannot be NULL.')
    }
    ref.data <- data.table(coh.data)
    ref.obs <- coh.obs
    ref.pyrs <- coh.pyrs
  }
  
  coh.data <- data.table(coh.data)
  ref.data <- data.table(ref.data)

  vl <- unique( c(coh.pyrs, coh.obs, adjust, print) )
  if( !is.null(mstate) )  {
    vl <- vl[which( vl != mstate )]
  }
  all_names_present(coh.data, vl )
  
  if ( !is.null(ref.pyrs) & !is.null(ref.obs) ) {
    all_names_present(ref.data, c(ref.pyrs, ref.obs, adjust))
  }
  
  # Melt lexpand data -------------------------------------------------------
  
  if( length(coh.obs) > 1 ) {
    if( is.null(mstate) ){
      stop('coh.obs length is > 1. Set variable name for mstate.')
    }
    if( !mstate %in% names(ref.data) ){
      warning('mstate variable name does not match names in ref.data.')
    }

    aggre <- unique(c(adjust, print, spline, coh.pyrs))
    aggre <- aggre[which(aggre != mstate)]

    coh.data <- melt( data = coh.data, id.vars = aggre, measure.vars = coh.obs, 
                      value.name = 'coh.observations', 
                      variable.name = mstate, variable.factor = FALSE)
    coh.obs <- 'coh.observations'
  
    # parse Y name form string 'formXtoY'
    q <- quote(
      robust_values(substr(get(mstate), 
                           start = regexpr( pattern = 'to', text = get(mstate) ) + 2, 
                           stop  = nchar(x = get(mstate) )))
      )
    coh.data[,(mstate) := eval(q) ]
    
    if( !(mstate %in% adjust)) {
      warning('Consider including mstate variable also in adjust. See help(sir) for details.')
    }
  }
  
  # prepare data steps, reduce dimensions -----------------------------------
  
  setnames(coh.data, c(coh.obs, coh.pyrs), c('coh.observations','coh.personyears'))  
  
  
  coh.data <- expr.by.cj(data = coh.data,
                         by.vars = unique( sort(c(adjust, print, spline)) ), 
                         expr = list(coh.observations = sum(coh.observations), 
                                     coh.personyears  = sum(coh.personyears)))
  coh.data <- na.omit(coh.data)
  # rates
  if( !is.null(ref.rate) ){
    setnames(ref.data, ref.rate, 'ref.rate')
    ref.data <- expr.by.cj(data = ref.data, by.vars = c(adjust), 
                           expr = list(ref.rate = mean(ref.rate)))
  } else {
    setnames(ref.data, c(ref.obs, ref.pyrs), c('ref.obs','ref.pyrs'))
    ref.data <- expr.by.cj(data = ref.data, by.vars = c(adjust), 
                           expr = list(ref.obs = sum(ref.obs),
                                       ref.pyrs= sum(ref.pyrs)))
    ref.data[, ref.rate := ref.obs / ref.pyrs ]
  }
  
  # Merge
  sir.table <- merge(coh.data, ref.data, by=c(adjust), all.x=TRUE)
  sir.table[, expected := ref.rate * coh.personyears]
  sir.table <- na2zero(sir.table)
  
  if ( !is.null(print) | !is.null(spline)){
    sir.table <- sir.table[ ,list(observed = sum(coh.observations), 
                                  expected = sum(expected),
                                  pyrs = sum(coh.personyears)), 
                           by = c(unique(c(print, spline)))]
    setkeyv(sir.table, c(print, spline))
  }
  else {
    sir.table <- sir.table[ ,list(observed = sum(coh.observations), 
                                  expected = sum(expected),
                                  pyrs = sum(coh.personyears))]
  }
  
  return(sir.table)
}




# Input: sir.table
# Output: list of data.tables and values
sir_est <- function( table,
                     print = NULL,
                     adjust = NULL,
                     test.type = 'homogeneity',
                     alpha = 0.95,
                     p.adj = NULL,
                     EAR = FALSE,
                     round.by = 2,
                     round.by.pvalue = 4) {
  pyrs <- NULL ## APPEASE R CMD CHECK
  
  # functions -----------------------------------------------------------
  
  # function to SIR p-value
  chi.p <- function(o, e) {
    pchisq( ( (abs(o - e) - 0.5)^2)/e, df=1, lower.tail=FALSE)
  }
  
  # Univariate and Total SIR ----------------------------------------------
  sir.table <- data.table(table)  
  sir.table[ ,':='(sir = observed / expected,
                   lower_ci = poisson.ci(observed, expected, conf.level=alpha)[,4],
                   upper_ci = poisson.ci(observed, expected, conf.level=alpha)[,5],
                   p_value  = chi.p(observed, expected)) ]
  # adjusted p-value
  if( !is.null(p.adj) ) {
    sir.table[ ,p_adj := p.adjust(p_value, method = p.adj)]
  }
  
  # total sir
  combined <- data.table(table)
  combined <- combined[,list(observed = sum(observed), 
                             expected = sum(expected),
                             pyrs = sum(pyrs))]
  combined[ ,':='(sir = observed/expected,
                  lower_ci = poisson.ci(observed, expected, conf.level=alpha)[,4],
                  upper_ci = poisson.ci(observed, expected, conf.level=alpha)[,5],
                  p_value  = chi.p(observed, expected))]
  
  # set proper column names to CI
  lower_name <- paste( (1 - alpha)/2*100, '%')
  upper_name <- paste( (1 - (1 - alpha)/2)*100, '%')
  setnames(sir.table, c('lower_ci','upper_ci'), c(lower_name, upper_name))
  setnames(combined, c('lower_ci','upper_ci'), c(lower_name, upper_name))
  
  
  # Poisson regression ------------------------------------------------------
  
  # write model
  fa <- a <- NULL
  
  #print <- c('fot','year')
  sir.table
  if(!is.null(print)){
    fa <- rev(print)
  }
    
  # model formula
  if (!is.null(fa)){
    a <- paste0('as.factor(',paste( fa, collapse = '):as.factor('),')')
    sir.formula <- paste('observed ~ 0 +', a)
  }
  else {
    sir.formula <- paste('observed ~ 1') 
  }

  # fit model if possible -----------------------------------------------------
  sir.multi <- data.table(table)
  multi <- NULL
  fit <- tryCatch(do.call("glm", list(formula = terms(as.formula(sir.formula),
                                                      keep.order = FALSE), 
                                      offset = log(sir.table[,expected]), 
                                      data = sir.table, family=poisson(log))), 
                  error=function(f) NULL )
  # LRT test (homogeneity) ----------------------------------------------------
  
  lrt_sig <- NULL  
  if( sir.formula != 'observed ~ 1' & !is.null(fit) ) {
    
    # fit null-model
    
    # homogeneity test
    if( test.type == 'homogeneity' ){    
      fit_full <- tryCatch(
        do.call("glm", list(formula = terms(as.formula( paste0('observed ~ 1 + ', a) )), 
                            offset = log(sir.table[,expected]), 
                            data = sir.table, family=poisson(log))), 
        error=function(f) NULL )
    }
    
    # trend test
    else if( test.type == 'trend' ){
      fit_full <- tryCatch(
        do.call("glm", list(formula = terms(as.formula( paste0('observed ~ 1 + ', paste(print, collapse=' + ') ) )), 
                            offset = log(sir.table[,expected]), 
                            data = sir.table, family=poisson(log))), 
        error=function(f) NULL )
    } 
    else {
      stop('Select test.type: "homogeneity" or "trend".')
    }        
    fit_null <- tryCatch(
      do.call("glm", list(formula = terms(as.formula('observed ~ 1') ), 
                          offset = log(sir.table[,expected]), 
                          data = sir.table, family=poisson(log))), 
      error=function(f) NULL )
    if(!is.null(fit_full) ) {
      lrt <- anova(fit_full, fit_null, test = 'Chisq')
      lrt_sig <- lrt[['Pr(>Chi)']][2]    
    }
  }
  
  # confidence intervals
  
  if (!is.null(fit)) {
    ci <- NULL
    ci <- suppressMessages( suppressWarnings( 
      tryCatch(exp(confint(fit, level=alpha)), error=function(e) NULL )
    ) )
    ci.info <- 'Confidence intervals calculated from profile-likelihood.'
    wald <- FALSE
    if(is.null(ci)) {
      ci <- cbind(lower_name=as.vector(exp(fit[[1]] - sqrt(diag(vcov(fit)))*qnorm((1-alpha)/2, mean = 0, sd = 1, lower.tail = FALSE))), 
                  upper_name=as.vector(exp(fit[[1]] + sqrt(diag(vcov(fit)))*qnorm((1-alpha)/2, mean = 0, sd = 1, lower.tail = FALSE))))
      ci.info <- 'Confidence intervals are normal-approximated (Wald)'
      wald <- TRUE
    }
    message(ci.info)
    
    # collect results
    
    if ( !is.null(print) ) {
      multi <- tryCatch( data.table(sir.multi, 
                                    sir = na.omit(exp(as.numeric(coef(fit)))),
                                    data.table(na.omit(ci)),
                                    p_value = as.vector(summary(fit)$coef[, "Pr(>|z|)"]) ),
                         error = function(e) NULL )
    } 
    else {
      multi <- tryCatch( data.table(sir.multi, 
                                    sir = na.omit(exp(as.numeric(coef(fit)))),
                                    lower_name = ci[1], upper_name = ci[2], 
                                    p_value = as.vector(summary(fit)$coef[, "Pr(>|z|)"]) ),
                         error = function(e) NULL )
      setnames(multi, c('lower_name','upper_name'), c(lower_name, upper_name))
      wald <- FALSE 
    }
    if ( wald ){
      setnames(multi, c('lower_name','upper_name'), c(lower_name, upper_name))
    }
  } 
  # model fit failed
  else {
    warning('Could not fit glm')
    Model <- NULL
  }
  
  # Round results -----------------------------------------------------------
  
  round2 <- function(x, by=round.by){
    round(x, by)
  }
  round3 <- function(x, by=round.by.pvalue){
    round(x, by)
  }
  round_cols <- function(data){
    if(!is.null(data)) {
      cols <- names(data)
      # drop p-values
      adjust.cols <- which(cols == 'observed'):length(cols)
      cols <- cols[adjust.cols]
      cols1 <- cols[ which(substr(cols, 1, 2) != 'p_')]
      cols2 <- cols[ which(substr(cols, 1, 2) == 'p_')]
      data[,(cols1) := lapply(.SD, round2), .SDcols=cols1]
      data[,(cols2) := lapply(.SD, round3), .SDcols=cols2]
    }
  }
  
  multi <- round_cols(multi)
  sir.table <- round_cols(sir.table)
  combined <- round_cols(combined)
  
  # check that univar and model SIR match -----------------------------------
  
  if(!is.null(multi)) {
    if( !identical( sir.table[,floor(sir)],  multi[,floor(sir)]) ) {
      warning("SIR's in model and univariate output doesn't match")
    }
  }
  
  # EAR
  if (EAR) {
    sir.table[,EAR := round((observed - expected)/pyrs * 1000, 2)]
    if(!is.null(multi)) {
      multi[,EAR := round((observed - expected)/pyrs * 1000, 2)]
    }
  }
  
  
  # Print results -----------------------------------------------------------
  results <- list(total = combined, 
                  univariate = sir.table, 
                  model = multi,
                  adjusted = adjust,
                  lrt.test = lrt_sig,
                  test.type = test.type)
  return(results)
}

# Input: sir.table
# Output: estimates and sequences for plotting splines
#' @import splines
#' @import data.table
#' @import stats
sir_spline <- function(  table,
                         print = NULL,
                         adjust = NULL,
                         spline,
                         knots = NULL,
                         reference.points = NULL,
                         dependent.splines = TRUE){
  knts <- 
    spline.seq.A <- 
    spline.seq.B <- 
    spline.seq.C <- 
    spline.est.A <- 
    spline.est.B <- 
    spline.est.C <- NULL
  
  if (!is.null(knots) & length(knots) != length(spline) ) {
    stop('Arguments spline and knots has to be same length.')
  }
  
  
  # Spline functions -------------------------------------------------------
  
  # function to get spline seq
  spline.seq <- function(data, spline.var=NULL) {
    # palauttaa jotaina
    if(is.na(spline.var)) {
      return(NULL)
    }
    spline.seq <- seq( min( data[,get(spline.var)] ), 
                       max( data[,get(spline.var)] ), length.out = 100)
    return(spline.seq)
  }
  
  # function to search optimal number of knots by AIC
  spline.knots <- function(data, knots = NULL, spline.vars = NULL){    
    # search optimal number of knots
    if( is.null(knots) ) {
      knts <- list()
      for (jj in 1:length(spline.vars)) {
        # reduce data to fit model
        data0 <- data[,list(observed=sum(observed), expected = sum(expected)), by = eval(spline.vars[jj])] 
        data0 <- data0[expected > 0]
        spline.fit <- glm(observed ~ 1, offset=log(expected), family=poisson(log), data = data0)
        aic0 <- summary(spline.fit)[['aic']]
        limit <- 20
        ii <- 2
        while(  ii < limit ){
          tmp.knots <- ii
          knts[jj] <- list( data0[ ,quantile( rep(get(spline.vars[jj]),observed), probs = seq(0,100,length.out = tmp.knots)/100)] )
          spline.fit <- glm(observed ~ Ns(get(spline.vars[jj]), knots = knts[[jj]]), offset=log(expected), family=poisson(log), data=data0)
          aic0 <- c(aic0, summary(spline.fit)[['aic']])
          ii <- ii + 1
        }
        tmp.knots <- which(aic0 == min(aic0))[1]
        if(tmp.knots == 1) {
          message(paste0('Null model better than spline in ', jj))
          tmp.knots <- 2
        }
        knts[jj] <- list(data0[ ,quantile( rep(get(spline.vars[jj]),observed), probs = seq(0,100,length.out = tmp.knots)/100)])
        rm(tmp.knots)
      }
      knots <- unlist(lapply(knts, length))
    }
    else {
      # knot predefined
      if( is.list(knots) ){
        knts <- knots
        knots <- unlist(lapply(knots, length))
      } 
      # knot number predefined
      else {
        if( any(knots < 2) ) { 
          message('Min knots number set to 2.') 
          knots[knots < 2] <- 2
        }
        knts <- list()
        for(i in 1:length(knots)) {
          knts[i] <- list( data[ ,quantile( rep(get(spline.vars[i]), observed), probs = seq(0,100,length.out = knots[i])/100)])
        }
      }
    }
    names(knts) <- spline.vars
    return(knts)
  }
  
  # function to estimate 2-3 dim splines in same model
  spline.estimates.dep <- function(sir.spline = sir.spline,
                                   spline.seq.A = spline.seq.A,
                                   spline.seq.B = spline.seq.B,
                                   spline.seq.C = spline.seq.C,
                                   reference.points = reference.points,
                                   knts = knts
  ){
    
    if( all(!is.null(reference.points), (length(reference.points) + 1) != length(spline)) ){
      stop('Parameter reference.points length should be length of spline - 1.')
    }
    
    
    form <- 'Ns(get(spline[[1]]), kn=knts[[1]])'
    nsA <- Ns( spline.seq.A, knots = knts[[1]])
    if ( length(spline) >= 2) {
      form <- paste0(form, ' + Ns(get(spline[[2]]), kn=knts[[2]])')
      nsB <- Ns( spline.seq.B, knots = knts[[2]])
    }
    if ( length(spline) == 3) {
      form <- paste0(form, ' + Ns(get(spline[[3]]), kn=knts[[3]])')
      nsC <- Ns( spline.seq.C, knots = knts[[3]])
    }
    
    form <- paste0('observed ~ ', form)
    spline.fit <- do.call("glm", list(formula = as.formula(form),
                                      offset = log(sir.spline[expected > 0,expected]),
                                      family = poisson,
                                      data = sir.spline[expected>0]))
    if( any( ci.exp(spline.fit)[,1] == 1) ){
      message("NA's in spline estimates.")
    }
    
    aic <- summary(spline.fit)[['aic']]
    
    rf.C <- rf.B <- NA
    # set assigned reference points or get minimum values
    if( !is.null(reference.points) ) {
      rf.B <- reference.points[1]
      rf.C <- reference.points[2]
    } 
    else {
      rf.B <- min( sir.spline[,get(spline[2])] )
      if(!is.na(spline[3])) {
        rf.C <- min( sir.spline[,get(spline[3])] )
      }
    }
    
    if( !is.na(rf.B) )  {
      B <- Ns( rep(rf.B, 100), knots = knts[[2]])
      if( findInterval(rf.B, range(sir.spline[,get(spline[2])])) != 1 ) {
        message("WARNING: reference point 2 doesn't fall into spline variable interval")
      }
    }
    
    if( !is.na(rf.C) ){
      C <- Ns( rep(rf.C, 100), knots = knts[[3]])
      if( findInterval(rf.C, range(sir.spline[,get(spline[3])])) != 1) {
        message("WARNING: reference point 3 doesn't fall into spline variable interval")
      }
    } 
    
    # make subset of model parameters
    if( !is.null(knts[2]) ) {
      sub.B <- which( grepl('spline[[2]]', names(spline.fit$coefficients),fixed = TRUE) )
    }
    if( !is.null(knts[3]) ) {
      sub.C <- which( grepl('spline[[3]]', names(spline.fit$coefficients),fixed = TRUE) )
    }
    if ( length(spline) == 2) {
      spline.est.A <- ci.exp(spline.fit, ctr.mat = cbind(1, nsA, nsB))
      spline.est.B <- ci.exp(spline.fit, subset = sub.B, ctr.mat = nsB - B)
      spline.est.C <- NULL      
    }
    if ( length(spline) == 3) {
      spline.est.A <- ci.exp(spline.fit, ctr.mat = cbind(1, nsA, nsB, nsC))
      spline.est.B <- ci.exp(spline.fit, subset= sub.B, ctr.mat = nsB - B)
      spline.est.C <- ci.exp(spline.fit, subset= sub.C, ctr.mat = nsC - C)
    }
    list(a = spline.est.A, 
         b = spline.est.B, 
         c = spline.est.C)
  }
  
  # function to estimate independet splines
  spline.estimates.uni <- function(data, spline.var, spline.seq, knots, knum) {  
    if(is.na(spline.var)) return(NULL)
    knots <- knots[[knum]]
    data <- data[,list(observed=sum(observed), expected = sum(expected)), by = eval(spline.var)][expected > 0]
    spline.uni <- glm(observed ~ Ns(get(spline.var), knots = knots), offset=log(expected), family=poisson(log), data = data)
    nsx <- Ns( spline.seq, knots = knots)
    spline.est <- ci.exp(spline.uni, ctr.mat = cbind(1, nsx))
    spline.est
  }
  
  
  
  # Poisson regression Splines -------------------------------------------------
  
  sir.spline <- data.table(table)
  
  # convert spline variables to numeric
  temp.fun <- function(x){
    as.numeric(as.character(x))
  }
  sir.spline[, (spline) := lapply(.SD, temp.fun), .SDcols = spline]
  
  
  
  # set knots
  knts <- spline.knots(data=sir.spline, knots = knots, spline.vars = spline)  
  
  # set sequences
  spline.seq.A <- spline.seq(data=sir.spline, spline.var=spline[1])
  spline.seq.B <- spline.seq(data=sir.spline, spline.var=spline[2])
  spline.seq.C <- spline.seq(data=sir.spline, spline.var=spline[3])
  
  if( length(spline) == 1 ) {
    dependent.splines <- FALSE
  }
  
  # convert print to factor
  print <- print[1]
  
  # loop for each level of print:
  if( !is.null(print) ) {
    prnt.levels <- sir.spline[,unique( get(print) )]
    sir.spline[,(print) := factor(get(print))]
  }
  else {
    print <- 'temp'
    sir.spline[,temp := 1]
    prnt.levels <- 1
  }
  
  spline.est.A <- NULL
  spline.est.B <- NULL
  spline.est.C <- NULL
  
  for(i in prnt.levels){
    if( dependent.splines ) {
      out <- spline.estimates.dep(sir.spline = sir.spline[get(print) == i],
                                  spline.seq.A = spline.seq.A,
                                  spline.seq.B = spline.seq.B,
                                  spline.seq.C = spline.seq.C,
                                  reference.points = reference.points,
                                  knts = knts)
      est.A <- out[['a']]
      est.B <- out[['b']]
      est.C <- out[['c']]
    }
    else{
      est.A <- spline.estimates.uni(data = sir.spline[get(print) == i], spline.var = spline[1], spline.seq = spline.seq.A, knots = knts, knum = 1)
      est.B <- spline.estimates.uni(data = sir.spline[get(print) == i], spline.var = spline[2], spline.seq = spline.seq.B, knots = knts, knum = 2)  
      est.C <- spline.estimates.uni(data = sir.spline[get(print) == i], spline.var = spline[3], spline.seq = spline.seq.C, knots = knts, knum = 3)
    }
    
    add_i <- function(est.x, i){
      if(is.null(est.x)) {
        return(NULL)
      }
      cbind(i, data.frame(est.x))
    }
    
    
    est.A <- add_i(est.A, i)
    est.B <- add_i(est.B, i)
    est.C <- add_i(est.C, i)
    
    spline.est.A <- rbind(spline.est.A, est.A)
    spline.est.B <- rbind(spline.est.B, est.B)
    spline.est.C <- rbind(spline.est.C, est.C)
  }
  
  # get p-value and anova-table
  anovas <- NULL
  p <- NULL
  if(dependent.splines) {
    form.a <- 'Ns(get(spline[[1]]), kn=knts[[1]]) + Ns(get(spline[[2]]), kn=knts[[2]])'
    form.b <- 'get(print):Ns(get(spline[[1]]), kn=knts[[1]]) + get(print):Ns(get(spline[[2]]), kn=knts[[2]])'
    if ( length(spline) == 3) {
      form.a <- paste0(form.a, ' + Ns(get(spline[[3]]), kn=knts[[3]])')
      form.b <- paste0(form.b, ' + get(print):Ns(get(spline[[3]]), kn=knts[[3]])')
    }
    
    fit.fun <- function( form.string ){
      do.call("glm", list(formula = as.formula( form.string ),
                          offset = log(sir.spline[expected > 0,expected]),
                          family = poisson,
                          data = sir.spline[expected>0]))
    }
    
    fit.1 <- fit.fun( paste0('observed ~ ', form.a) )
    fit.2 <- fit.fun( paste0('observed ~ ', 'get(print)+', form.a))
    fit.3 <- fit.fun( paste0('observed ~ ', form.b))
    fit.4 <- fit.fun( paste0('observed ~ ', 'get(print)+', form.b) )
    
    global.p<- anova(fit.4, fit.1, test='LRT')
    level.p <- anova(fit.2, fit.1, test='LRT')
    #shape.p <- anova(fit.4, fit.3, test='LRT')
    
    anovas <- list(global.p = global.p, level.p = level.p)
    p <- rbind(global.p[['Pr(>Chi)']][2], level.p[['Pr(>Chi)']][2]) # , shape.p, 
  }
  else {    
    lrt.uni <- function(data=sir.spline, spline.var=spline[1], print=print, knots=knts, knum = 1) {
      if (is.na(spline.var)) return (NULL)
      data <- data.table(data)
      knots <- knots[[knum]]
      fit0 <- glm(observed ~ get(print)+Ns(get(spline.var), knots = knots), offset=log(expected), family=poisson(log), data = data[expected>0])
      fit1 <- glm(observed ~ Ns(get(spline.var), knots = knots), offset=log(expected), family=poisson(log), data = data[expected>0])
      fit2 <- glm(observed ~ get(print)*Ns(get(spline.var), knots = knots), offset=log(expected), family=poisson(log), data = data[expected>0])
      anova(fit2,fit1,fit0, test='Chisq') # [['Pr(>Chi)']][2]
    }
    
    var1.p <- lrt.uni(spline.var = spline[1], print=print, knots=knts, knum = 1)
    var2.p <- lrt.uni(spline.var = spline[2], print=print, knots=knts, knum = 2)
    var3.p <- lrt.uni(spline.var = spline[3], print=print, knots=knts, knum = 3)
    
    p <- list(spline.a = var1.p[['Pr(>Chi)']][2], 
              spline.b = var2.p[['Pr(>Chi)']][2], 
              spline.c = var3.p[['Pr(>Chi)']][2])
    anovas <- list(spline.a = var1.p, spline.b = var2.p, spline.c = var3.p)
  }
  
  output <- list( spline.est.A = spline.est.A,
                  spline.est.B = spline.est.B,
                  spline.est.C = spline.est.C,
                  spline.seq.A = spline.seq.A,
                  spline.seq.B = spline.seq.B,
                  spline.seq.C = spline.seq.C,
                  adjust = adjust,
                  print = print,
                  spline = spline,
                  anovas = anovas,
                  knots = knts,
                  spline.dependent = dependent.splines,
                  p.values = p)
  output
}

# input data and argument list. replaces print in upper environment with name a vector.
data_list <- function( data, arg.list, env ) {
  if(missing(env)){
    arg.list <- substitute(arg.list)
    env <- parent.frame()
  }
  d <- data.table(data)
  
  l <- eval(arg.list, envir = d, enclos = parent.frame()) 
  
  if( is.list( l ) ) {
    n <- intersect(names(l), names(d))
    if(length(n)>0){
      d[,(n) := NULL]
    }
    #     if(is.null(names(l))) {
    #       v <- 1:length(l)
    #       setnames(l, v, paste0('V', v))
    #     }
    l <- as.data.table(l)
    l <- data.table(l)
    assign('print', colnames(l), envir = env) # set names to parent environment
    if( ncol(d) > 0) {
      l <- data.table(d, l)
    }
    return(l)
  } else { 
    return(data) 
  }
}


globalVariables(c('observed','expected','p_adj','p_value','temp','coh.observations','coh.personyears'))

