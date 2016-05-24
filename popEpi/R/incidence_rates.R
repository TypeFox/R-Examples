#' @title Direct-Standardised Incidence/Mortality Rates
#' @author Matti Rantanen, Joonas Miettinen
#'
#' @description \code{rate} calculates adjusted rates using
#' preloaded weights data or user specified weights.
#'
#' @param data aggregated data (see e.g. \code{\link{lexpand}}, 
#' \code{\link{aggre}} if you have subject-level data)
#' @param pyrs person-years variable name in data.
#' \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{pyrs = pyrs}.
#' @param obs observations variable name in data.
#' \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{obs = obs}.
#' @param adjust variable for adjusting the rates.
#' \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{adjust = agegroup}.
#' @param print variable name to stratify the rates.
#' \link[=flexible_argument]{Flexible input}, typically e.g.
#' \code{print = sex} or \code{print = list(sex, area)}.
#' @param weights typically a list of weights or a \code{character} string
#' specifying an age group standardization scheme; see
#' the \link[=direct_standardization]{dedicated help page} 
#' and examples.
#' 
#' @param subset a logical expression to subset data.
#' 
#' @details Input data needs to be in aggregated format with observations 
#' and person-years. For individual data use \code{\link{lexpand}}, or
#' \code{\link{ltable}} and merge person-years manually.
#' 
#' 
#' @return Returns a \code{data.table} with observations, person-years, rates and
#' adjusted rates, if availble. Results are stratified by \code{print}.
#' Adjusted rates are identified with suffix \code{.adj} and  
#' \code{.lo} and \code{.hi} are for confidence intervals lower and upper 
#' 95\% bounds, respectively.
#' The prefix \code{SE.} stands for standard error.
#' 
#' @seealso \code{\link{lexpand}}, \code{\link{ltable}}
#' 
#' @examples 
#' ## Prepare data with lexpand and then reformat agegroup.
#' data(sibr)
#' x <- lexpand(sibr, birth = bi_date, entry = dg_date, exit = ex_date,  
#'              breaks = list(per = c(1990,2000,2010,2020), age = c(0:17*5,Inf)),
#'              aggre = list(agegroup = age, year.cat = per),
#'              status =  status != 0)
#'
#' x$agegroup <- cut(x$agegroup,  c(0:17*5,Inf), right = FALSE)
#'
#' ## calculate rates for selected periods with Nordic 2000 weights:
#' r1 <- rate( data = x, obs = from0to1, pyrs = pyrs, print = year.cat, 
#'             adjust = agegroup, weights = 'nordic')
#' r1
#'
#' ## use total person-years by stratum as weights (some have zero)
#' w <- ltable(x, by.vars = "agegroup", expr = sum(pyrs))
#' w[is.na(w$V1),]$V1 <- 0
#' 
#' r2 <- rate( data = x, obs = from0to1, pyrs = pyrs, print = year.cat, 
#'             adjust = agegroup,
#'             weights = w$V1)
#' r2
#' 
#' ## use data.frame of weights:
#' names(w) <- c("agegroup", "weights")
#' r2 <- rate( data = x, obs = from0to1, pyrs = pyrs, print = year.cat, 
#'             adjust = agegroup,
#'             weights = w)
#' r2
#' 
#' ## internal weights (same result as above)
#' r3 <- rate( data = x, obs = from0to1, pyrs = pyrs, print = year.cat, 
#'             adjust = agegroup,
#'             weights = "internal")
#' r3
#'
#' @import data.table
#' @export rate

rate <- function( data,
                  obs = NULL,
                  pyrs = NULL,
                  print = NULL,
                  adjust = NULL, 
                  weights = NULL,
                  subset = NULL
) {
  
  PF <- parent.frame(1L)
  TF <- environment()
  
  ## subsetting -----------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data = data, substiset = subset, enclos = PF)
  data <- data[subset,]
  setDT(data)
  
  # evalPopArg
  obs <- substitute(obs)
  inc.obs <- evalPopArg(data = data, arg = obs, enclos = PF)
  if (!length(inc.obs)) {
    stop("No observations given.")
  }
  obsNames <- copy(names(inc.obs))
  tmpObsNames <- makeTempVarName(data = data, pre = "obs")
  setnames(inc.obs, obsNames, tmpObsNames)
  
  pyrs <- substitute(pyrs)
  inc.pyr <- evalPopArg(data = data, arg = pyrs, enclos = PF)
  if (!length(inc.pyr)) {
    stop("No pyrs given.")
  }
  pyrNames <- copy(names(inc.pyr))
  tmpPyrNames <- makeTempVarName(data = data, pre = "pyr")
  setnames(inc.pyr, pyrNames, tmpPyrNames)
  
  print <- substitute(print)
  inc.pri <- evalPopArg(data = data, arg = print, enclos = PF)
  prNames <- tmpPrNames <- NULL
  if (length(inc.pri)) {
    prNames <- copy(names(inc.pri))
    tmpPrNames <- makeTempVarName(data = data, 
                                  pre = paste0("print", seq_along(prNames)))
    setnames(inc.pri, prNames, tmpPrNames)
  }
  
  adjust <- substitute(adjust)
  inc.adj <- evalPopArg(data = data, arg = adjust, enclos = PF)
  adNames <- tmpAdNames <- NULL
  if (length(inc.adj)) {
    adNames <- copy(names(inc.adj))
    tmpAdNames <- makeTempVarName(data = data, 
                                  pre = paste0("adjust", seq_along(adNames)))
    setnames(inc.adj, adNames, tmpAdNames)
  }
  
  ## collect data --------------------------------------------------------------
  data <- cbind(inc.obs, inc.pyr)
  if (!is.null(prNames)) data <- cbind(data, inc.pri) 
  if (!is.null(adNames)) data <- cbind(data, inc.adj)
  
  
  ## handle weights ------------------------------------------------------------
  weights <- substitute(weights)
  weights <- eval(weights, envir = PF)
  weights <- copy(weights)
  if (length(inc.adj)) {
    ## rename adjust variables in inc.adj back to original names
    ## for more human-readable errors in checkWeights if any occur
    setnames(inc.adj, tmpAdNames, adNames)
  }
  
  checkWeights(weights, inc.adj)
  if (is.list(weights) && !is.data.frame(weights)) {
    ## ensure weights list / DF names match to temp adjust var names
    weights <- weights[adNames]
    names(weights) <- tmpAdNames
  } else if (is.data.frame(weights)) {
    setnames(weights, adNames, tmpAdNames)
  }
  
  ## form table with weights ---------------------------------------------------
  NA.msg <- "Data contains %%NA_COUNT%% NA values."
  data <- makeWeightsDT(data, 
                        values = list(tmpObsNames, tmpPyrNames),
                        print = tmpPrNames, 
                        adjust = tmpAdNames, 
                        weights = weights, 
                        internal.weights.values = tmpPyrNames,
                        NA.text = NA.msg)
  
  ## estimate standardized rates -----------------------------------------------
  data <- rate_est(data = data,
                   obs = tmpObsNames,
                   pyrs = tmpPyrNames,
                   print = tmpPrNames,
                   weights = "weights")
  
  ## final touch ---------------------------------------------------------------
  setDT(data)
  setattr(data, "class", c("rate", "data.table", "data.frame"))
  setattr(data, name = 'rate.meta', value = list(obs = obsNames, 
                                                 pyrs = pyrNames, 
                                                 weights = weights,
                                                 adjust = adNames,
                                                 print = prNames,
                                                 call = match.call(),
                                                 NAs = NA))
  setnames(data, c(tmpObsNames, tmpPyrNames, tmpPrNames),
           c(obsNames, pyrNames, prNames))
  
  # data.frame output option  
  if (!getOption("popEpi.datatable")) {
    setDFpe(data)
  }
  
  return(data[])
}



stdr.weights <- function(wp = 'world00_1') {
  
  ## This one returns the standard population
  ## output: data.table with colnames: agegroup, reference
  ## standard populations are from datasets: stdpop18 and stdpop101
  allow.pop <- c("world_1966_18of5", 
                 "europe_1976_18of5", 
                 "nordic_2000_18of5", 
                 "world_2000_18of5", 
                 "world_2000_20of5", 
                 "world_2000_101of1")
  wp <- match.arg(wp, allow.pop)
  
  if (length(wp) > 1) {
    stop('Standard population name is not a scalar (length != 1).')
    
  } else if (wp %in% allow.pop[1:3]) {
    
    # get standard pop
    sr <- data.table(popEpi::stdpop18)
    setnames(sr, 1:4, c("agegroup",allow.pop[1:3]))
    sr[, agegroup := 1:18]
    sr[, setdiff(allow.pop[1:3], wp) := NULL]
    
    setnames(sr, wp, 'reference')
    
  } else if (wp %in% allow.pop[4:6]) {
    
    sr <- data.table(popEpi::stdpop101)
    if (wp == "world_2000_18of5") {
      sr[,agegroup := cut(agegroup, breaks=c(0:17*5,Inf), right=FALSE, labels=FALSE)]
      sr <- sr[,list(world_std = sum(world_std)), by="agegroup"]
    }
    if (wp == 'world_2000_20of5') {
      sr[,agegroup := cut(agegroup, breaks=c(0:19*5,Inf), right=FALSE, labels=FALSE)]
      sr <- sr[,list(world_std = sum(world_std)), by="agegroup"]
    }
    else {
      sr <- sr[,list(world_std = sum(world_std)), by="agegroup"]
    }
    setnames(sr, "world_std", "reference")
  }
  else {
    stop("Invalid standard population name.")
  }
  sr[]
}
globalVariables(c('stdpop18','stdpop101','agegroup','world_std'))


# rate_table <- function(data, 
#                        obs = 'obs',
#                        pyrs = 'pyrs',
#                        adjust = NULL,
#                        print = NULL,
#                        weight.data = 'world66_5',
#                        weights = NULL
# ) {
#   ## This one fetches and merges the standard population
#   ## or transforms the population already in the data to standard format.
#   
#   colsum1 <- function(c) c/sum(c)
#   # merge WHO weights to data
#   # Everything should sum to one on each level of print
#   data <- data.table(data)
#   if (!is.null(weights) && all(weights %in% colnames(data)) ) {
#     ## use predefined weights
#     if ( !is.null(print)) {
#       data[, reference := colsum1(.SD), .SDcols = weights, by = c(print)]
#     }
#     else {
#       data[, reference := colsum1(.SD), .SDcols = weights]
#     }
#     data[, (weights) := NULL]
#   }
#   else if ( !is.null(adjust) ) { # add: if ( !is.null(adjust) )
#     ## aggregate data before adding weights
#     eval0 <- paste0('list(obs = sum(',obs,',na.rm=TRUE),pyrs = sum(',pyrs,', na.rm=TRUE))')
#     eval0 <- parse(text = eval0)
#     data <- data[, eval(eval0), by = c(adjust, print) ] 
#     setnames(data, c('obs','pyrs'), c(obs, pyrs))
#     stdr.list <- c('world66_5','europe','nordic',"world00_1",
#                    "world00_20of5","world00_18of5")
#     
#     if ( !is.null(weight.data) && weight.data %in% stdr.list) {
#       ## get preloaded WHO std data
#       
#       if (length(adjust) > 1) stop('Set only one variable name for indicating age group')
# 
#       wd <- stdr.weights(wp = weight.data)
#       wd <- wd[, reference := colsum1(.SD), .SDcols = 'agegroup']
#       setnames(wd, 'agegroup', adjust)
#     }
#   
#     else if (is.null(weight.data) || weight.data=='cohort') {
#       # get cohort std
#       p1 <- paste0('sum(',pyrs,', na.rm=TRUE)')
#       p1 <- parse(text = p1)
#       wd <- data[,list( pyrs = eval(p1)), by = c(unique(adjust))]
#       wd[, reference := colsum1(.SD), .SDcols = 'pyrs']
#       wd[,c('pyrs') := NULL]
#     }
# 
#     data <- merge(x  = data, y  = wd[, c('reference', adjust) , with=FALSE],
#                   by = adjust, all.x = TRUE)
#   }
#   else {
#     if (!is.null(print)) {
#       data <- data[, list(obs = sum(get(obs), na.rm=TRUE), pyrs = sum(get(pyrs), na.rm=TRUE)), by = c(print) ]
#       setnames(data, c('obs','pyrs'), c(obs, pyrs))
#     }
#   }
#   return(data[])
# }
# 
# globalVariables(('reference'))

rate_est <- function(data = data, 
                     obs = 'obs', 
                     pyrs = 'pyrs', 
                     print = NULL, 
                     weights = NULL
) {
  ## This one estimates the rates and calculates CI's and SE's.
  
  badVars <- paste0("Internal error: missing following variable names in ",
                    "working data: %%VARS%%. Complain to the pkg maintainer ",
                    "if you see this.")
  all_names_present(data, c(obs, pyrs, print), msg = badVars)
  
  data <- data.table(data)
  if ( is.null(weights) |  !weights %in% colnames(data)) {
    weights <- NULL
  }
  
  if (all(!is.null(weights), !is.null(obs), !is.null(pyrs))) {
    # rate.adj
    
    f2 <- function(list) list[[1]]/list[[2]]*list[[3]]
    funx <- function(n,d,w,fun)  eval(parse(text=fun))
    

    # variance rate.adj for each strata A
    fun1 <- '(._d_/._n_^2) * ._w_^2'
    fun2 <- '._d_ / ._n_ * ._w_'
    
    make_fun <- function(n = NA, d = NA, w = NA, fun) {
      fun <- gsub(pattern = "._n_", replacement = n, x = fun)
      fun <- gsub(pattern = "._d_", replacement = d, x = fun)
      fun <- gsub(pattern = "._w_", replacement = w, x = fun)
      parse(text = fun)
    }
    eval.me1 <- make_fun(d = obs, n = pyrs, w=weights, fun = fun1)
    eval.me2 <- make_fun(d = obs, n = pyrs, w=weights, fun = fun2)
    data[, var.temp := eval(eval.me1)]
    data[, lam.temp := eval(eval.me2)]
    # add std weighted rates and variances
    #data[, ':='(var.temp = funx(d=get(obs), n=get(pyrs), w=get(weights), fun = fun1),
    #            lam.temp = funx(d=get(obs), n=get(pyrs), w=get(weights), fun = fun2)) ]
    data[, rate.adj := f2(.SD), .SDcols= c(obs, pyrs, weights)]

    # aggregate data
    ie <- paste0('list(', obs, '=sum(',obs,',na.rm=TRUE), ', pyrs, '=sum(',pyrs,',na.rm=TRUE),',
                 'rate.adj=sum(rate.adj,na.rm=TRUE),' ,'lam.temp=sum(lam.temp,na.rm=TRUE), var.temp=sum(var.temp,na.rm=TRUE))') 
    l <- parse(text = ie)

    data <- data[, eval(l), by=print]
    # rate.adj: S.E.
    data[, SE.rate.adj := exp( sqrt((1/lam.temp)^2 * var.temp)) ]
    # rate.adj: CI
    data[, ':='(rate.adj.lo = exp(log(rate.adj)-log(SE.rate.adj)*1.96),
                rate.adj.hi = exp(log(rate.adj)+log(SE.rate.adj)*1.96)) ]
    data[,c('lam.temp','var.temp') := NULL]
  }
  
  else {
    ie <- paste0('list(', obs, '=sum(',obs,'), ', pyrs, '=sum(',pyrs,'))') 
    l <- parse(text = ie)
    data <- data[, eval(l), by=print]
  }
  ia <- paste0('rate := ',obs,'/', pyrs)
  k <- parse(text = ia)
  data[, eval(k), by = print]
  eval.me3 <- paste('exp(1/',obs,')')
  eval.me3 <- parse(text = eval.me3)
  data[, SE.rate := eval(eval.me3)]
  data[, ':='(rate.lo = exp(log(rate)-log(SE.rate)*1.96),
              rate.hi = exp(log(rate)+log(SE.rate)*1.96)) ]
  return(data[])
}

globalVariables(c('var.temp','lam.temp','rate.adj','SE.rate.adj','SE.rate'))
