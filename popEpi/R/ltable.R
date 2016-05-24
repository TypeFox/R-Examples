#' @title Tabulate Counts and Other Functions by Multiple Variables into a 
#' Long-Format Table
#' @author Joonas Miettinen, Matti Rantanen
#' @description \code{ltable} makes use of \code{data.table} 
#' capabilities to tabulate frequencies or 
#' arbitrary functions of given variables into a long format 
#' \code{data.table}/\code{data.frame}. \code{expr.by.cj} is the 
#' equivalent for more advanced users.
#' @param data a \code{data.table}/\code{data.frame}
#' @param by.vars names of variables that are used for categorization, 
#' as a character vector, e.g. \code{c('sex','agegroup')}
#' @param expr object or a list of objects where each object is a function 
#' of a variable (see: details)
#' @param subset a logical condition; data is limited accordingly before
#' evaluating \code{expr} - but the result of \code{expr} is also
#' returned as \code{NA} for levels not existing in the subset. See Examples.
#' @param use.levels logical; if \code{TRUE}, uses factor levels of given 
#' variables if present;  if you want e.g. counts for levels
#' that actually have zero observatios but are levels in a factor variable, 
#' use this
#' @param na.rm logical; if \code{TRUE}, drops rows in table that have 
#' \code{NA} as values in any of \code{by.vars} columns
#' @param robust logical; if \code{TRUE}, runs the outputted data's 
#' \code{by.vars} columns through \code{robust_values} before outputting
#' @param .SDcols advanced; a character vector of column names 
#' passed to inside the data.table's brackets 
#' \code{DT[, , ...]}; see \code{\link{data.table}}; if \code{NULL},
#' uses all appropriate columns. See Examples for usage.
#' @param enclos advanced; an environment; the enclosing
#' environment of the data.
#' @param ... advanced; other arguments passed to inside the 
#' data.table's brackets \code{DT[, , ...]}; see \code{\link{data.table}}
#' 
#' @import data.table 
#' 
#' @details 
#' 
#' Returns \code{expr} for each unique combination of given \code{by.vars}.
#' 
#' By default makes use of any and all \code{\link{levels}} present for 
#' each variable in  \code{by.vars}. This is useful,
#' because even if a subset of the data does not contain observations 
#' for e.g. a specific age group, those age groups are 
#' nevertheless presented in the resulting table; e.g. with the default 
#' \code{expr = list(obs = .N)} all age group levels
#' are represented by a row and can have  \code{obs = 0}.
#' 
#' The function differs from the
#' vanilla \code{\link{table}} by giving a long format table of values
#' regardless of the number of \code{by.vars} given.
#' Make use of e.g. \code{\link{cast_simple}} if data needs to be 
#' presented in a wide format (e.g. a two-way table).
#' 
#' The rows of the long-format table are effectively Cartesian products 
#' of the levels of each variable in  \code{by.vars},
#' e.g. with  \code{by.vars = c("sex", "area")} all levels of  
#' \code{area} are repeated for both levels of  \code{sex}
#' in the table.
#' 
#' The \code{expr} allows the user to apply any function(s) on all 
#' levels defined by  \code{by.vars}. Here are some examples:
#' \itemize{
#'   \item .N or list(.N) is a function used inside a \code{data.table} to 
#'   calculate counts in each group
#'   \item list(obs = .N), same as above but user assigned variable name
#'   \item list(sum(obs), sum(pyrs), mean(dg_age)), multiple objects in a list
#'   \item list(obs = sum(obs), pyrs = sum(pyrs)), same as above with user 
#'   defined var names
#' }
#' 
#' If  \code{use.levels = FALSE}, no \code{levels} information will
#'  be used. This means that if e.g. the  \code{agegroup}
#' variable is a factor and has 18 levels defined, but only 15 levels
#'  are present in the data, no rows for the missing
#' levels will be shown in the table.
#' 
#' \code{na.rm} simply drops any rows from the resulting table where 
#' any of the  \code{by.vars} values was \code{NA}. 
#' 
#' @seealso
#' \code{\link{table}}, \code{\link{cast_simple}}, \code{\link{melt}}
#' 
#' @export ltable
#' 
#' @examples
#' sr <- copy(sire)
#' sr$agegroup <- cut(sr$dg_age, breaks=c(0,45,60,75,85,Inf))
#' ## counts by default
#' ltable(sr, "agegroup")
#' 
#' ## any expression can be given
#' ltable(sr, "agegroup", list(mage = mean(dg_age)))
#' ltable(sr, "agegroup", list(mage = mean(dg_age), vage = var(dg_age)))
#' 
#' ## also returns levels where there are zero rows (expressions as NA)
#' ltable(sr, "agegroup", list(obs = .N, 
#'                             minage = min(dg_age), 
#'                             maxage = max(dg_age)), 
#'        subset = dg_age < 85)
#'        
#' #### expr.by.cj
#' expr.by.cj(sr, "agegroup")
#' 
#' ## any arbitrary expression can be given
#' expr.by.cj(sr, "agegroup", list(mage = mean(dg_age)))
#' expr.by.cj(sr, "agegroup", list(mage = mean(dg_age), vage = var(dg_age)))
#' 
#' ## only uses levels of by.vars present in data
#' expr.by.cj(sr, "agegroup", list(mage = mean(dg_age), vage = var(dg_age)), 
#'            subset = dg_age < 70)
#'            
#' ## .SDcols trick
#' expr.by.cj(sr, "agegroup", lapply(.SD, mean), 
#'            subset = dg_age < 70, .SDcols = c("dg_age", "status"))

ltable <- function(data, 
                   by.vars = NULL, 
                   expr = list(obs = .N), 
                   subset = NULL, 
                   use.levels = TRUE, 
                   na.rm = FALSE,
                   robust = TRUE) {
  
  PF <- parent.frame()
  TF <- environment()
  
  e <- substitute(expr)
  
  ## eval subset ---------------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data, subset, enclos = PF)
  
  ## create table --------------------------------------------------------------
  res <- expr.by.cj(data = data,
                    by.vars = by.vars,
                    expr = e,
                    subset = subset,
                    use.levels = use.levels,
                    na.rm = na.rm,
                    robust = robust)
  
  
  ## final touch ---------------------------------------------------------------
  
  if (!getOption("popEpi.datatable")) {
    setDFpe(res)
  }
  res
  
}




#' @describeIn ltable Somewhat more streamlined \code{ltable} with 
#' defaults for speed. Explicit determination of enclosing environment
#' of data.
#' @export expr.by.cj

expr.by.cj <- function(data, 
                       by.vars = NULL, 
                       expr = list(obs = .N), 
                       subset = NULL, 
                       use.levels = FALSE, 
                       na.rm = FALSE,
                       robust = FALSE,
                       .SDcols = NULL,
                       enclos = parent.frame(1L),
                       ...) {
  
  PF <- enclos
  TF <- environment()
  
  
  ## checks --------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be data.frame (data.table is fine too)")
  }
  
  stopifnot(is.environment(enclos))
  stopifnot(is.logical(na.rm))
  stopifnot(is.logical(use.levels))
  
  stopifnot(is.character(by.vars) || is.null(by.vars))
  all_names_present(data, c(by.vars))
  
  stopifnot(is.character(.SDcols) || is.null(.SDcols))
  all_names_present(data, .SDcols)
  
  tab <- data.table(data[1:min(10, nrow(data)),])
  e <- substitute(expr)
  e <- tab[, evalRecursive(e, env = .SD, enc = PF)$argSub]
  
  ## eval subset ---------------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data, subset, enclos = PF)
  
  ## retrieve data to use without taking copy ----------------------------------
  
  tabVars <- unique(c(by.vars, all.vars(e), .SDcols))
  tabVars <- intersect(names(data), tabVars)
  
  tab <- mget(tabVars, envir = as.environment(data))
  setDT(tab)
  
  tmpDum <- makeTempVarName(data, pre = "dummy_")
  if (!length(by.vars)) {
    if (!length(tab)) {
      ## no by.vars nor variables in expr
      tab <- data.table(rep(1L, nrow(data)))
      setnames(tab, "V1", tmpDum)
    } else {
      tab[, c(tmpDum) := 1L]
    }
    by.vars <- tmpDum
  }
  
  ## create joining table ------------------------------------------------------
  lev_fun <- function(x) {
    if (use.levels && is.factor(x)) {
      factor(levels(x), levels = levels(x))
    } else {
      sort(unique(x))
    }
  }
  
  cj <- lapply(as.list(tab)[by.vars], lev_fun)
  cj <- do.call(CJ, c(cj, unique = FALSE, sorted = FALSE))
  if (na.rm) cj <- na.omit(cj)
  
  ## eval expression -----------------------------------------------------------
  exprVars <- setdiff(names(tab), by.vars)
  if (!length(.SDcols)) .SDcols <- exprVars
  if (!length(.SDcols)) .SDcols <- tabVars
  res <- tab[subset][cj, eval(e, envir = .SD), 
                     on = by.vars, 
                     by = .EACHI, 
                     .SDcols = .SDcols, ...]
  
  setcolsnull(res, delete = tmpDum, soft = TRUE)
  by.vars <- setdiff(by.vars, tmpDum)
  
  ## final touch ---------------------------------------------------------------
  if (length(res)) setcolorder(res, c(by.vars, setdiff(names(res), by.vars)))
  if (length(by.vars)) setkeyv(res, by.vars)
  if (!getOption("popEpi.datatable")) {
    setDFpe(res)
  }
  res
  
}

