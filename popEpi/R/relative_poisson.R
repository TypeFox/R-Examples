#' @title Excess hazard Poisson model
#' @author Joonas Miettinen, Karri Seppa
#' @description Estimate a Poisson piecewise constant excess
#' hazards model
#' @param data a dataset splitted with e.g. \code{\link{lexpand}};
#' must have expected hazard merged within
#' @param formula a formula which is passed on to \code{glm}; see Details
#' @param fot.breaks optional; a numeric vector of [a,b) breaks to specify
#' survival intervals over the follow-up time; if \code{NULL}, the 
#' existing breaks along the mandatory \code{fot} time scale in \code{data}
#' are used (e.g. the breaks for \code{fot} supplied to \code{lexpand})
#' @param subset a logical vector or condition; e.g. \code{subset = sex == 1};
#' limits the data before estimation
#' @param check logical; if \code{TRUE}, tabulates excess cases by all
#' factor variables in the formula to check for negative / \code{NA} 
#' excess cases before fitting the glm
#' @param ... any argument passed on to \code{glm}
#' @import stats
#' @details
#' 
#' \strong{Basics}
#' 
#' \code{relpois} employs a custom link function of the Poisson variety
#' to estimate piecewise constant parametric excess hazards. The pieces
#' are determined by \code{fot.breaks}. A \code{log(person-years)} offset
#' is passed automatically to the \code{glm} call.
#' 
#' \strong{Formula usage}
#' 
#' The formula can be used like any ordinary \code{glm} formula. The user must
#' define the outcome in some manner, which is usually \code{lex.Xst} after splitting
#' with e.g. \code{lexpand}. The exception is the possibility of including 
#' the baseline excess hazard terms by including the 
#' reserved term \code{FOT} in the formula.
#' 
#' For example, \code{lex.Xst != 0 ~ FOT + agegr} estimates a model with constant
#' excess hazards at the follow-up intervals as specified by 
#' the pertinent breaks used in splitting \code{data},
#' as well as for the different age groups.
#' \code{FOT} is created ad hoc if it is used in the formula.
#' If you leave out \code{FOT}, the hazard is effectively
#' assumed to be constant across the whole follow-up time. 
#' 
#' You can also simply use your own follow-up time interval variable that
#' you have created before calling \code{relpois}. However, when using 
#' \code{FOT}, \code{relpois} automatically checks for e.g. 
#' negative excess cases in follow-up intervals,
#' allowing for quickly finding splitting breaks
#' where model estimation is possible. It also drops any data outside the
#' follow-up time window.
#' 
#' \strong{Splitting and merging population hazard}
#' 
#' The easiest way to both split and to include population hazard information is 
#' by using \code{\link{lexpand}}. You may also fairly easily do it by hand
#' by splitting first and then merging in your population hazard information.
#' 
#' 
#' \strong{Data requirements}
#' 
#' The population hazard information must be available for each record and named
#' \code{pop.haz}. The follow-up time variable must be named \code{"fot"} e.g.
#' as a result of using \code{lexpand}. The \code{lex.dur} variable must also
#' be present, containing person-year information. 
#' 
#' 
#' @return
#' A \code{glm} object created using a custom Poisson family construct. Some
#' \code{glm} methods are applicable.
#' 
#' @seealso
#' \code{\link{lexpand}}, \code{\link{poisson}}, \code{\link{glm}}
#' 
#' @export relpois
#' 
#' @examples
#' ## use the simulated rectal cancer cohort
#' sr <- copy(sire)
#' sr$agegr <- cut(sr$dg_age, c(0,45,60,Inf), right=FALSE)
#' 
#' ## usable straight away after splitting
#' fb <- c(0,3/12,6/12,1,2,3,4,5)
#' x <- lexpand(sr, birth = bi_date, entry = dg_date,
#'              exit = ex_date, status=status,
#'              breaks = list(fot=fb), pophaz=popmort)
#' rpm <- relpois(x, formula = lex.Xst %in% 1:2 ~ FOT + agegr)
#'  
#' ## some methods for glm work. e.g. test for interaction
#' # rpm2 <- relpois(x, formula = lex.Xst %in% 1:2 ~ FOT*agegr)
#' # anova(rpm, rpm2, test="LRT")
#' # AIC(rpm, rpm2)

relpois <- function(data, 
                    formula, 
                    fot.breaks = NULL, subset = NULL, check=TRUE, ...) {
  ## prep arguments ------------------------------------------------------------
  excess_cases <- fot <- pop.haz <-  NULL ## appease R CMD CHECK
  
  ## somehow the class of the data is being altered by this function
  oldClass <- class(data)
  
  if (missing(formula) || !inherits(formula, "formula")) stop("formula not defined")
  
  form_vars <- all.vars(formula)
  dataname <- as.name(deparse(substitute(data)))
  
  if (!inherits(data, "Lexis")) {
    stop("data is not a Lexis object; data must be a result of splitting or using Lexis")
  }
  
  if ("FOT" %in% names(data)) {
    stop("FOT is a reserved name but you have a variable with that name in data; rename/delete it first")
  }
  
  # wasDF <- FALSE
  if (!is.data.table(data)) {
    data <- copy(data)
    setDT(data)
    message("Took a copy of your data because it was a data.frame and not a data.table. This may take up a lot of memory.")
    message("It is recommended to convert your data to data.table before using this function using as.data.table or setDT")
  }
  
  req_vars <- unique(c("lex.id", "fot", "lex.dur", "pop.haz", setdiff(form_vars, "FOT")))
  all_names_present(data, req_vars)
  
  surv.breaks <- attr(data, "breaks")$fot
  if (is.null(surv.breaks)) {
    stop("did not find any breaks information in data attributes named 'fot';
            probable reason: splitted data was edited after splitting - don't do that")
  } else {
    if (!is.null(fot.breaks)) {
      if (any(!fot.breaks %in% surv.breaks)) {
        stop("fot.breaks must be a subset of the breaks for 'fot' used in splitting;
           type attr(data, 'breaks')$fot to see the breaks you used in splitting")
      } else {
        surv.breaks <- fot.breaks
      }
    }
    
  }
  
  ## prep & subset data --------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data, subset)
  
  if (any(is.na(data[subset, ]$pop.haz))) {
    stop("some pop.haz are NA")
  }
  
  on.exit({
    setcolsnull(data, c("FOT", tmpdexp), soft = TRUE)
  }, add = TRUE)
  
  if ("FOT" %in% form_vars)  {
    data[, "FOT" := cut(fot, breaks = surv.breaks, right = FALSE)]
    # set(data, j = "FOT", value = cut(data$fot, breaks = surv.breaks, right=FALSE))
    subset <- subset & !is.na(data$FOT)
  }
  tmpdexp <- makeTempVarName(data, pre = "TEMP_d.exp_")
  data[, c(tmpdexp) := pop.haz*lex.dur]
  # set(data, j = tmpdexp, value = data$pop.haz * data$lex.dur)
  
  
  if (check) {
    ## test for negative excess cases in factor variable combinations ------------
    
    ## determine factor variables for cross-tabulating
    
    fac_vars <- colnames(attr(terms.formula(formula), "factors"))
    
    fac_list <- paste0(fac_vars, collapse=", ")
    fac_list <- paste0("list(", fac_list, ")")
    fac_list <- parse(text=fac_list)
    
    wh_fac <- as.data.table(data)[subset, unlist(lapply(eval(fac_list), is.factor))]
    fac_vars <- fac_vars[wh_fac]
    
    if (length(fac_vars) == 0) fac_vars <- NULL
    fac_list <- paste0(fac_vars, collapse=", ")
    fac_list <- paste0("list(", fac_list, ")")
    fac_list <- parse(text=fac_list)
    
    
    ## test negativity of excess cases
    LHS <- as.character(formula)
    LHS <- LHS[2]
    LHS <- parse(text = LHS)
    
    excas <- as.data.table(data)[subset, list(excess_cases = sum(eval(LHS)-get(tmpdexp))), keyby=eval(fac_list)]
    setnames(excas, 1:ncol(excas), c(fac_vars, "excess_cases"))
    
    
    if (any(is.na(excas$excess_cases))) {
      stop("some excess cases were NA; is pop.haz available for all records?")
    }
    excas <- excas[excess_cases <= 0]
    if (any(excas$excess_cases<=0)) {
      print(excas)
      warning("negative excess cases found in some combinations of factor variables;
           see printed table and try e.g. wider FOT intervals")
    }
  }
  
  
  ## custom poisson family -----------------------------------------------------
  RPL <- copy(poisson())
  RPL$link <- "glm relative survival model with Poisson error"
  RPL$linkfun <- function(mu, d.exp = data[subset, ][[tmpdexp]]) log(mu - d.exp)
  RPL$linkinv <- function(eta, d.exp = data[subset, ][[tmpdexp]]) d.exp + exp(eta)
  
  
  init_maker <- quote( {
    if (any(y < 0)) stop(paste("Negative values not allowed for", 
                               "the Poisson family"))
    n <- rep.int(1, nobs)
    mustart <- pmax(y, rep(1,times=length(y))) + 0.1
  } )
  
  RPL$initialize <- init_maker
  
  
  ## glm call ------------------------------------------------------------------
  ## update() won't work
  ## anova() works
  
  
  ml <- glm(formula = formula, data=data[subset,], offset=log(lex.dur),
            #             subset = subset, ## Error in xj[i] : invalid subscript type 'closure' 
            family = RPL, ...)
  
  ## final touches -------------------------------------------------------------
  ml$d.exp <- data[subset, ][[tmpdexp]]
  ml$FOT <- data[subset, ]$FOT
  ml$fot.breaks <- surv.breaks
  ml$call$data <- dataname
  ml$call$formula <- formula
  setattr(ml, "class", c("relpois", "glm", "lm"))
  
  setattr(data, "class", oldClass) ## see beginning of function
  
  
  ml
}



#' @title Excess hazard Poisson model
#' @author Joonas Miettinen, Karri Seppa
#' @description Estimate a Poisson Piecewise Constant Excess
#' Hazards Model
#' @param formula a formula with the counts of events as the response.
#' Passed on to \code{glm}. May contain usage of the \code{offset()} function
#' instead of supplying the offset for the Poisson model via the argument
#' \code{offset}.
#' @param data an \code{aggre} object (an aggregated data set; 
#' see \code{\link{as.aggre}} and \code{\link{aggre}})
#' @param d.exp the counts of expected cases. Mandatory.
#' E.g. \code{d.exp = EXC_CASES}, where \code{EXC_CASES} is a column in data.
#' @param offset the offset for the Poisson model, supplied as e.g.
#' \code{offset = log(PTIME)}, where \code{PTIME} is a subject-time
#' variable in data. Not mandatory, but almost always should be supplied.
#' @param breaks optional; a numeric vector of [a,b) breaks to specify
#' survival intervals over the follow-up time; if \code{NULL}, the 
#' existing breaks along the mandatory time scale mentioned in \code{formula}
#' are used
#' @param subset a logical vector or condition; e.g. \code{subset = sex == 1};
#' limits the data before estimation
#' @param piecewise \code{logical}; if \code{TRUE}, and if any time scale
#' from data is used (mentioned) in the formula, the time scale is 
#' transformed into a factor variable indicating intervals on the time scale.
#' Otherwise the time scale left as it is, usually a numeric variable.
#' E.g. if \code{formula = counts ~ TS1*VAR1}, \code{TS1} is transformed
#' into a factor before fitting model.
#' @param ... any other argument passed on to \code{\link[stats]{glm}} such as 
#' \code{control} or \code{weights}
#' @import stats
#' 
#' @return
#' A \code{relpois} object created using a custom Poisson family construct.
#' 
#' @seealso
#' \code{\link{lexpand}}, \code{\link{poisson}}, \code{\link{glm}}
#' 
#' @examples
#' ## use the simulated rectal cancer cohort
#' sr <- copy(sire)
#' sr$agegr <- cut(sr$dg_age, c(0,45,60,Inf), right=FALSE)
#' 
#' ## usable straight away after splitting
#' fb <- c(0,3/12,6/12,1,2,3,4,5)
#' x <- lexpand(sr, birth = bi_date, entry = dg_date,
#'              exit = ex_date, status=status %in% 1:2,
#'              breaks = list(fot=fb), 
#'              pophaz=popmort, pp = FALSE,
#'              aggre = list(agegr, fot))
#'              
#' rpm <- relpois_ag(formula = from0to1 ~ fot + agegr,  data = x,
#'                   d.exp = d.exp, offset = log(pyrs))
#' summary(rpm)
#'  
#' ## the usual functions for handling glm models work
#' rpm2 <- update(rpm, . ~ fot*agegr)
#' anova(rpm, rpm2, test="LRT")
#' AIC(rpm, rpm2)
#' @export

relpois_ag <- function(formula, data, d.exp, offset = NULL, breaks = NULL, subset = NULL, piecewise = TRUE, ...) {
  
  TF <- environment()
  PF <- parent.frame(1L)
  original_formula <- formula
  
  if (!inherits(data, "aggre")) stop("data is not an aggre object. Please aggregate your data first using e.g. lexpand(). If your data is aggregated, use as.aggre() to mark it as such.")
  
  formula <- evalRecursive(formula, env = TF, enc = PF)$arg
  if (missing(formula) || !inherits(formula, "formula")) stop("formula not defined")
  
  ## detect survival time scale ------------------------------------------------
  oldBreaks <- attr(data, "breaks")
  if (is.null(oldBreaks)) stop("data does not have breaks information. Is it a result of using aggre() or as.aggre()?")
  survScale <- intersect(all.vars(formula), names(oldBreaks))
  if (length(survScale) > 1L) stop("Found several used time scales in formula, which is not supported (found ", paste0("'", survScale, "'", collapse = ", "), ")")
  
  ## pre-find args -------------------------------------------------------------
  sub_d.exp <- evalRecursive(substitute(d.exp), env = data[1L, ], enc = PF)$argSub
  sub_offset <- evalRecursive(substitute(offset), env = data[1L, ], enc = PF)$argSub
  
  ## prep & subset data --------------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data, subset)
  
  av <- c(all.vars(formula), all.vars(sub_d.exp), all.vars(sub_offset))
  av <- intersect(names(data), av)
  
  x <- subsetDTorDF(data, subset = subset, select = av)
  
  setDT(x)
  setattr(x, "class", c("aggre", "data.table", "data.frame"))
  
  ## handle breaks -------------------------------------------------------------
  if (!is.null(breaks)) {
    if (!piecewise) stop("Supplied breaks but piecewise = FALSE. Please select piecewise = TRUE if you want piecewise estimates defined by the breaks.")
    breaks <- sort(unique(breaks))
    if (!all(breaks %in% oldBreaks[[survScale]])) stop("Supplied breaks are not a subset of the breaks in data. See the breaks in data using e.g. attr(data, 'breaks')")
    
  } else breaks <- oldBreaks[[survScale]]
  
  if (piecewise && length(survScale) > 0L) {
    breaks <- breaks - .Machine$double.eps^0.5 ## ensures cutting correctness
    
    set(x, j = survScale, value = cut(x[[survScale]], breaks = breaks, right = FALSE, labels = FALSE))
    
    breaks <- round(breaks, 2L)
    pieces <- paste0("[", breaks[-length(breaks)], ", ", breaks[-1L], ")")
    set(x, j = survScale, value = pieces[x[[survScale]]])
    
  }
  
  ## eval value args -----------------------------------------------------------
  d.exp <- evalPopArg(x, sub_d.exp, enclos = PF, 
                      DT = TRUE, recursive = TRUE)
  if (is.null(d.exp)) stop("argument d.exp was not supplied")
  d.exp <- rowSums(d.exp)
  if (length(d.exp) == nrow(data)) d.exp <- d.exp[subset]
  
  offset <- evalPopArg(x, sub_offset, enclos = PF, 
                       DT = TRUE, recursive = TRUE)
  if (!is.null(offset)) offset <- rowSums(offset)
  if (length(offset) == nrow(data)) offset <- offset[subset]
  
  
  ## custom poisson family -----------------------------------------------------
  RPL <- copy(poisson())
  RPL$link <- "glm relative survival model with Poisson error"
  RPL$linkfun <- function(mu, d.exp = TF$d.exp) {
    log(mu - d.exp)
  }
  RPL$linkinv <- function(eta, d.exp = TF$d.exp) {
    d.exp + exp(eta)
  }
  
  
  init_maker <- quote( {
    if (any(y < 0)) stop(paste("Negative values not allowed for", 
                               "the Poisson family"))
    n <- rep.int(1, nobs)
    mustart <- pmax(y, rep(1,times=length(y))) + 0.1
  } )
  
  RPL$initialize <- init_maker
  
  
  ## glm call ------------------------------------------------------------------
  
  ## NOTE: parent.frame(3L) to find this (this function's) environment
  ml <- glm(formula = formula, data=x, offset = parent.frame(3L)$offset, family = RPL, ...)
  
  ## final touches -------------------------------------------------------------
  
  ml$call <- match.call()
  setattr(ml, "class", c("relpois", "glm", "lm"))
  
  ml
}





