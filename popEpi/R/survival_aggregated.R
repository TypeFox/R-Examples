#' @title Estimate Survival Time Functions
#' @aliases survtab
#' @author Joonas Miettinen, Karri Seppa
#' @description Estimate survival time functions (survival, relative/net
#' survival, CIFs) using aggregated data (\code{survtab_ag}) or 
#' \code{\link[Epi]{Lexis}} data (\code{survtab}).
#' @param formula a \code{formula} object. For \code{survtab_ag}, the response 
#' must be the time scale to compute survival time function estimates
#' over, e.g. \code{fot ~ sex + adjust(agegr)}. For \code{survtab} it can
#' also be the just e.g. \code{fot ~ sex + adjust(agegr)}, whereupon it is
#' assumed that \code{lex.Xst} in your data is the status variable in the
#' intended format. To be explicit, use \code{\link[survival]{Surv}}
#' containing the 
#' survival time scale and the event indicator as the right-hand-side, e.g. 
#' \code{Surv(fot, lex.Xst) ~ sex}.
#' Stratifying variables can be included on the right-hand side
#' separated by '\code{+}'. May contain usage of \code{adjust()} 
#' --- see Details.
#' @param data for \code{survtab_ag}, a data set of aggregated counts, 
#' subject-times, etc., as created
#' by \code{\link{aggre}}; for pre-aggregated data see \code{\link{as.aggre}}.
#' For \code{survtab}, a \code{Lexis} object.
#' @param adjust can be used as an alternative to passing variables to 
#' argument \code{formula} within a call to \code{adjust()}; e.g.
#' \code{adjust = "agegr"}. \link[=flexible_argument]{Flexible input}.
#' @param weights typically a list of weights or a \code{character} string
#' specifying an age group standardization scheme; see
#' the \link[=direct_standardization]{dedicated help page} 
#' and examples.
#' @param surv.breaks (\code{survtab_ag} only) a vector of breaks on the survival time scale. Only used
#' to compute estimates using a subset of the intervals in data or larger intervals
#' than in data. E.g. one might use \code{surv.breaks = 0:5} when the aggregated
#' data has intervals with the breaks \code{seq(0, 10, 1/12)}.
#' @param breaks (\code{survtab} only) a named list of breaks, e.g.
#' \code{list(FUT = 0:5)}. If data is not split in advance, \code{breaks}
#' must at the very least contain a vector of breaks to split the survival time 
#' scale (mentioned in argument \code{formula}). If data has already been split
#' (using e.g. \code{\link{splitMulti}}) along at least the used survival time
#' scale, this may be \code{NULL}.
#' @param pophaz (\code{survtab} only) a \code{data.frame} containing
#' expected hazards for the event of interest to occur. See the
#' \link[=pophaz]{dedicated help page}. Required when
#' \code{surv.type = "surv.rel"} or \code{"cif.rel"}. \code{pophaz} must
#' contain one column named \code{"haz"}, and any number of other columns
#' identifying levels of variables to do a merge with split data within
#' \code{survtab}. Some columns may be time scales, which will
#' allow for the expected hazard to vary by e.g. calendar time and age.
#' 
#' @param n variable containing counts of subjects at-risk at the start of a time 
#' interval; e.g. \code{n = "at.risk"}. 
#' Required when \code{surv.method = "lifetable"}.
#' \link[=flexible_argument]{Flexible input}.

#' @param d variable(s) containing counts of subjects experiencing an event. With
#' only one type of event, e.g. \code{d = "deaths"}. With multiple types of 
#' events (for CIF or cause-specific survival estimation), supply e.g.
#' \code{d = c("canD", "othD")}. If the survival time function to be estimated
#' does not use multiple types of events, supplying more than one variable
#' to \code{d} simply causes the variables to be added together. 
#' Always required. \link[=flexible_argument]{Flexible input}.

#' @param n.cens variable containing counts of subjects censored during a 
#' survival time interval; E.g. \code{n.cens = "alive"}.
#' Required when \code{surv.method = "lifetable"}. 
#' \link[=flexible_argument]{Flexible input}.

#' @param pyrs variable containing total subject-time accumulated within a survival
#' time interval; E.g. \code{pyrs = "pyrs"}. 
#' Required when \code{surv.method = "hazard"}. Flexible input.

#' @param d.exp variable denoting total "expected numbers of events" 
#' (typically computed \code{pyrs * pop.haz}, where 
#' \code{pop.haz} is the expected hazard level) 
#' accumulated within a survival time interval; E.g. \code{pyrs = "pyrs"}.
#' Required when computing EdererII relative survivals or 
#' CIFs based on excess counts of events. Flexible input.

#' @param n.pp variable containing total Pohar-Perme weighted counts of
#' subjects at risk in an interval,
#' supplied as argument \code{n} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == "at-risk")}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.
#' 
#' @param d.pp variable(s) containing total Pohar-Perme weighted counts of events,
#' supplied as argument \code{d} is supplied. Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == some_event)}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param d.pp.2 variable(s) containing total Pohar-Perme 
#' "double-weighted" counts of events,
#' supplied as argument \code{d} is supplied. Computed originally on the subject
#' level as analogous to \code{pp * pp * as.integer(status == some_event)}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param n.cens.pp variable containing total Pohar-Perme weighted counts 
#' censorings,
#' supplied as argument \code{n.cens} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == "censored")}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param pyrs.pp variable containing total Pohar-Perme weighted subject-times,
#' supplied as argument \code{pyrs} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * pyrs}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param d.exp.pp variable containing total Pohar-Perme weighted counts 
#' of excess events,
#' supplied as argument \code{pyrs} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * d.exp}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.
#' 
#' @param surv.type one of \code{'surv.obs'},
#' \code{'surv.cause'}, \code{'surv.rel'}, 
#' \code{'cif.obs'} or \code{'cif.rel'}; 
#' defines what kind of survival time function(s) is/are estimated; see Details
#' @param surv.method either \code{'lifetable'} or \code{'hazard'}; determines
#' the method of calculating survival time functions, where the former computes
#' ratios such as \code{p = d/(n - n.cens)} 
#' and the latter utilizes subject-times 
#' (typically person-years) for hazard estimates such as \code{d/pyrs} 
#' which are used to compute survival time function estimates.
#' The former method requires argument \code{n.cens} and the latter 
#' argument \code{pyrs} to be supplied. 
#' @param relsurv.method  either \code{'e2'} or \code{'pp'}; 
#' defines whether to compute relative survival using the
#' EdererII method or using Pohar-Perme weighting;
#' ignored if \code{surv.type != "surv.rel"}
#'  
#' @param subset a logical condition; e.g. \code{subset = sex == 1}; 
#' subsets the data before computations
#' 
#' @param conf.level confidence level used in confidence intervals; 
#' e.g. \code{0.95} for 95 percent confidence intervals
#' @param conf.type character string; must be one of \code{"plain"}, 
#' \code{"log-log"} and \code{"log"}; 
#' defines the transformation used on the survival time
#' function to yield confidence 
#' intervals via the delta method
#' @param verbose logical; if \code{TRUE}, the function is chatty and
#'  returns some messages and timings along the process
#' 
#' @details
#' 
#' \strong{Basics}
#' 
#' \code{survtab_ag} computes estimates of survival time functions using 
#' pre-aggregated data. When using subject-level data directly, use 
#' \code{survtab}. For aggregating data, see \code{\link{lexpand}}
#' and \code{\link{aggre}}. Data sets can be coerced into \code{Lexis}
#' objects using \code{\link[Epi]{Lexis}}.
#' 
#' By default
#' \code{survtab_ag} makes use of the exact same breaks that were used in 
#' splitting the original data (with e.g. \code{lexpand}), so it is not 
#' necessary to specify any
#' \code{surv.breaks}. If specified, the 
#' \code{surv.breaks} must be a subset of the pertinent 
#' pre-existing breaks. Interval lengths (\code{delta} in output) are 
#' also calculated based on existing breaks or \code{surv.breaks}, 
#' so the upper limit of the breaks should
#' therefore be meaningful and never e.g. \code{Inf}. 
#' 
#' \code{survtab} may be a split or unsplit \code{Lexis} data set, but it 
#' is recommended to supply the \code{breaks} argument anyway.
#' 
#' if \code{surv.type = 'surv.obs'}, only 'raw' observed survival 
#' is calculated over the chosen time intervals. With
#' \code{surv.type = 'surv.rel'}, also relative survival estimates 
#' are supplied in addition to observed survival figures. 
#' 
#' \code{surv.type = 'cif.obs'} requests cumulative incidence functions (CIF) 
#' to be estimated, where all unique events (supplied via \code{d})
#' are seen as competing risks. 
#' CIFs are estimated for each competing risk based 
#' on a survival-interval-specific proportional hazards
#' assumption as described by Chiang (1968).  
#' With \code{surv.type = 'cif.rel'}, a CIF is estimated with using 
#' excess cases as the ''cause-specific'' cases.
#' 
#' if \code{surv.type = 'surv.cause'}, cause-specific survivals are estimated
#' separately for each unique value of \code{event.values}.
#' 
#' The vignette \href{../doc/survtab_examples.html}{survtab_examples} 
#' has some practical examples.
#' 
#' \strong{Relative / net survival}
#'  
#' When \code{surv.type = 'surv.rel'}, the user can choose 
#' \code{relsurv.method = 'pp'}, whereupon
#' Pohar-Perme weighting is used.
#' By default \code{relsurv.method = 'e2'}.
#'
#' \strong{Adjusted estimates}
#' 
#' Adjusted estimates in this context mean computing estimates separately
#' by the levels of adjusting variables and returning weighted averages
#' of the estimates. For example, computing estimates separately by
#' age groups and returning a weighted average estimate (age-adjusted estimate).
#' 
#' Adjusting requires specification of both the adjusting variables and
#' the weights for all the levels of the adjusting variables. The former can be
#' accomplished by using \code{adjust()} with the argument \code{formula},
#' or by supplying variables directly to argument \code{adjust}. E.g. the
#' following are all equivalent (using \code{survtab_ag}):
#' 
#' \code{formula = fot ~ sex + adjust(agegr, area)}
#' 
#' \code{formula  = fot ~ sex} and \code{adjust = c("agegr", "area")}
#' 
#' \code{formula  = fot ~ sex} and \code{adjust = list(agegr, area)}
#' 
#' When using \code{survtab}, the response must be a 
#' \code{\link[survival]{Surv}} object, e.g. 
#' \code{Surv(time = fot, event = lex.Xst)}, but otherwise the 
#' syntax is the same.
#' 
#' The adjusting variables must match with the variable names in the
#' argument \code{weights}, which may be supplied as a \code{list} or
#' a \code{data.frame}. The former can be done by e.g.
#' 
#' \code{weights = list(agegr = VEC1, area = VEC2)},
#' 
#' where \code{VEC1} and \code{VEC2} are vectors of weights (which do not
#' have to add up to one). When passing a \code{list} of weights, the order
#' of the weights must match with the order of the levels of the variable:
#' For \code{factor} variables, they must correspond to the \code{factor}'s
#' levels. Otherwise they must match to the sorted levels of the variable,
#' i.e. if variable \code{agegr} has three levels \code{c(1, 2, 3)},
#' the weights in e.g. \code{list(agegr = c(0.1, 0.4, 0.5))} would pass 
#' the weight
#' \code{0.1} for level \code{1} and so on, regardless of the order of values in
#' the data.
#' 
#' 
#' See 
#' \href{../doc/survtab_examples.html}{survtab_examples} 
#' for an example of using a \code{data.frame} to pass weights.
#' 
#' 
#' 
#' \strong{Period analysis and other data selection schemes}
#' 
#' If one wishes to calculate e.g. period analysis (delayed entry estimates), 
#' one should limit the data when/before supplying to this function. 
#' See e.g. \code{\link{lexpand}}.
#' 
#' @return
#' Returns a table of life time function values and other 
#' information with survival intervals as rows.
#' Returns some of the following estimates of survival time functions:
#' 
#' \itemize{
#'  \item \code{surv.obs} - observed (raw) survival
#'  \item \code{CIF_k} - cumulative incidence function for cause \code{k}
#'  \item \code{CIF.rel} - cumulative incidence function using excess cases
#'  \item \code{r.e2} -  relative survival, EdererII
#'  \item \code{r.pp} -  relative survival, Pohar-Perme weighted
#' }
#' The suffix \code{.as} implies adjusted estimates, and \code{.lo} and
#' \code{.hi} imply lower and upper confidence limits, respectively. 
#' The prefix \code{SE.} stands for standard error.
#' 
#' @import data.table
#' 
#' @seealso
#' \code{\link{splitMulti}}, \code{\link{lexpand}}, 
#' \code{\link{ICSS}}, \code{\link{sire}},
#' \code{\link{survtab}}
#' \href{../doc/survtab_examples.html}{The survtab_examples vignette}
#' 
#' @references
#' 
#' Perme, Maja Pohar, Janez Stare, and Jacques Estève. 
#' "On estimation in relative survival." Biometrics 68.1 (2012): 113-120.
#' 
#' Hakulinen, Timo, Karri Seppa, and Paul C. Lambert. 
#' "Choosing the relative survival method for cancer survival estimation." European Journal of Cancer 47.14 (2011): 2202-2210.
#'  
#' Seppa, Karri, Timo Hakulinen, and Arun Pokhrel. 
#' "Choosing the net survival method for cancer survival estimation." European Journal of Cancer (2013).
#' 
#' CHIANG, Chin Long. Introduction to stochastic processes in biostatistics. 1968.
#'  
#' 
#' @examples
#' ## see more examples with explanations in vignette("survtab_examples")
#' 
#' #### survtab_ag usage
#' 
#' data(sire)
#' ## prepare data for e.g. 5-year "period analysis" for 2008-2012
#' ## note: sire is a simulated cohort integrated into popEpi.
#' BL <- list(fot=seq(0, 5, by = 1/12),
#'            per = c("2008-01-01", "2013-01-01"))
#' x <- lexpand(sire, birth = bi_date, entry = dg_date, exit = ex_date,
#'              status = status %in% 1:2,
#'              breaks = BL,
#'              pophaz = popmort,
#'              aggre = list(fot))
#'              
#' ## calculate relative EdererII period method 
#' st <- survtab_ag(fot ~ 1, data = x)
#' 
#' #### survtab usage
#' library(Epi)
#' library(survival)
#'
#' ## NOTE: recommended to use factor status variable
#' x <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)), 
#'            exit = list(CAL = get.yrs(ex_date)), 
#'            data = sire[sire$dg_date < sire$ex_date, ],
#'            exit.status = factor(status, levels = 0:2, 
#'                                 labels = c("alive", "canD", "othD")), 
#'            merge = TRUE)
#' 
#' ## phony group variable
#' set.seed(1L)
#' x$group <- rbinom(nrow(x), 1, 0.5)
#' 
#' \dontrun{
#' ## observed survival. explicit supplying of status:
#' st <- survtab(Surv(time = FUT, event = lex.Xst) ~ group, data = x, 
#'               surv.type = "surv.obs",
#'               breaks = list(FUT = seq(0, 5, 1/12)))
#' ## this assumes the status is lex.Xst (right 99.9 % of the time)
#' st <- survtab(FUT ~ group, data = x, 
#'               surv.type = "surv.obs",
#'               breaks = list(FUT = seq(0, 5, 1/12)))
#'
#'#### using dates with survtab
#' x <- Lexis(entry = list(FUT = 0L, AGE = dg_date-bi_date, CAL = dg_date),
#'            exit = list(CAL = ex_date),
#'            data = sire[sire$dg_date < sire$ex_date, ],
#'            exit.status = factor(status, levels = 0:2, 
#'                                 labels = c("alive", "canD", "othD")), 
#'            merge = TRUE)
#' ## phony group variable
#' set.seed(1L)
#' x$group <- rbinom(nrow(x), 1, 0.5)
#' 
#' st <- survtab(Surv(time = FUT, event = lex.Xst) ~ group, data = x, 
#'               surv.type = "surv.obs",
#'               breaks = list(FUT = seq(0, 5, 1/12)*365.25))    
#'                   
#' ## NOTE: population hazard should be reported at the same scale
#' ## as time variables in your Lexis data.
#' data(popmort, package = "popEpi")
#' pm <- data.frame(popmort)
#' names(pm) <- c("sex", "CAL", "AGE", "haz")
#' ## from year to day level
#' pm$haz <- pm$haz/365.25 
#' pm$CAL <- as.Date(paste0(pm$CAL, "-01-01")) 
#' pm$AGE <- pm$AGE*365.25 
#' 
#' st <- survtab(Surv(time = FUT, event = lex.Xst) ~ group, data = x, 
#'               surv.type = "surv.rel", relsurv.method = "e2",
#'               pophaz = pm,
#'               breaks = list(FUT = seq(0, 5, 1/12)*365.25))  
#' }
#' @export
survtab_ag <- function(formula = NULL,
                       
                       data, 
                       
                       adjust = NULL,
                       weights = NULL,
                       
                       surv.breaks = NULL, 
                       
                       n = "at.risk",
                       d = "from0to1",
                       n.cens = "from0to0",
                       pyrs = "pyrs",
                       d.exp = "d.exp",
                       
                       n.pp = NULL,
                       d.pp = "d.pp",
                       d.pp.2 = "d.pp.2",
                       n.cens.pp = "n.cens.pp",
                       pyrs.pp = "pyrs.pp",
                       d.exp.pp = "d.exp.pp",
                       
                       surv.type="surv.rel", 
                       surv.method="hazard", 
                       relsurv.method="e2",  
                       
                       subset = NULL,
                       
                       conf.level = 0.95, 
                       conf.type = "log-log",
                       
                       verbose=FALSE) {
  
  if (verbose) starttime <- proc.time()
  
  Tstop <- delta <- Tstart <- surv.int <- n.eff <- n.eff.pp <- surv.obs <- 
     lag1_surv.obs <- p.obs <- CIF.rel <- NULL ## APPEASE R CMD CHECK
  
  TF <- environment()
  PF <- parent.frame(1L)
  
  this_call <- match.call()
  used_args <- as.list(this_call)[-1L]
  fl <- formals("survtab_ag")
  used_args <- c(used_args, fl[!names(fl) %in% names(used_args)])
  used_args <- used_args[names(fl)]
  rm(fl)
  
  # check data -----------------------------------------------------------------
  if (missing(data) || nrow(data) == 0) stop("data missing or has no rows")
  
  if (!inherits(data, "aggre")) stop("Data must be an aggre object; see ?aggre")
  
  # check arguments ------------------------------------------------------------
  
  surv.type <- match.arg(surv.type, c("surv.obs","surv.rel","surv.cause", "cif.obs", "cif.rel"))
  surv.method <- match.arg(surv.method, c("lifetable","hazard"))
  relsurv.method <- match.arg(relsurv.method, c("e2", "pp", "EdererII", "Pohar-Perme", "pohar-perme", "edererII", "ederer2"))
  if (relsurv.method %in% c("EdererII", "edererII", "ederer2")) relsurv.method <- "e2"
  if (relsurv.method %in% c("Pohar-Perme", "pohar-perme")) relsurv.method <- "pp"
  relsurv.method <- match.arg(relsurv.method, c("e2", "pp"))
  conf.type <- match.arg(conf.type, c("log","log-log","plain"))
  
  
  # handle breaks in attributes ------------------------------------------------
  
  found_breaks <- NULL
  attrs <- attributes(data)
  if (is.null(attrs$breaks)) {
    stop("Data does not have breaks information and surv.breaks not defined; ",
         "this would be true if data is output from aggre() or lexpand(). ",
         "If it is and you did not tamper with it, complain to the ",
         "package maintainer.")
  } 
  
  breakScales <- names(attrs$breaks)
  
  ## argument 'formula' pre-check ----------------------------------------------
  if (!inherits(formula, "formula")) {
    stop("Argument 'formula' does not appear to be a formula object. ",
         "Usage: e.g. fot ~ sex")
  }
  if (length(formula) != 3L) {
    stop("formula does not appear to be two-sided; supply it as e.g. fot ~ sex")
  }
  surv.scale <- deparse(formula[[2L]])
  if (!surv.scale %in% names(data)) {
    stop("Supplied time scale '", surv.scale, "' is not a name of a time ",
         "scale by which data has been aggregated (no column with ",
         "that name in data)")
  }
  if (!surv.scale %in% breakScales) {
    stop("Supplied time scale '", surv.scale, "' is not a name of a time ",
         "scale by which data has been split AND aggregated by (could not ",
         "find breaks for that time scale in data's attributes)")
  }
  
  ## check breaks --------------------------------------------------------------
  
  found_breaks <- attrs$breaks[[ surv.scale ]]
  
  
  if (is.null(surv.breaks) && !is.null(found_breaks)) {
    surv.breaks <- found_breaks
  } else if (any(!surv.breaks %in% found_breaks)) {
    stop("given surv.breaks is not a subset of the breaks used to ",
         "split data; cannot proceed.")
  }
  surv.breaks <- sort(unique(surv.breaks))
  
  # data prep & subsetting -----------------------------------------------------
  subset <- substitute(subset)
  subset <- evalLogicalSubset(data, subset)
  
  origData <- data
  
  data <- data[subset, ]
  setDT(data)
  
  # handle count etc. variables ------------------------------------------------
  
  valVars <- c("d")
  valVars <- c(valVars, if (surv.method == "hazard") "pyrs" else c("n", "n.cens"))
  
  valVars <- c(valVars, if (surv.type == "surv.rel" && relsurv.method == "e2")  "d.exp" else NULL)
  
  valVars <- c(valVars, if (surv.type == "cif.rel")  "d.exp" else NULL)
  
  ppVars <- c("d.pp", "d.exp.pp", "d.pp.2", 
              if (surv.method == "hazard") "pyrs.pp" else c("n.cens.pp", "n.pp"))
  valVars <- c(valVars, if (surv.type == "surv.rel" && relsurv.method == "pp") ppVars else NULL)
  
  fo <- formals("survtab_ag")
  mc <- as.list(match.call())[-1]
  mc <- c(mc, fo[!names(fo) %in% names(mc)])
  
  mc <- mc[valVars]
  
  mc <- lapply(mc, function(elem) {
    evalPopArg(data = data, arg = elem, DT = TRUE, enclos = PF, recursive = TRUE)
    })
  
  ## NOTE: this does not delete but sets the value to NULL.
  mc[unlist(lapply(mc, function(x) {
    NROW(x) == 0L || is.null(x) || is.language(x) || inherits(x, "try-error")
    }))] <- NULL
  
  lackVars <- setdiff(valVars, names(mc[!unlist(lapply(mc, is.null))]))
  if (length(lackVars) > 0) {
    stop("Following arguments were NULL or could not be evaluated but are ",
         "required: ", paste0("'", lackVars, "'", collapse = ", "), ". ",
         "Usual suspects: arguments are NULL or refer to variables that ",
         "cannot be found in data.")
  }
  
  eventVars <- NULL
  ## NOTE: not sure if other arguments than 'd' should be allowed to be of 
  ## length > 1 (cause-specific 'd'); restricted for now to 'd' but easy to
  ## allow in the procedure below.
  ## nl will contain the names of the variables corresponding to each argument,
  ## e.g. d = c("d.1", "d.2"), etc.
  mc[[1]] <- data.table(mc[[1L]]) ## this avoids an exotic error in set().
  nl <- lapply(mc, names)
  for (k in 1:length(mc)) {
    jay <- argName <- names(mc[k])
    cn <- names(mc[[k]])
    
    if (length(cn) > 1) jay <- paste0(jay, ".", cn) ## e.g. d.1, d.2, ...
    if (argName %in% c("d")) {
      eventVars <- jay
      if (surv.type %in% c("surv.cause") && length(cn) == 1L) {
        stop("surv.type = 'surv.cause', but only one type of event supplied ",
             "via argument 'd'. If you want to compute cause-specific ",
             "survivals, please supply multiple types of events via ",
             "'d'; otherwise use surv.type = 'surv.obs'") 
      } else  if (length(cn) > 1 && !argName %in% c("d","d.pp", "d.pp.2", "n.pp")) {
        stop("'", argName, "' has/evaluates to ", length(cn), 
             " columns; only 'd', 'd.pp', and 'd'pp.2', 'n.pp' may evaluate ",
             "to more than one column of the value arguments")
      }
      
    }
    
    setnames(mc[[k]], cn, jay)
    set(mc[[1]], j = jay, value = mc[[k]])
    nl[[argName]] <- jay
  }
  mc <- mc[[1]]
  
  if (!is.null(eventVars)) {
    set(mc, j = "d", value = rowSums(mc[, mget(eventVars)]))
    valVars <- unique(c(valVars, "d", eventVars))
  }
  
  ## sum e.g. d.pp.1 + d.pp.2 = d.pp
  dna <- names(nl)[names(nl) %in% c("d.pp", "d.pp.2")]
  if (length(dna)) dna <- dna[unlist(lapply(nl[dna], function(x) length(x) > 1L))]
  if (length(dna)) {
    
    for (k in dna) {
      set(mc, j = k, value = mc[, rowSums(.SD), .SDcols = nl[[k]]])
    }
    setcolsnull(mc, unlist(nl[dna]))
    valVars <- setdiff(valVars, unlist(nl[dna]))
    valVars <- c(valVars, dna)
    valVars <- unique(valVars)
  }
  
  
  all_names_present(mc, valVars)
  setcolorder(mc, valVars)
  
  ## addition: internal weights use n at beginning of first interval
  
  if (is.character(weights)) {
    checkWeights(weights)
    if (!"n" %in% valVars) {
      n <- substitute(n)
      mc$n <- evalPopArg(n, data = data, enclos = PF)
      
      valVars <- unique(c(valVars, "n"))
      
      if (is.null(mc$n)) {
        
        stop("Requested internal weights to be computed and used to standardize ", 
             "estimates, but argument 'n' not supplied. This is currently ",
             "required for computing internal weights (the values of 'n' ", 
             "in the first interval will be used for this). Please supply 'n' ",
             "or supply hand-made weights (preferred for your clarity).")
      }
    } 
    
    data[, c("n") := mc$n]
    
  }
  
  
  # making weighted table of aggregated values ---------------------------------
  ## NOTE: at-risk counts require special treatment when surv.breaks
  ## are a subset of the available breaks: cannot sum at-risk figures!
  ## instead should simply pick the value at the start of the
  ## (now larger) interval. Will accomplish this by setting values not
  ## at the start of an interval to zero and summing anyway.
  if (surv.method == "lifetable") {
    wh_internal <- list(surv.breaks)
    names(wh_internal) <- surv.scale
    wh_internal <- data[wh_internal, on = eval(surv.scale), which = TRUE]
    wh_internal <- setdiff(1:nrow(data), wh_internal)
    mc[wh_internal, intersect(c("n", "n.pp"), names(mc)) := 0L]
  }
  
  ## NOTE: while ssSub will pass the whole column of e.g. fot values, which will
  ## not limit the data to e.g. up 5 years of follow-up if original data went 
  ## further, surv.breaks may be only up to 5 years and will limit the data
  ## in makeWeightsDT using a CJ-merge-trick appropriately (via custom.levels).
  bl <- list(surv.breaks)
  setattr(bl, "names", surv.scale)
  
  adjust <- evalPopArg(data, adjust, enclos = PF, naming = "model")
  
  iws <- if ("n" %in% names(data)) "n" else "pyrs"
  
  data <- makeWeightsDT(data = data, values = list(mc), enclos = PF,
                        print = NULL, formula = formula, adjust = adjust, 
                        by.other = surv.scale, Surv.response = FALSE,
                        custom.levels = bl, weights = weights,
                        internal.weights.values = iws,
                        custom.levels.cut.low = surv.scale)
  
  allVars <- attr(data, "makeWeightsDT")
  allVars[] <- lapply(allVars, function(x) if (length(x) == 0L) NULL else x)
  prVars <- allVars$prVars
  adVars <- allVars$adVars
  # boVars <- allVars$boVars ## this is surv.scale
  valVars <- allVars$vaVars
  
  ## to avoid e.g. 'factor(sex, 1:2)' going bonkers
  prVars_orig <- prVars
  if (length(prVars) > 0L) {
    prVars <- makeTempVarName(names = c(names(data), adVars), 
                              pre = paste0("print_", 1:length(prVars)))
  }
  adVars_orig <- adVars
  if (length(adVars) > 0L) {
    adVars <- makeTempVarName(names = c(names(data), prVars), 
                              pre = paste0("print_", 1:length(adVars)))
  }
  if (length(c(prVars, adVars))) setnames(data, c(prVars_orig, adVars_orig), c(prVars, adVars))
  byVars <- c(prVars, adVars)
  
  # formulate some needed variables --------------------------------------------
  setkeyv(data, c(byVars, surv.scale))
  data[, Tstop := surv.breaks[-1]]
  setnames(data, surv.scale, "Tstart")
  data[, delta := Tstop - Tstart]
  data[, surv.int := 1:.N, by = eval(byVars)]
  setcolorder(data, c(byVars, "surv.int", "Tstart", "Tstop", "delta", valVars, intersect(names(data), "weights")))
  
  if (surv.method == "lifetable") {
    testEvents <- data[, n - shift(n, n = 1, type = "lead", fill = NA), by = eval(byVars)]$V1
    testEvents <- data$n.cens + data$d - testEvents
    
    if (sum(abs(testEvents), na.rm = TRUE)) {
      on.exit({
        
        data[, "n.cens + d - (n-lead1_n)" := testEvents]
        wh <- testEvents != 0L
        wh <- wh & !is.na(wh)
        if (interactive()) {
          printSD <- c(byVars, "Tstop", "d", "n", "n.cens", 
                       "n.cens + d - (n-lead1_n)")
          print(data[wh, .SD, .SDcols = printSD], top = 5, nrow = 10)
          
        }
        
      }, add = TRUE)
      
      stop("Supplied n.cens and d do not sum to total number of events and ",
           "censorings based on n alone. Note that lifetable analysis ",
           "is currently not supported for period analysis (or other ",
           "comparable limitations of data).",
           if (interactive())" See table below and check your variables.")
    }
    rm(testEvents)
    data[, n.eff := n - n.cens/2L]
  }
  
  
  # compute observed survivals  ------------------------------------------------
  if (verbose) ostime <- proc.time()
  
  if (surv.method=="lifetable") {
    comp.st.surv.obs.lif(surv.table = data, surv.by.vars = byVars)
  }
  if (surv.method=="hazard") {
    comp.st.surv.obs.haz(surv.table = data, surv.by.vars = byVars)
  }
  
  data <- comp.st.conf.ints(data, al=1-conf.level, surv="surv.obs", transform = conf.type)
  
  if (verbose) cat("Time taken by computing observed survivals:", timetaken(ostime), "\n")
  
  
  ## empty surv.int checking ---------------------------------------------------
  testVar <- if (surv.method == "lifetable") "n" else "pyrs"
  ## sum over adjusting variables
  data <- test_empty_surv_ints(data, by = c(prVars, adVars), 
                               show.by = c(prVars_orig, adVars_orig),
                               sum.over = adVars,
                               test.var = testVar)
  
  ## sum over nothing
  if (length(adVars) > 0L) {
    data <- test_empty_surv_ints(data, by = c(prVars, adVars), 
                                 show.by = c(prVars_orig, adVars_orig),
                                 sum.over = NULL, test.var = testVar)
  }
  
  ## if adjusting, crop all estimates by adjusting variables
  ## to shortest estimate
  if (length(adVars)) {
    adLe <- data[, list(min = min(surv.int), max = max(surv.int)), keyby = eval(adVars)]
    adLe <- c(max(adLe$min), min(adLe$max))
    data <- data[surv.int %in% `:`(adLe[1L], adLe[2L])]
  }
  
  # create and print table of bad surv.ints ------------------------------------
  
  badObsSurv <- data$surv.obs == 0 | is.na(data$surv.obs)
  if (sum(badObsSurv)) {
    
    zerotab <- data[badObsSurv, 
                    list(first.bad.surv.int = min(as.integer(surv.int)), 
                         last.bad.surv.int = max(as.integer(surv.int)), 
                         surv.obs=min(surv.obs)), keyby = eval(byVars)]
    
    
    message("Some cumulative surv.obs were zero or NA:")
    if (length(byVars)) setnames(zerotab, c(prVars, adVars), c(prVars_orig, adVars_orig))
    print(zerotab)
    if (surv.method == "lifetable" && data[surv.obs == 0, .N] > 0) {
      message("NOTE: Zero surv.obs leads to zero relative survivals as well. Adjusting with weights WILL use the zero surv.obs / relative survival values.")
    }
    
  }
  rm(badObsSurv)
  
  # compute cause-specific survivals  ------------------------------------------
  if (surv.type == "surv.cause") {
    
    ## NOTE: these related to adjusting life-table estimates for delayed entry...
    #       data[, n.eff := n - n.cens/2 + n.de/2 + n.de.cens/4] # + d.de/2
    #       n.cens_1 := n.cens + (d-d_1)
    #       n.de.cens := n.de.cens + (d.de - d.de_1)
    
    if (surv.method == "lifetable") {
      for (k in eventVars) {
        k <- gsub(pattern = "d_", replacement = "", x = k)
        d_k <- paste0("d_", k)
        # d.de_k <- paste0("d.de_",k)
        
        n.eff_k <- paste0("n.eff_",k)
        
        ## old: " := n - (n.cens + (d-", d_k,")/2 + n.de/2 + (n.de.cens + d.de - ", d.de_k,")/4 )"
        # expr <- paste0(n.eff_k, " := n - (n.cens + (d-", d_k,")/2 )")
        
        set(data, j = c(n.eff_k), value = data$n.eff + (data$d - data[[d_k]])/2L ) # + d.de/2
        # data[,  eval(parse(text = expr), envir = .SD)]
        
      }
      
    }
    
    surv_names <- names(data)[grep("surv.obs", names(data))]
    surv_names <- c("d", if (surv.method == "lifetable") "n.eff" else NULL, surv_names)
    setnames(data, surv_names, paste0(surv_names, ".orig"))
    
    for (k in eventVars) {
      
      k <- gsub(pattern = "d.", replacement = "", x = k)
      setnames(data, paste0("d.",k), "d")
      
      if (surv.method=="hazard") {
        comp.st.surv.obs.haz(surv.table = data, surv.by.vars = byVars)
      } else {
        setnames(data, paste0("n.eff_", k), "n.eff")
        comp.st.surv.obs.lif(surv.table = data, surv.by.vars = byVars)
      }
      os.table <- comp.st.conf.ints(data, al=1-conf.level, surv="surv.obs", transform = conf.type)
      
      new_surv_names <- setdiff(surv_names, c("d", if (surv.method == "lifetable") "n.eff" else NULL))
      new_surv_names <- gsub("surv.obs", paste0("surv.obs.", k), new_surv_names)
      new_surv_names <- c(paste0(c("d.", if (surv.method == "lifetable") "n.eff." else NULL), k), new_surv_names)
      setnames(data, surv_names, new_surv_names)
      
      
    }
    setnames(data, paste0(surv_names, ".orig"), surv_names)
  }
  
  # compute cause-specifc/excess-case CIFs -------------------------------------
  if (surv.type %in% c("cif.obs", "cif.rel")) {
    
    data[, lag1_surv.obs := shift(surv.obs, n = 1L, type = "lag", fill = 1), by = eval(byVars)]
    data[, p.obs := surv.obs/lag1_surv.obs]
    
    if (surv.type == "cif.obs") {
      for (k in eventVars) {
        
        k <- gsub("d.", "", x = k)
        d.k <- paste0("d.", k)
        
        d.var <- paste0("d.",k)
        q.var <- paste0("q.", k)
        CIF_var <- paste0("CIF_", k)
        data[, (q.var)   := (1-p.obs)*get(d.var)/d]
        data[get(d.var) == 0L | d == 0L, (q.var) := 0]
        data[, (CIF_var) := cumsum(lag1_surv.obs*get(q.var)), by = eval(byVars)]
      }
    }
    
    if (surv.type == "cif.rel") {
      ## assuming d.exp in data
      data[, CIF.rel := (1-p.obs)*(d-d.exp)/d]
      data[d.exp>d, CIF.rel := NA]
      data[, CIF.rel := cumsum(lag1_surv.obs*CIF.rel), by = eval(byVars)]
    }
    
    ## SEs currently not known for CIFs; impute 0 to make adjusting work
    CIF_vars <- names(data)[substr(names(data),1,3) == "CIF"]
    data[, c(paste0("SE.", CIF_vars)) := 0L]
    
    setcolsnull(data, c("lag1_surv.obs", "p.obs", paste0("q.", substr(eventVars, 3, nchar(eventVars)))))
    
  }
  
  
  # relative survivals ---------------------------------------------------------
  if (surv.type == "surv.rel" & relsurv.method == "e2") {
    
    # compute r.e2 -------------------------------------------------------------
    comp.st.rs <- function(rs.table, rs.by.vars = byVars) {
      
      p.exp <- delta <- surv.exp <- surv.obs <- n.eff.pp <- 
        surv.obs <- NULL ## APPEASE R CMD CHECK
      ## EdererII
      
      ##-------------
      if (surv.method == "hazard") {
        rs.table[, p.exp := exp(-delta*d.exp/pyrs)] 
        rs.table[, surv.exp := cumprod(p.exp), by = eval(rs.by.vars)]
        comp.st.r.e2.haz(surv.table = rs.table, surv.by.vars = rs.by.vars)
      } else {
        rs.table[, p.exp := 1 - d.exp/n]
        rs.table[, surv.exp := cumprod(p.exp), by = eval(rs.by.vars)]
        
        if (rs.table[, min(surv.obs, na.rm=T) == 0]) {
          rs.table[surv.obs == 0, surv.exp := 1]
        }
        
        comp.st.r.e2.lif(surv.table = rs.table, surv.by.vars = rs.by.vars)
        
        if (rs.table[, min(surv.obs, na.rm=T) == 0]) {
          rs.table[surv.obs == 0, intersect(c("surv.exp","r.e2","SE.r.e2","r.e2.lo","r.e2.hi"), names(rs.table)) := 0]
        }
      }
      
      ## ------------
      
      rs.table <- comp.st.conf.ints(rs.table, al=1-conf.level, surv="r.e2", transform = conf.type)
      
      return(rs.table)
    }
    
    data <- comp.st.rs(rs.table = data)
    
    
  }
  
  # compute r.pp ---------------------------------------------------------------
  if (surv.type == "surv.rel" & relsurv.method == "pp") {
    
    all_names_present(data, c("d.pp", "d.exp.pp", "d.pp.2"))
    ## pohar perme: analysis weighted by expected cumulative survival
    comp.st.pp <- function(pp.table, by.vars = byVars) {
      ## relative survival
      if (surv.method == "hazard") {
        all_names_present(data, c("pyrs.pp"),
                          msg = paste0("internal error: work data did not have",
                                       " variable named pyrs.pp. Complain ",
                                       "to package maintainer if you see this."))
        comp.st.r.pp.haz(surv.table = pp.table, surv.by.vars = by.vars)
      } else {
        data[, n.eff.pp := n.pp - 0.5*n.cens.pp]
        all_names_present(data, c("n.pp", "n.cens.pp", "n.eff.pp"),
                          msg = paste0("internal error: work data did not have",
                                       " variable named n.eff.pp. Complain ",
                                       "to package maintainer if you see this."))
        comp.st.r.pp.lif(surv.table = pp.table, surv.by.vars = by.vars)
        
        if (pp.table[, min(surv.obs, na.rm=T) == 0]) {
          pp.table[surv.obs == 0, intersect(c("r.pp","SE.r.pp","r.pp.lo","r.pp.hi"), names(pp.table)) := 0]
        }
      }
      
      pp.table <- comp.st.conf.ints(pp.table, al=1-conf.level, surv="r.pp", transform = conf.type )
      
      return(pp.table)
    }
    data <- comp.st.pp(pp.table = data)
  }
  
  # compute adjusted estimates -------------------------------------------------
  if ("weights" %in% names(data)) {
    adEsts <- names(data)[substr(names(data), 1, 8) == "surv.obs"]
    adEsts <- c(adEsts, "r.e2", "r.pp")
    adEsts <- c(adEsts,  names(data)[substr(names(data),1,3)=="CIF"])
    adEsts <- intersect(adEsts, names(data))
    adEsts <- adEsts[unlist(lapply(adEsts, function(x) !substr(x, nchar(x)-2L, nchar(x)) %in% c(".lo", ".hi")))]
    adSEs <- paste0("SE.", adEsts)
    
    data.w <- data[, lapply(mget(c(adEsts, adSEs)), function(x) sum(x*weights)), keyby = c(prVars, "surv.int")]
    data <- data[, lapply(mget(valVars), sum), keyby = c(prVars, "surv.int", "Tstart", "Tstop", "delta")]
    data <- merge(data, data.w, by = c(prVars, "surv.int"), all = TRUE)
    setnames(data, c(adEsts, adSEs), paste0(c(adEsts, adSEs), ".as"))
    
    for (var in paste0(adEsts, ".as")) {
      data <- comp.st.conf.ints(data, al=1-conf.level, surv=var, transform =conf.type)
    }
    
    
  }
  
  # clean-up -------------------------------------------------------------------
  ## back to original names of print / adjust (used to avoid e.g. 
  ## 'factor(V1, 1:2)' going bonkers in data.table)
  if (length(c(prVars))) setnames(data, c(prVars), c(prVars_orig))
  prVars <- prVars_orig
  adVars <- adVars_orig
  
  ## reorder table, format numeric values, etc.
  
  miscVars <- intersect(names(data), c("surv.int", "Tstart", "Tstop", "delta"))
  
  survVars <- c("surv.obs.lo","surv.obs","surv.obs.hi","SE.surv.obs",
                "r.e2.lo","r.e2","r.e2.hi","SE.r.e2",
                "r.pp.lo","r.pp","r.pp.hi","SE.r.pp",
                paste0("CIF.rel.", c("lo", "", "hi")), "SE.CIF.rel",
                "surv.obs.as.lo","surv.obs.as","surv.obs.as.hi","SE.surv.obs.as",
                "r.e2.as.lo","r.e2.as","r.e2.as.hi","SE.r.e2.as",
                "r.pp.as.lo","r.pp.as","r.pp.as.hi","SE.r.pp.as",
                paste0("CIF.rel.as.", c("lo", "", "hi")), "SE.CIF.rel.as"
  )
  survVars <- intersect(survVars, names(data))
  
  ## which variables are estimates, SEs, CIs, etc.
  survVars.ca <- setdiff(names(data), c(prVars, valVars, miscVars, survVars))
  CIF_vars <- survVars.ca[substr(survVars.ca, 1,3)=="CIF" | substr(survVars.ca, 1,6)=="SE.CIF"]
  survVars <- c(survVars, CIF_vars)
  
  surv.obs.vars <- survVars.ca[substr(survVars.ca, 1,8) == "surv.obs" | substr(survVars.ca, 1,11) == "SE.surv.obs"]
  survVars <- c(survVars, surv.obs.vars)
  
  survVars <- unique(intersect(survVars, names(data)))
  
  ## remove some unuseful variables
  setcolsnull(data, c("SE.A", "SE.B"))
  setcolsnull(data, survVars[substr(survVars, 1, 6) == "SE.CIF"]) ## since they are zero for now
  survVars <- intersect(survVars, names(data))
  
  SEVars <- survVars[substr(survVars, 1, 3) == "SE."]
  CIVars <- survVars[substr(survVars, nchar(survVars) - 2L, nchar(survVars)) %in% c(".lo", ".hi")]
  estVars <- setdiff(survVars, c(SEVars, CIVars))
  
  order <- unique(c(prVars, miscVars, valVars, survVars))
  order <- intersect(order, names(data))
  
  setcolsnull(data, setdiff(names(data), order))
  setcolorder(data,order)
  
  setkeyv(data, c(prVars, "surv.int"))
  
  # attributes -----------------------------------------------------------------
  setkeyv(data, c(prVars, "surv.int"))
  setattr(data, "class", c("survtab", "data.table", "data.frame"))
  if (!getOption("popEpi.datatable")) setDFpe(data)
  if (length(prVars) == 0) prVars <- NULL ## might be character(0) 
  
  used_args$data <- origData
  used_args$formula <- formula
  used_args$weights <- evalRecursive(arg = weights, env = PF)$weights
  
  arglist <- list(call = this_call, 
                  arguments = used_args,
                  surv.breaks = surv.breaks,
                  print.vars = prVars,
                  adjust.vars = adVars,
                  value.vars = valVars,
                  misc.vars = miscVars,
                  surv.vars = survVars,
                  est.vars = estVars,
                  SE.vars = SEVars,
                  CI.vars = CIVars)
  varsArgs <- substr(names(arglist), nchar(names(arglist))-4L, nchar(names(arglist))) == ".vars"
  varsArgs <- names(arglist)[varsArgs]
  arglist[varsArgs] <- lapply(arglist[varsArgs], function(x) if (length(x) == 0L) NULL else x)
                  
  setattr(data, "survtab.meta", arglist)
  
  if (verbose) cat("Time taken by whole process: ", timetaken(starttime), "\n")
  data[]
}


# ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#               status = status %in% 1:2, pophaz = popmort, pp = TRUE,
#               aggre = list(sex, fot), fot = seq(0, 5, 1/12))
# ag[, d.exp := pmax(0L, from0to1 - 3L)]
# st <- survtab_ag(ag, surv.type = "surv.obs", surv.method = "hazard")
# st <- survtab_ag(ag, surv.type = "surv.cause", surv.method = "hazard", d = list(a = from0to1-3, b = 3))

# sire <- copy(sire)
# sire$sex <- rbinom(nrow(sire), size = 1, prob = 0.5)
# ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#               status = status %in% 1:2, pophaz = popmort, pp = TRUE,
#               aggre = list(sex, agegr = cut(dg_age, c(0,60,70,80, Inf), labels = FALSE), fot), 
#               fot = seq(0, 5, 1/12))
# ag <- lexpand(sire, birth = "bi_date", entry = "bi_date", exit = "ex_date",
#               status = status %in% 1:2,
#               aggre = list(sex, age), 
#               age = seq(0, 100, 1))
# wdt <- data.table(agegr = 1:4, weights = c(0.2, 0.4, 0.3, 0.1))
# wli <- list(agegr = c(0.2, 0.4, 0.3, 0.1))
# st <- survtab_ag(fot ~ sex + adjust(agegr), data = ag, surv.type = "surv.obs", surv.method = "hazard", weights = wli)
# st <- survtab_ag(fot ~ sex + adjust(agegr), data = ag, surv.type = "surv.rel", 
#                  d.pp = "from0to1.pp", d.pp.2 = "from0to1.pp.2", 
#                  d.exp.pp = "d.exp.pp", pyrs.pp = "ptime.pp",
#                  surv.method = "hazard", weights = wli,
#                  relsurv.method = "pp")
# ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
#               status = status, pophaz = popmort, pp = TRUE,
#               aggre = list(sex, agegr = cut(dg_age, c(0,60,70,80, Inf), labels = FALSE), fot), 
#               fot = seq(0, 5, 1/12))
# st <- survtab_ag(fot ~ sex + adjust(agegr), data = ag, 
#                  d = list(cand = from0to1, othd = from0to2),
#                  surv.type = "surv.cause", weights = wli)
# st <- survtab_ag(fot ~ sex, data = ag, surv.type = "surv.obs", surv.method = "hazard", adjust = "agegr", weights = wli)
# st <- survtab_ag(fot ~ adjust(agegr), data = ag, surv.type = "surv.obs", weights = wli)
# st <- survtab_ag(fot ~ 1, data = ag, adjust = "agegr", surv.type = "surv.obs", weights = wli)
# st <- survtab_ag(fot ~ 1, data = ag, adjust = "agegr", surv.type = "surv.obs", weights = wli)
# st <- survtab_ag(fot ~ 1, data = ag, surv.type = "surv.obs")

# wli2 <- wli
# wli$sex <- c(0.4, 0.6)
# st <- survtab_ag(fot ~ adjust(sex, agegr), data = ag, surv.type = "surv.obs", weights = wli)
# st <- survtab_ag(fot ~ adjust(agegr), data = ag, surv.type = "surv.obs", weights = wli["agegr"])
# ag[, d.exp := pmax(from0to1 - 1, 0L)]
# st <- survtab_ag(fot ~ adjust(sex, agegr), data = ag, surv.type = "surv.rel", weights = wli)
# st <- survtab_ag(fot ~ adjust(sex, agegr), data = ag, surv.type = "surv.cause", weights = wli)
# ag[, othd := pmax(from0to1 - 1L, 0L)]
# st <- survtab_ag(fot ~ adjust(sex, agegr), data = ag, d = list(cand = from0to1, othd = pmax(from0to1-1L, 0L)), surv.type = "surv.cause", weights = wli)
