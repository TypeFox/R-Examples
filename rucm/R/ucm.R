#'Unobserved components methods for a time series
#'@import KFAS
#'
#'
#'@description Function \code{ucm} decomposes a time series into components such as trend, seasonal, cycle, and the regression effects due to predictor series using Unobserved Components Model (UCM).
#'
#'@rdname ucm
#'@name ucm
#'@param formula an object of class \code{formula} containing the symbolic description of the model with dependent and independent terms. If there are no independent terms, replace rhs with 0. 
#'@param data a required data frame or list containing variables in the model.
#'@param irregular logical; if irregular component is to be included in the model. Defaults to \code{TRUE}.
#'@param irregular.var value to fix variance of irregular component.
#'@param level logical; if level is to be included in the model. Defaults to \code{TRUE}.
#'@param level.var value to fix variance of level component.
#'@param slope logical; if slope is to be included in the model along with level. Defaults to \code{FALSE}.
#'@param slope.var value to fix variance of the slope component.
#'@param season logical; if seasonal component is to be included in the model. Defaults to \code{FALSE}.
#'@param season.length value of length of seasonal component. Required when \code{season} is included.
#'@param season.var value to fix variance of seasonal component.
#'@param cycle logical; if cyclical component is to be included in the model. Defaults to \code{FALSE}.
#'@param cycle.period length of cyclical component. Required when \code{cycle} is included. 
#'@param cycle.var value to fix variance of cyclical component.
#'
#'@details Formula of the model can be of the forma as in \code{lm} with response variable on rhs and predictor variables or 0 (if no predictor variables) on the rhs. 
#'
#'@return object of class \code{ucm}, which is a list with the following components:
#'\item{est}{Estimates of predictor variables, if present.}
#'\item{irr.var}{Estimated variance of irregular component, if present.}
#'\item{est.var.level}{Estimated variance of the level component, if present.}
#'\item{est.var.slope}{Estimated variance of slope of the level, if present.}
#'\item{est.var.season}{Estimated variance of the seasonal component, if present.}
#'\item{est.var.cycle}{Estimated variance of the cyclical component, if present.}
#'\item{s.level}{An object of the same class as of dependent variable containing the time varying level values, if level is present.}
#'\item{s.lope}{An object of the same class as of dependent variable containing the time varying slope values, if slope is present.}
#'\item{s.season}{An object of the same class as of dependent variable containing the time varying seasonal values, if season is present.}
#'\item{s.cycle}{An object of the same class as of dependent variable containing the time varying cyclical values, if cycle is present.}
#'\item{vs.level}{A vector containing time varying estimated variance of level, if level is present.}
#'\item{vs.slope}{A vector containing time varying estimated variance of slope, if slope is present.}
#'\item{vs.season}{A vector containing time varying estimated variance of seasonal component, if season is present.}
#'\item{vs.cycle}{A vector containing time varying estimated variance of cyclical component, if cycle is present.}
#'\item{call}{Original call of the function.}
#'\item{model}{The original model of class \code{\link{SSModel}} from \code{\link{KFAS}} package.}
#'
#'@seealso \code{\link{KFAS}}, \code{\link{SSModel}} for a detailed discussion on State Space Models.
#'
#'@examples
#'modelNile <- ucm(Nile~0, data = Nile, slope = TRUE)
#'modelNile
#'modelNile$s.level
#'@export

ucm <- function (formula, data, irregular = TRUE, irregular.var = NA, 
    level = TRUE, level.var = NA, slope = FALSE, slope.var = NA, 
    season = FALSE, season.length = NA, season.var = NA, cycle = FALSE, 
    cycle.period = NA, cycle.var = NA) 
{
    if (missing(data)) {
      data <- environment(formula)
    }
    if (!(level) && slope) 
        stop("Level to be included to have slope")
    if (season && is.na(season.length)) 
        stop("Specify season length")
    if (cycle && is.na(cycle.period)) 
        stop("Specify cycle period")
    if (irregular) 
        H <- irregular.var
    if (level) {
        comp_trend <- deparse(substitute(SSMtrend(degree = degree, 
            Q = Q), list(degree = 1, Q = level.var)))
    }
    if (slope) {
        comp_trend <- deparse(substitute(SSMtrend(degree = degree, 
            Q = Q), list(degree = 2, Q = list(level.var, slope.var))))
    }
    if (season) {
        comp_sea <- deparse(substitute(SSMseasonal(period = period, 
            Q = Q), list(period = season.length, Q = season.var)))
    }
    if (cycle) {
        comp_cyc <- deparse(substitute(SSMcycle(period = period, 
            Q = Q), list(period = cycle.period, Q = cycle.var)))
    }
    all_terms <- terms(formula)
    indep.var <- attr(all_terms, "term.labels")
    dep.var <- all_terms[[2]]
    comp <- "0"
    if (length(indep.var) > 0) 
        comp <- paste(indep.var, collapse = "+")
    if (exists("comp_trend")) 
        comp <- paste(comp, comp_trend, sep = "+")
    if (exists("comp_sea")) 
        comp <- paste(comp, comp_sea, sep = "+")
    if (exists("comp_cyc")) 
        comp <- paste(comp, comp_cyc, sep = "+")
    ssm.formula <- paste0(dep.var, "~", comp)
    init.times <- 1
    if (irregular) 
        init.times <- init.times + 1
    if (level) 
        init.times <- init.times + 1
    if (slope) 
        init.times <- init.times + 1
    if (season) 
        init.times <- init.times + 1
    if (cycle) 
        init.times <- init.times + 2
    if (missing(data) || is.ts(data)) parm <- log(eval(call("var", dep.var)))
    else parm <- log(var(data[, as.character(dep.var)]))
    inits <- rep(parm, times = init.times)
    if (irregular) 
        modelH <- SSModel(as.formula(ssm.formula), H = H, data = data)
    else modelH <- SSModel(as.formula(ssm.formula), data = data)
    modelH <- fitSSM(inits = inits, model = modelH, method = "BFGS")$model
    out <- KFS(modelH, filtering = "state", smoothing = "state")
    irr.var <- modelH$H
    names(irr.var) <- "Irregular_Variance"
    if (length(indep.var) > 0) {
        find.indep.var <- match(indep.var, rownames(modelH$T))
        est <- out$alphahat[nrow(out$alphahat), find.indep.var]
        se <- sapply(find.indep.var, function(x) mean(sqrt(out$V[x, 
            x, ])))
    }
    else {
        est <- NULL
        se <- NULL
    }
    if (level) {
        find.level <- match("level", rownames(modelH$T))
        s.level <- out$alphahat[, find.level]
        vs.level <- out$V[find.level, find.level, ]
        est.var.level <- modelH$Q[find.level - length(indep.var), 
            find.level - length(indep.var), ]
        names(est.var.level) <- "Level_Variance"
    }
    else {
        s.level <- NULL
        vs.level <- NULL
        est.var.level <- NULL
    }
    if (slope) {
        find.slope <- match("slope", rownames(modelH$T))
        s.slope <- out$alphahat[, find.slope]
        vs.slope <- out$V[find.slope, find.slope, ]
        est.var.slope <- modelH$Q[find.slope - length(indep.var), 
            find.slope - length(indep.var), ]
        names(est.var.slope) <- "Slope_Variance"
    }
    else {
        s.slope <- NULL
        vs.slope <- NULL
        est.var.slope <- NULL
    }
    if (season) {
        find.season <- match("sea_dummy1", rownames(modelH$T))
        s.season <- out$alphahat[, find.season]
        vs.season <- out$V[find.season, find.season, ]
        est.var.season <- modelH$Q[find.season - length(indep.var), 
            find.season - length(indep.var), ]
        names(est.var.season) <- "Season_Variance"
    }
    else {
        s.season <- NULL
        vs.season <- NULL
        est.var.season <- NULL
    }
    if (cycle) {
        find.cycle <- match("cycle", rownames(modelH$T))
        s.cycle <- out$alphahat[, find.cycle]
        vs.cycle <- out$V[find.cycle, find.cycle, ]
        est.var.cycle <- modelH$Q[match(1, modelH$R[find.cycle, ,1]), match(1, modelH$R[find.cycle, ,1]), ]
        names(est.var.cycle) <- "Cycle_Variance"
    }
    else {
        s.cycle <- NULL
        vs.cycle <- NULL
        est.var.cycle <- NULL
    }
    if (is.ts(data)) {
        df <- length(data) - length(indep.var)
    }
    else df <- nrow(data) - length(indep.var)
    mc <- match.call(expand.dots = TRUE)
    res.out <- list(est = est, irr.var = irr.var, est.var.level = est.var.level, 
        est.var.slope = est.var.slope, est.var.season = est.var.season, 
        est.var.cycle = est.var.cycle, s.level = s.level, s.slope = s.slope, 
        s.season = s.season, s.cycle = s.cycle, vs.level = vs.level, 
        vs.slope = vs.slope, vs.season = vs.season, vs.cycle = vs.cycle)
    res.out$call <- mc
    res.out$model <- modelH
    class(res.out) <- "ucm"
    attr(res.out, "se") <- se
    attr(res.out, "df") <- df
    invisible(res.out)
}
