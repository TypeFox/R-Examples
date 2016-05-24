###############################################################################
## package 'secrdesign'
## select.stats.R
## 2014-11-25 moved from methods.R
## 2015-02-17 user-specified 'true' values
## 2016-02-04 fix weighted problem for param != 'D'
###############################################################################

weighted <- function (onescenario, param) {
  ## 2016-02-04 rename argument scenario to onescenario to ensure OK for param != 'D'
    with(onescenario, {
        if (param == 'D')
            sum(D)
        else if (param == 'pmix') {
            warning("assuming first scenario row corresponds to requested pmix")
            D[1]/sum(D)
        }
        else {
            wt <- D/sum(D)
            sum(onescenario[,param] * wt)
        }
    })
}

select.stats <- function (object, parameter = 'D', statistics, true) {
    
    if (!inherits(object, 'estimatetables'))
        stop ("select.stats requires input of class estimatetables")
    if (is.na(object$outputtype))
        stop ("cannot select.stats output of unknown type")
    if (object$outputtype == 'secrfit')
        stop ("cannot select.stats secr fit - use predict() first")

    #if (missing(estname) | missing(SEname)) {
    estname <- ""
    SEname <- ""
    if (object$outputtype %in% c('predicted', 'derived', 'regionN')) {
        estname <- 'estimate'
        SEname <- 'SE.estimate'
    }
    else if (object$outputtype == 'coef') {
        estname <- 'beta'
        SEname <- 'SE.beta'
    }
    #}
    
    typical <- object$output[[1]][[1]]  ## first scenario, first replicate
    stat0 <- names(typical)[sapply(typical, is.numeric)]
    if (missing(statistics)) {
        stat1 <- stat0
        stat2 <- c('RB','RSE','COV')
        ## stat2 <- c('true','RB','RSE','COV','ERR')
    }
    else {
        stat1 <- statistics[statistics %in% stat0]
        stat2 <- statistics[statistics %in% c('true','RB','RSE','COV','ERR')]
    }
    if (any(stat2 %in% c('true','RB','RSE','COV','ERR')) & (estname == '')) {
        stat2 <- character(0)
        warning ("cannot compute requested statistics with your data")
    }

    extractfn <- function (out, true, estimated) {
        getij <- function(df, i, j) {
            if (nrow(df) == 0)
                rep(NA, length(j))
            else
                df[i,j]
        }
        if (length(stat1)>0) {
            tmp <- lapply(out, getij, parameter, stat1)
            tmp <- do.call (rbind, tmp)
            rownames(tmp) <- 1:nrow(tmp)
        }
        else
            tmp <- matrix(nrow = length(out), ncol = 0)

        for (st in stat2) {
            if (st == 'true') {
                tmp <- cbind(tmp, rep(true, nrow(tmp)))
            }
            if (st == 'RB') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, (est - true) / true)
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'RSE') {
                est <- sapply(out, getij, parameter, estname)
                SE.est <- sapply(out, getij, parameter, SEname)
                tmp <- cbind(tmp, SE.est/est)
            }
            if (st == 'ERR') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, abs(est-true))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'COV') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    lcl <- sapply(out, getij, parameter, 'lcl')
                    ucl <- sapply(out, getij, parameter, 'ucl')
                    tmp <- cbind(tmp, as.numeric((true>lcl) & (true<ucl)))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
        }
        colnames(tmp) <- c(stat1, stat2)
        tmp
    }

    uniqueScenarioIndex <- match(unique(object$scenarios$scenario), object$scenarios$scenario)
    fitIndex <- object$scenarios$fitindex[uniqueScenarioIndex]
    estimated <- sapply(fitIndex,
                        function(x) {
                            method <- if ('method' %in% names(object$fit.args))
                                object$fit.args$method
                            else
                                object$fit.args[[x]]$method
                            no <- object$fit & (method != 'none')
                            ifelse(length(no) == 0, TRUE, no)
                        }
                        )
    if (missing(true)) {
    splitScenarios <- split(object$scenarios, object$scenarios$scenario)
    trueD <- sapply(splitScenarios, weighted, param='D')   ## vector length = number of scenarios

    if (object$outputtype == 'regionN') {
        true <- trueD
        true <- true * attr(object, 'regionsize')  ## for each scenario
    }
    else if (object$outputtype %in% c('predicted','derived')){
        if (parameter == 'D')
            true <- trueD
        else
            ## true <- object$scenarios[,parameter]
            true <- sapply(splitScenarios, weighted, parameter)
    }
    else true <- NA
}
else if (length(true) != length(object$output))
    stop ("specify one 'true' value for each scenario")
    object$output <- mapply(extractfn,
                            object$output,
                            true,
                            estimated,
                            SIMPLIFY = FALSE)
    object$outputtype <- 'numeric'
    class(object) <- c('selectedstatistics', 'secrdesign', 'list')
    attr(object, 'parameter') <- parameter
    object
}
