
p3state <- function (data, coxdata = NULL, formula = NULL, regression = NULL)
{
    require(survival)
    if (missing(data))
        stop("Argument 'data' is missing with no default")
    if (!is.data.frame(data))
        stop("Argument 'data' must be a data.frame")
    if (any(names(data)[1:5] != c("times1", "delta", "times2",
        "time", "status")))
        stop("'data' must contain the right variables")
    if (any(data[, 2] != 0 & data[, 2] != 1))
        stop("The variable 'delta' in the argument 'data' must be 0 or 1")
    if (any(data[, 5] != 0 & data[, 5] != 1))
        stop("The variable 'delta' in the argument 'data' must be 0 or 1")
    if (any(is.na(data[, 1])))
        stop("Missing values are not allowed in variable 'times1'")
    if (any(is.na(data[, 3])))
        stop("Missing values are not allowed in variable 'times2'")
    if (mode(data[, 1]) != "numeric")
        stop("The variable 'times1' must be numeric")
    if (mode(data[, 3]) != "numeric")
        stop("The variable 'times2' must be numeric")
    if (!is.numeric(data[, 1]))
        stop("The variable 'times1' must be numeric")
    if (!is.numeric(data[, 3]))
        stop("The variable 'times2' must be numeric")
    if (any(data[, 1] + data[, 3] != data[, 4]))
        stop("The variable 'time' in the argument 'data' must be equal to times1+times2")
    if (any(data[, 2] == 0 & data[, 3] > 0))
        stop("The variable 'times2' in the argument 'data' must be equal to 0 when delta=0")
    if (any(data[, c(1, 3)] < 0))
        stop("The variables 'times1' and 'times2' in the argument 'data' cannot be negative")
    if (sum(data[, 2]) == 0)
        stop("This is not a multi-state model")
    if (missing(formula) & missing(regression))
        regression <- FALSE
    else regression <- TRUE
    if (missing(formula) & (regression == TRUE))
        stop("Argument 'formula' is missing with no default")
    estimate <- TRUE
    if (missing(coxdata) & regression == TRUE)
        coxdata <- data.creation.reg(data)
    nt12 <- sum(data$delta == 1)
    nt13 <- sum(data$delta == 0 & data$status == 1)
    nt11 <- nrow(data) - nt12 - nt13
    nt23 <- sum(data$delta == 1 & data$status == 1)
    nt22 <- sum(data$delta == 1 & data$status == 0)
    descriptives <- c(nt12, nt13, nt11, nt23, nt22)
    if (estimate == TRUE) {
        mydata <- data
        m <- nrow(mydata)
        mat <- matrix(data = NA, ncol = 7, nrow = m)
        mat[, 1] <- mydata[, 1]
        mat[, 2] <- mydata[, 2]
        mat[, 3] <- mydata[, 3]
        mat[, 4] <- mydata[, 4]
        mat[, 5] <- mydata[, 5]
        fit <- survfit(Surv(mat[, 4], mat[, 5]) ~ 1)
        fit2 <- survfit(Surv(mat[, 1], mat[, 2] + (1 - mat[,
            2]) * mat[, 5]) ~ 1)
        mat <- rbind(c(0, 0, 0, 0, 0, 0, 1), mat)
        unctimes <- sort(unique(mydata[mydata[,5]==1,4]))
        KMW <- vector(length = length(unctimes))
        for (k in 1:length(summary(fit)$surv)) KMW[k] <- (c(1,
            summary(fit)$surv)[k] - summary(fit)$surv[k])/summary(fit)$n.event[k]
        for (k in 1:length(summary(fit)$surv)) {
            p <- which(mat[, 4] == unctimes[k] & mat[, 5] ==
                1)
            mat[p, 6] <- KMW[k]
        }
        mat2 <- matrix(data = NA, ncol = 7, nrow = (m + 1))
        mat2[, 1] <- sort(mat[, 1])
        mat2[, 2] <- mat[order(mat[, 1]), 2]
        mat2[, 3] <- mat[order(mat[, 1]), 3]
        mat2[, 4] <- mat[order(mat[, 1]), 4]
        mat2[, 5] <- mat[order(mat[, 1]), 5]
        mat2[, 6] <- mat[order(mat[, 1]), 6]
        mat2[, 7] <- mat[, 7]
        mat2[1, 7] <- 1
        for (k in 1:(length(summary(fit2)$time))) {
            p <- which(mat2[, 1] == summary(fit2)$time[k])
            mat2[p, 7] <- summary(fit2)$surv[k]
        }
        q <- which(is.na(mat2[, 7]))
        lq <- length(q)
        for (kk in 1:lq) mat2[q[kk], 7] <- mat2[q[kk] - 1, 7]
        mydata <- data.frame(mat2)
        colnames(mydata) <- c("times1", "delta", "times2", "time",
            "status", "Wi", "H")
        toplot <- mydata
    }
    if (regression == TRUE) {
        fmla <- attr(terms(formula), "term.labels")
        ncov <- length(fmla)
        predictor <- fmla
        colvar <- rep(0, ncov)
        for (k in 1:ncov) {
            if (fmla[k] %in% names(data))
                colvar[k] <- which(names(data) == fmla[k])
            else {
                for (j in 1:ncol(data)) {
                  if (any(grep(names(data)[j], fmla[k])))
                    colvar[k] <- j
                }
            }
        }
        if (any(colvar == 0))
            stop("'formula' must contain the right variables")
        fmla2 <- c(fmla, "treat")
        fmla3 <- c(fmla, "start")
        covar <- as.formula(paste(" Surv(coxdata[,2],coxdata[,3],coxdata[,4])~ ",
            paste(fmla, collapse = "+")))
        covar2 <- as.formula(paste(" Surv(coxdata$start,coxdata$stop,coxdata$event)~ ",
            paste(fmla2, collapse = "+")))
        covar3 <- as.formula(paste(" Surv(coxdata$stop,coxdata$event)~ ",
            paste(fmla3, collapse = "+")))
        tdcm <- coxph(covar2, data = coxdata)
        if (nt13 > 0) {
            cmm13 <- coxph(covar, data = coxdata, subset = (coxdata[,
                5] == 0), na.action = na.exclude)
        }
        covar <- as.formula(paste(" Surv(coxdata[,2],coxdata[,3],1-coxdata[,4]-coxdata[,6])~ ",
            paste(fmla, collapse = "+")))
        cmm12 <- coxph(covar, data = coxdata, subset = (coxdata[,
            5] == 0), na.action = na.exclude)
        covar <- as.formula(paste(" Surv(coxdata[,2],coxdata[,3],coxdata[,4])~ ",
            paste(fmla, collapse = "+")))
        cmm23 <- coxph(covar, data = coxdata, subset = (coxdata[,
            5] == 1), na.action = na.exclude)
        covar <- as.formula(paste(" Surv(coxdata[,3]-coxdata[,2],coxdata[,4])~ ",
            paste(fmla, collapse = "+")))
        TMA <- coxph(Surv(stop, event) ~ start, data = coxdata,
            subset = (coxdata[, 5] == 1))
        TMA2 <- coxph(covar3, data = coxdata, subset = (coxdata[,
            5] == 1))
        csmm23 <- coxph(covar, data = coxdata, subset = (coxdata[,
            5] == 1), na.action = na.exclude)
    }
    if (estimate == TRUE & regression == FALSE)
        object <- list(descriptives = descriptives, datafr = toplot)
    if (estimate == FALSE & regression == TRUE & nt13 > 0)
        object <- list(descriptives = descriptives, tdcm = tdcm,
            msm13 = cmm13, msm12 = cmm12, cmm23 = cmm23, csmm23 = csmm23,
            tma = TMA, tma2 = TMA2)
    if (estimate == FALSE & regression == TRUE & nt13 == 0)
        object <- list(descriptives = descriptives, tdcm = tdcm,
            msm12 = cmm12, cmm23 = cmm23, csmm23 = csmm23, tma = TMA,
            tma2 = TMA2)
    if (estimate == TRUE & regression == TRUE & nt13 > 0)
        object <- list(descriptives = descriptives, datafr = toplot,
            tdcm = tdcm, msm13 = cmm13, msm12 = cmm12, cmm23 = cmm23,
            csmm23 = csmm23, tma = TMA, tma2 = TMA2)
    if (estimate == TRUE & regression == TRUE & nt13 == 0)
        object <- list(descriptives = descriptives, datafr = toplot,
            tdcm = tdcm, msm12 = cmm12, cmm23 = cmm23, csmm23 = csmm23,
            tma = TMA, tma2 = TMA2)
    class(object) <- "p3state"
    object
}
