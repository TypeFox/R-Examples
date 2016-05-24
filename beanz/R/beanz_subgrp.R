

##-----------------------------------------------------------------------------------
##                           SUBGROUP ANALYSIS
##-----------------------------------------------------------------------------------

#' Get subgroup treatment effect estimation and variance
#'
#' Compute subgroup treatment effect estimation and variance from subject level data.
#'
#' @param data.all subject level dataset
#'
#' @param var.resp column name in \code{data.all} for response
#'
#' @param var.trt column name  in \code{data.all} for treatment assignment
#'
#' @param var.cov array of column names in \code{dat.all} that corresponds to binary or
#'     ordinal baseline covaraites
#'
#' @param var.censor column name in \code{data.all} for censoring if the
#'     response is time to event data
#'
#' @param resptype type of response. The options are \code{binary},
#'     \code{continuous} or \code{survial}
#'
#' @return A dataframe with treatment effect estimation and variance for each subgroup
#'
#' @importFrom survival coxph Surv
#'
#'
#' @examples
#'
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' subgrp.effect <- r.get.subgrp.raw(solvd.sub,
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);}
#'
#' @export
#'

r.get.subgrp.raw <- function(data.all, var.resp, var.trt, var.cov, var.censor,
                             resptype=c("continuous", "binary", "survival")) {

    resptype <- match.arg(resptype);

    ss     <- r.get.subgroup(data.all, var.cov);
    subgrp <- ss$subgrp;
    xgrp   <- ss$xgrp;

    y      <- data.all[,var.resp];
    trt    <- data.all[,var.trt];

    if (!is.null(var.censor)) {
        censor <- data.all[, var.censor];
    }

    ##definition
    ng  <- nrow(subgrp);
    rst <- 1:ng;
    rst <- cbind(rst, subgrp);

    y.stat <- NULL;
    for (k in 1:ng) {
        inx     <- which(k == xgrp);
        cur.y   <- y[inx];
        cur.trt <- trt[inx];

        if (resptype == "survival") {
            ##survival
            cur.censor <- censor[inx]
            cur.mdl    <- coxph(Surv(cur.y, cur.censor) ~ cur.trt);
            eff        <- coef(summary(cur.mdl))[1];
            eff.sd     <- coef(summary(cur.mdl))[3];
        } else {
            ##continuous or binary
            if (resptype == "continuous") {
                cur.fm <- gaussian();
            } else if (resptype == "binary") {
                cur.fm <- binomial();
            }

            cur.mdl <- glm(cur.y~cur.trt, family=cur.fm);
            eff     <- coef(summary(cur.mdl))[2, 1];
            eff.sd  <- coef(summary(cur.mdl))[2, 2];
        }

        ##add result
        y.stat  <- rbind(y.stat, c(eff, eff.sd^2, length(inx)));
    }

    rst           <- cbind(rst, y.stat);
    colnames(rst) <- c("Subgroup", colnames(subgrp), "Estimate", "Variance",  "N");
    rownames(rst) <- NULL;

    ##return
    rst <- as.data.frame(rst);
}

#' Get subgroup treatment effect estimation and variance
#'
#' Compute subgroup treatment effect estimation and variance for subgroup effect
#' summary data. The estimation and variance are combined if there are multiple
#' record of the same subgroup, defined by the covariates, in the data.
#'
#' @param data.all subject level dataset
#'
#' @param var.ey column name in \code{data.all} for estimated treatment effect
#'
#' @param var.variance column name in \code{data.all} for variance of subgroup
#'     treatment assignment
#'
#' @param var.cov array of column names in \code{dat.all} that corresponds to
#'     binary or ordinal baseline covaraites
#'
#' @return A dataframe with treatment effect estimation and variance for each subgroup
#'
#' @export
r.subgrp.effect <- function(data.all, var.ey, var.variance, var.cov) {
    ey       <- data.all[,var.ey];
    variance <- data.all[,var.variance];

    ss     <- r.get.subgroup(data.all, var.cov);
    subgrp <- ss$subgrp;
    xgrp   <- ss$xgrp;
    ng     <- nrow(subgrp);
    y.stat <- NULL;
    for (k in 1:ng) {
        inx     <- which(k == xgrp);
        ##cur.y   <- mean(ey[inx]);
        ##cur.var <- sum(variance[inx])/length(inx)^2;
        comb    <- get.merge.subg(ey[inx], variance[inx]);
        y.stat  <- rbind(y.stat, comb);
    }

    rst           <- cbind(1:ng, subgrp, y.stat);
    colnames(rst) <- c("Subgroup", colnames(subgrp), "Estimate", "Variance");
    rownames(rst) <- NULL;

    ##return
    rst <- as.data.frame(rst);
}


##get subgroup of each subject and the definitions
r.get.subgroup <- function(data.all, var.cov) {
    covs  <- data.all[,var.cov, drop=FALSE];
    u.cov <- unique(covs);
    u.cov <- u.cov[do.call(order, u.cov),,drop=FALSE];

    ##u.cov <- get.s(covs, 1 == length(var.cov));
    ##colnames(u.cov) <- var.cov;
    ng    <- nrow(u.cov);
    xgrp  <- rep(-1, nrow(data.all));
    for (i in 1:ng) {
        cur.inx <- which(apply(covs, 1,
                               function(x) {
                                   all(x == u.cov[i,])
                               }));
        xgrp[cur.inx] <- i;
    }

    list(subgrp=u.cov, xgrp=xgrp);
}


##merge subgroups
get.merge.subg <- function(estimates, variances) {
    ## cochran anova estimates for inter-group variances
    tau2 <- max(0, var(estimates)-mean(variances), na.rm=TRUE);

    wi <- 1/(variances + tau2);

    com.mean <- sum(wi * estimates)/sum(wi);
    com.var  <- 1/sum(wi);

    c(com.mean, com.var);
}


##get subgroups
get.s <- function(covs, is.vec=FALSE) {
    u.cov <- unique(covs);
    if (is.vec) {
        u.cov      <- sort(u.cov);
        dim(u.cov) <- c(length(u.cov), 1);
    } else {
        u.cov <- u.cov[do.call(order, as.data.frame(u.cov)),];
    }
    u.cov
}
