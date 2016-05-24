##############################################################################
##
##   R package dynsurv by Xiaojing Wang, Jun Yan, and Ming-Hui Chen
##   Copyright (C) 2011
##
##   This file is part of the R package dynsurv.
##
##   The R package dynsurv is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package dynsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package dynsurv. If not, see <http://www.gnu.org/licenses/>.
##
##############################################################################

##############################################################################
# Utility functions
##############################################################################
# Expand row, require package plyr
expand <- function(data, id="id", time="time", status="status") {

    pos <- match(c(id, time, status), names(data))
    if (length(pos) != 3)
        stop("Variable names not match!\n")

    eventTime <- sort(unique(data[, time][data[, status] == 1]))

    foo <- function(x) {
        tStop <- union(subset(eventTime, eventTime <= max(x[, time])), max(x[, time]))
        tStart <- c(0, head(tStop, -1))
        st <- rep(0, length(tStop))
        st[tStop %in% x[x[, status] == 1, time]] <- 1

        cbind(x[1, -pos[2:3]], tStart=tStart, tStop=tStop, st=st, row.names=NULL)
    }

    res <- ddply(data, id, foo)
    names(res)[ncol(res)] <- status
    res
}

control_sfun <- function(df=5, knots=NULL, boundary=NULL) {
    list(df=df, knots=knots, boundary=boundary)
}

const <- function(x) {
    x
}

##############################################################################
# Fit a time-varying coefficient Cox model, using B-splines
##############################################################################
splineCox <- function(formula, data, control=list()) {

    Call <- match.call()
    control <- do.call("control_sfun", control)

    ########################################################
    # Model matrix after expanding the factor covariate
    mm <- model.matrix(formula, data)
    N <- nrow(mm)
    nBeta <- ncol(mm) - 1
    is.tv <- rep(TRUE, nBeta)
    is.tv[grep("const[(]([A-Za-z0-9._]*)[)]", colnames(mm)[-1])] <- FALSE

    cov.names <- gsub("const[(]([A-Za-z0-9._]*)[)]", "\\1", colnames(mm)[-1])

    ########################################################
    # Model frame before expanding the factor covariate
    mf <- model.frame(formula, data)
    nCov <- ncol(mf) - 1

    FtNms <- rep("Ft:", nCov)
    FtNms[grep("const[(]([A-Za-z0-9._]*)[)]", names(mf)[-1])] <- ""

    names(mf) <- gsub("const[(]([A-Za-z0-9._]*)[)]", "\\1", names(mf))

    # First 3 columns = c("id", "time", "status")
    DF <- cbind(id=1:N, mf[, 1][, 1:2], mf[, -1, drop=FALSE])

    ########################################################
    # Prepare B-spline paramters
    boundary <- control$boundary
    if (is.null(boundary))
        boundary <- range(DF$time)

    knots <- control$knots
    df <- control$df
    if (is.null(control$knots)) {
        insideTime <- subset(DF$time, DF$time >= boundary[1] & DF$time <= boundary[2])

        # number of interior knots = df - degree(3) - intercept(1)
        sq <- seq.int(from=0, to=1, length.out=df-2)[-c(1, df-2)]
        knots <- stats::quantile(insideTime, sq)
    }

    basis <- list(df=control$df, knots=knots, intercept=TRUE,
                  Boundary.knots=boundary)

    ########################################################
    # Call coxph to fit the expanded data
    newDF <- expand(DF, id="id", time="time", status="status")

    # B-spline basis matrix
    Ft <- do.call("bs", c(list(x=newDF$tStop), basis))

    newFml <- as.formula(paste("survival::Surv(tStart, tStop, status) ~ ",
                               paste(paste(FtNms, names(mf)[-1], sep=""),
                                     collapse="+"), "+ cluster(id)"))
    
    fit <- coxph(newFml, newDF)

    rl <- list(call=Call, control=control, bsp.basis=basis,
               N=N, nBeta=nBeta, cov.names=cov.names, is.tv=is.tv,
               coxph.fit=fit)

    class(rl) <- "splineCox"
    rl
}
