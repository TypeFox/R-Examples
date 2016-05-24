

# dput(brewer.pal(7, "Set1"))
.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
           "#FF7F00", "#FFF33", "#A65628")[-6]

                                        # dput(brewer.pal(8,"Paired"))
.paircols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
               "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")

figure_MatMult_i7 <- function(D) {
    if (missing(D)) D <- getBenchmarkData("i7_920")

    D <- D[ D$type=='matmult', -c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="Matrix Multiplication")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="Matrix Multiplication", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)
    invisible(NULL)
}

figure_MatMult_xeon <- function(D) {
    if (missing(D)) D <- getBenchmarkData("xeon_X5570")

    D <- D[ D$type=='matmult', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="Matrix Multiplication")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="Matrix Multiplication", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)
    invisible(NULL)
}

figure_QR_i7 <- function(D) {
    if (missing(D)) D <- getBenchmarkData("i7_920")

    D <- D[ D$type=='qr', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="QR Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="QR Decomposition", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)
    invisible(NULL)
}

figure_QR_xeon <- function(D) {
    if (missing(D)) D <- getBenchmarkData("xeon_X5570")

    D <- D[ D$type=='qr', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="QR Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="QR Decomposition", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)
    invisible(NULL)
}

figure_SVD_i7 <- function(D) {
    if (missing(D)) D <- getBenchmarkData("i7_920")

    D <- D[ D$type=='svd', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="SVD Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob","gpu")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="SVD Decomposition", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto","GPU"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)

    invisible(NULL)
}

figure_SVD_xeon <- function(D) {
    if (missing(D)) D <- getBenchmarkData("xeon_X5570")

    D <- D[ D$type=='svd', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="SVD Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="SVD Multiplication", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)

    invisible(NULL)
}

figure_LU_i7 <- function(D) {
    if (missing(D)) D <- getBenchmarkData("i7_920")

    D <- D[ D$type=='lu', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="LU Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="LU Multiplication", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)

    invisible(NULL)
}

figure_LU_xeon <- function(D) {
    if (missing(D)) D <- getBenchmarkData("xeon_X5570")

    D <- D[ D$type=='lu', ,-c(1:2,5)]
    op <- par(mfrow=c(1,2))
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension", ylab="Time in seconds", main="LU Decomposition")
    legend("topleft", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    matplot(x=D[,"nobs"], y=D[,c("ref","atlas","atl39","mkl","gotob")], type='l', lty=1, lwd=3, col=.cols,
            xlab="Matrix dimension (in logs)", ylab="Time in seconds (in logs)", main="LU Multiplication", log="xy")
    legend("bottomright", legend=c("Ref","Atlas","Atl39","MKL","Goto"), bty="n", col=.cols, lty=1, lwd=3)
    par(op)

    invisible(NULL)
}

figure_LogLogSlopes <- function() {
    LL <- loglogAnalysis()
    DF <- LL[["slope"]]

    DF[,"method"] <- ordered(as.character(DF[,"method"]),
                             levels=c("gpu", "goto", "mkl", "atl39", "atlas", "ref"))
    DF[,"type"] <- ordered(as.character(DF[,"type"]),
                           levels=c("matmult", "qr", "svd", "lu"))

    sb <- trellis.par.get("strip.background")
    sb[["col"]][1:2] <- c("gray80","gray90")
    trellis.par.set("strip.background", sb)

    ss <- trellis.par.get("superpose.symbol")
    ss[["col"]][1:6] <- .cols
    ss[["cex"]] <- rep(1.0, 7)
    ss[["pch"]] <- rep(19, 7)
    ss[["alpha"]] <- rep(0.75, 7)
    trellis.par.set("superpose.symbol", ss)

    with(DF, print(dotplot(method  ~ value | type, group=host,
                           xlab="Slope of elapsed times to matrix dimension in log/log chart",
                           key=simpleKey(text=c("i7","xeon"), space="bottom", column=2))))
    invisible(NULL)
}

figure_LogLogIntercept <- function() {
    LL <- loglogAnalysis()
    DF <- LL[["intercept"]]

    DF[,"method"] <- ordered(as.character(DF[,"method"]),
                             levels=c("gpu", "goto", "mkl", "atl39", "atlas", "ref"))
    DF[,"type"] <- ordered(as.character(DF[,"type"]),
                           levels=c("matmult", "qr", "svd", "lu"))

    sb <- trellis.par.get("strip.background")
    sb[["col"]][1:2] <- c("gray80","gray90")
    trellis.par.set("strip.background", sb)

    ss <- trellis.par.get("superpose.symbol")
    ss[["col"]][1:6] <- .cols
    ss[["cex"]] <- rep(1.0, 7)
    ss[["pch"]] <- rep(19, 7)
    ss[["alpha"]] <- rep(0.75, 7)
    trellis.par.set("superpose.symbol", ss)

    with(DF, print(dotplot(method  ~ value | type, group=host,
                           xlab="Intercept of elapsed times to matrix dimension in log/log chart",
                           key=simpleKey(text=c("i7","xeon"), space="bottom", column=2))))
    invisible(NULL)
}

figure_LogLogLattice <- function(titles=TRUE) {
    D <- rbind(getBenchmarkData("i7_920"),
               getBenchmarkData("xeon_X5570"))

    DM <- melt(D, id.vars=c("host", "type", "datum", "nobs", "nrun"))
    DM[,"type"] <- ordered(as.character(DM[,"type"]),
                           levels=c("matmult", "qr", "svd", "lu"))

    DM[,"host"] <- ordered(DM[,"host"], levels=c("i7_920", "xeon_X5570"))
    levels(DM[,"host"]) <- c("i7", "xeon")

    sb <- trellis.par.get("strip.background")
    sb[["col"]][1:2] <- c("gray80","gray90")
    trellis.par.set("strip.background", sb)

    sl <- trellis.par.get("superpose.line")
    sl[["col"]] <- .cols
    trellis.par.set("superpose.line", sl)

    op <- options(scipen=5)
    with(DM,print(xyplot(value ~ nobs| type+host,
                         group=variable, lwd=2,
                         scales=list(x=list(log=TRUE,at=c(100,400,1500,5000),labels=c(100,400,1500,5000)),
                                     y=list(log=TRUE,at=c(0.0001,0.01,1,100),labels=c(0.0001,0.01,1,100))),
                         panel=function(x,subscripts,groups,...) {
                             panel.superpose(x,subscripts,groups,type='l',...)
                         },
                         key=simpleKey(text=c("ref","atlas","atl93","mkl","goto","gpu"),
                                       space="right", lines=TRUE, points=FALSE),
                         xlab="Matrix dimension (on logarithmic axis)",
                         ylab="Elapsed time in seconds (on logarithmic axis)",
                         main=ifelse(titles,paste("Benchmarking BLAS and GPU:",
                                    "Comparing six implementations on four methods across two architectures"),
                              ""),
                         sub=ifelse(titles,paste("Benchmark setup, code, data and analysis are available in the R package",
                                   "gcbd (Eddelbuettel, 2010) via every CRAN mirror"),
                              "")
                         )))
    options(op)
    invisible(NULL)
}

figure_Lattice <- function(titles=TRUE) {
    D <- rbind(getBenchmarkData("i7_920"),
               getBenchmarkData("xeon_X5570"))

    DM <- melt(D, id.vars=c("host", "type", "datum", "nobs", "nrun"))
    DM[,"type"] <- ordered(as.character(DM[,"type"]),
                           levels=c("matmult", "qr", "svd", "lu"))

    DM[,"host"] <- ordered(DM[,"host"], levels=c("i7_920", "xeon_X5570"))
    levels(DM[,"host"]) <- c("i7", "xeon")

    sb <- trellis.par.get("strip.background")
    sb[["col"]][1:2] <- c("gray80","gray90")
    trellis.par.set("strip.background", sb)

    sl <- trellis.par.get("superpose.line")
    sl[["col"]] <- .cols
    trellis.par.set("superpose.line", sl)

    op <- options(scipen=5)
    with(DM,print(xyplot(value ~ nobs| type+host,
                         group=variable, lwd=2,
                         ylim=c(0,30),
                         #scales=list(x=list(log=TRUE,at=c(100,400,1500,5000),labels=c(100,400,1500,5000)),
                         #            y=list(log=TRUE,at=c(0.0001,0.01,1,100),labels=c(0.0001,0.01,1,100))),
                         panel=function(x,subscripts,groups,...) {
                             panel.superpose(x,subscripts,groups,type='l',...)
                         },
                         key=simpleKey(text=c("ref","atlas","atl93","mkl","goto","gpu"),
                                       space="right", lines=TRUE, points=FALSE),
                         xlab="Matrix dimension",
                         ylab="Elapsed time in seconds (capped at 30 seconds)",
                         main=ifelse(titles,paste("Benchmarking BLAS and GPU: Comparing six",
                                     "implementations on four methods across two architectures"), ""),
                         sub=ifelse(titles,paste("Benchmark setup, code, data and analysis are available",
                                   "in the R package gcbd (Eddelbuettel, 2010) via every CRAN mirror"), "")
                         )))
    options(op)
    invisible(NULL)
}

figure_LatticeByArch <- function(titles=TRUE) {
    D <- rbind(getBenchmarkData("i7_920"),
               getBenchmarkData("xeon_X5570"))

    DM <- melt(D, id.vars=c("host", "type", "datum", "nobs", "nrun"))
    DM[,"type"] <- ordered(as.character(DM[,"type"]),
                           levels=c("matmult", "qr", "svd", "lu"))

    DM[,"host"] <- ordered(DM[,"host"], levels=c("i7_920", "xeon_X5570"))
    levels(DM[,"host"]) <- c("i7", "xeon")

    sb <- trellis.par.get("strip.background")
    sb[["col"]][1:2] <- c("gray80","gray90")
    trellis.par.set("strip.background", sb)

    sl <- trellis.par.get("superpose.line")
    sl[["col"]] <- .cols
    trellis.par.set("superpose.line", sl)

    op <- options(scipen=5)
    with(DM,print(xyplot(value ~ nobs | type+variable,
                         group=host, lwd=2,
                         ylim=c(0,30),
                         panel=function(x,subscripts,groups,...) {
                             panel.superpose(x,subscripts,groups,type='l',...)
                         },
                         xlab="Matrix dimension",
                         ylab="Elapsed time in seconds",
                         auto.key=TRUE,
                         main=ifelse(titles,paste("Benchmarking BLAS and GPU:",
                                    "Comparing six implementations on four methods across two architectures"),
                                    ""),
                         sub=ifelse(titles,paste("Benchmark setup, code, data and analysis are available in the R package",
                                   "gcbd (Eddelbuettel, 2010) via every CRAN mirror"),
                                    "")
                         )))
    options(op)
    invisible(NULL)
}

