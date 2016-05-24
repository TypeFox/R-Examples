generateInitpts <- function (bestProfiledOut, vertices, ui, ci, hrepr, # hrep becomes useless 10/2015
                             fixedlist,precision, max.only = TRUE) {
    candidates <- list()
    trouble <- FALSE
    bary <- verticesBarycenter(vertices)
    parnbr <- length(bestProfiledOut)
    safety <- 1e-08
    if (precision == "rational") {
      bPcheck <- q2d(qmq(qmatmult(ui, matrix(d2q(bestProfiledOut))),
                         ci))
      if (all(bPcheck > 0)) {
        initlogl <- tofKpredict.nohull(bestProfiledOut, fixedlist = fixedlist)
      } else {
        trouble <- TRUE ## by default assumes a  potential problem when precision is double...
        initlogl <- NA
      }
    } else {
        bPcheck <- ui %*% bestProfiledOut - ci
        if (all(bPcheck > -safety/100)) {
            unsafe <- which(bPcheck < safety)
            nun <- length(unsafe)
            if (nun > 0) {
                wei <- max((rep(safety, nun) - bPcheck[unsafe])/((ui %*%
                  (bary - bestProfiledOut))[unsafe]))
                bestProfiledOut <- (1 - wei) * bestProfiledOut +
                  wei * bary
                rebPcheck <- ui %*% bestProfiledOut - ci
                if (all(rebPcheck > safety/10)) {
                  initlogl <- tofKpredict.nohull(bestProfiledOut,
                    fixedlist = fixedlist)
                }
                else {
                  initlogl <- NA
                }
            }
            else initlogl <- tofKpredict.nohull(bestProfiledOut,
                fixedlist = fixedlist)
        }
        else {
            trouble <- T
            initlogl <- NA
        }
    }
    if (!is.na(initlogl)) {
        tmplist <- list(list(par = bestProfiledOut, value = initlogl))
        names(tmplist) <- paste("input")
        candidates <- c(candidates, tmplist)
    }
    eps <- apply(vertices, 1, function(v) {
        1/min(1e+12, max(c(1e+06, abs(v - bary))))
    })
    bordel <- vertices + eps * t(bary - t(vertices))
    pred <- apply(bordel, 1, tofKpredict.nohull, fixedlist = fixedlist)
    if (any(is.na(pred)))
        stop.redef("(!) From generateInitpts(): Suspicious NA returned by tofKpredict.nohull.")
    snif <- bordel[which.max(pred), ]
    names(snif) <- colnames(bordel)
    tmplist <- list(list(par = snif, value = max(pred)))
    names(tmplist) <- paste("envelope")
    if (max.only && length(candidates) == 1) {
        if (tmplist[[1]]$value > candidates[[1]]$value)
            candidates <- tmplist
    }
    else candidates <- c(candidates, tmplist)
    if (TRUE) {
      if (trouble) {
        ptnbr <- floor(60^((parnbr/3)^(1/2)))
      } else {
        ptnbr <- floor(17^((parnbr/3)^(1/2)))
      }
      candid <- rvolTriangulation(ptnbr, volTriangulation(vertices)) ## 10/2015 appears faster
      blob <- apply(candid, 1, tofKpredict.nohull, fixedlist = fixedlist)
    } 
    if (length(na.omit(blob)) > 0) {
        snif <- candid[which.max(blob), ]
        names(snif) <- colnames(candid)
        tmplist <- list(rhull=list(par = snif, value = max(blob, na.rm = TRUE)))
        if (max.only && length(candidates) == 1) {
            if (tmplist[[1]]$value > candidates[[1]]$value)
                candidates <- tmplist ## list with on element '$rhull'
        }
        else candidates <- c(candidates, tmplist) ## list with several elements
    }
    return(candidates)
}
