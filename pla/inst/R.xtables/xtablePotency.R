xtablePotency <- function(Fits, digits = 4, showRelative = FALSE) {

    K <- Fits@pheur$K
    KP <- Fits@pheur$KP
    if (is.null(dim(K)))
        nT <- 1
    else
        nT <- dim(K)[[1]]
    if (is.null(dim(KP)))
        nF <- 1
    else {
        if ((dim(KP)[[1]] %% nT) == 0)
            nF <- dim(KP)[[1]] / nT
        else
            nF <- 1
    }
    if ((nT == 1) & (nF == 1))
        X <- rbind(K, c(NA, KP[-c(1, 5)], NA))
    else
        X <- rbind(K, cbind(NA, KP[, -c(1, 5)], NA))
    ## warning(paste0("#Treatments: ", nT, "  #Factors: ", nF))
    if (nT == 1)
        if (nF == 1)
            dimnames(X)[[1]] <- c("Log-transformed", "")
        else
            dimnames(X)[[1]] <- c("Log-transformed", dimnames(Fits@pheur$KP)[[1]])
    else
        if (nF == 1) {
            A <- c("Log-transformed", "")
            B <- dimnames(K)[[1]]
            dimnames(X)[[1]] <- paste(rep(A, rep(length(B), 2)), rep(B, 2))
        } else { }
    ## warning(dimnames(X)[[2]])
    ## warning(names(Fits@pheur$KP))
    ## warning(dimnames(Fits@pheur$KP)[[2]])
    nmsX <- dimnames(X)[[2]]
    if ((nT == 1) & (nF == 1))
        nmsX <- names(Fits@pheur$KP)
    else
        nmsX <- dimnames(Fits@pheur$KP)[[2]]
    ## warning(nmsX)
    nmsX[1] <- "Width"
    nmsX[5] <- "C*M"
    if (nF == 1)
        if (nT == 1)
            pos <- list(-1, nT + nT)
        else
            pos <- list(-1, nT, nT + nT)
    else
        pos <- list(-1, nT, nT + (1:nT) * nF)
    labelsX <- labelsXfun(nmsX, pos)
    strCaption <- paste0("\\textbf{Potency}")
    print(xtable(as.table(X), digits = digits,
                 caption = strCaption, label = "Potency"),
          size = "footnotesize",
          include.colnames = FALSE,
          hline.after = NULL,
          caption.placement = "top",
          add.to.row = list(pos = pos, command = labelsX))

    if (showRelative & (max(dim(Fits@relPotency)) > 0)) {
        X <- Fits@relPotency
        nmsX <- dimnames(X)[[2]]
        if (nF == 1)
            pos <- list(-1, nT + nT)
        else
            pos <- list(-1, (1:(nT-0)) * nF)
        labelsX <- labelsXfun(nmsX, pos)
        strCaption <- paste0("\\textbf{Potency relative estimates of source}")
        print(xtable(as.table(X), digits = digits, caption = strCaption,
                     label = "Relative Potency"),
              size = "footnotesize",
              include.colnames = FALSE,
              hline.after = NULL,
              caption.placement = "top",
              add.to.row = list(pos = pos, command = labelsX))
    }
}
