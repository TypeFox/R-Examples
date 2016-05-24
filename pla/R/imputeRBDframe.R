.imputeRBDframe <- function(data,
                            response = "Response",
                            replicate = "Replicate",
                            sample = "Sample",
                            step = "Dilution",
                            sampleStep = "SampleStep",
                            epsilon = 1e-15,
                            maxit = 25,
                            trace = FALSE, ...) {
    data <- cbind(data, Index = 1:dim(data)[1])
    Response <- data[response]
    if (any(is.na(Response))) {
        if (!any(dimnames(data)[[2]] == sampleStep)) {
            SampleStep <- paste0(as.character(unlist(data[sample])), ":",
                                 as.character(unlist(data[step])))
            names(SampleStep) <- "SampleStep"
            data <- cbind(data, SampleStep = SampleStep)
        }
        Asr <- split(data[, c(response, replicate)], data[, sampleStep])
        ncol = length(Asr)
        Bsr <- unlist(lapply(Asr, FUN = function(x)
                             split(x[, response], x[, replicate])))
        tableSR <- matrix(Bsr, byrow = TRUE, ncol = ncol)
        nrow = nrow(tableSR)
        cells <- data[is.na(data[response]),
                      c(sampleStep, replicate, sample, step, response, "Index")]
        iR <- unlist(cells[replicate])
        iS <- as.numeric(unlist(cells[sampleStep]))
        l <- dim(cells)[1]
        gm <- mean(tableSR, na.rm = TRUE)
        for (i in 1:l)
            tableSR[iS[i], iR[i]] <- gm
        j <- 0
        delta <- 1
        while ((epsilon < delta | j < 2) & j < maxit) {
            j <- j + 1
            table0 <- tableSR
            for (i in 1:l) {
                B <- sum(tableSR[,iR[i]]) - tableSR[iS[i], iR[i]]
                T <- sum(tableSR[iS[i],]) - tableSR[iS[i], iR[i]]
                G <- sum(tableSR        ) - tableSR[iS[i], iR[i]]
                yDot  <- (nrow * B + ncol * T - G) / ((nrow-1)*(ncol-1))
                cells[i, response] <- yDot
                tableSR[iS[i], iR[i]] <- yDot
            }
            delta <- max(abs(table0 - tableSR) / tableSR, na.rm = TRUE)
            if (trace) {
                print(cells)
                print(prettyNum(c(j, delta)))
            }
        }
        for (i in 1:l)
            data[cells[i, "Index"], response] <- cells[i, response]
    }
    return(data)
}
