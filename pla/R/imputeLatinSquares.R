.imputeLatinSquare <- function(data,
                               response = "Response",
                               row = "Row",
                               column = "Column",
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
        Arc <- split(data[, c(response, column)], data[, row])
        Brc <- unlist(lapply(Arc, FUN = function(x)
                             split(x[, response], x[, column])))
        k <- sqrt(length(Brc))
        tableRC <- matrix(Brc, byrow = TRUE, ncol = k)
        Asc <- split(data[, c(response, column)], data[, sampleStep])
        Bsc <- unlist(lapply(Asc, FUN = function(x)
                             split(x[, response], x[, column])))
        tableSC <- matrix(Bsc, byrow = TRUE, ncol = k)
        cells <- data[is.na(data[response]),
                      c(sampleStep, row, column, sample, step, response,
                        "Index")]
        iR <- unlist(cells[row])
        iC <- unlist(cells[column])
        iS <- as.numeric(unlist(cells[sampleStep]))
        l <- dim(cells)[1]
        gm <- mean(tableRC, na.rm = TRUE)
        for (i in 1:l) {
            tableRC[iR[i], iC[i]] <- gm
            tableSC[iS[i], iC[i]] <- gm
        }
        j <- 0
        delta <- 1
        while ((epsilon < delta | j < 2) & j < maxit) {
            j <- j + 1
            table0 <- tableRC
            for (i in 1:l) {
                B <- sum(tableRC[iR[i],]) - tableRC[iR[i], iC[i]]
                C <- sum(tableRC[,iC[i]]) - tableRC[iR[i], iC[i]]
                T <- sum(tableSC[iS[i],]) - tableSC[iS[i], iC[i]]
                G <- sum(tableSC        ) - tableSC[iS[i], iC[i]]
                yDot <- (k*(B+C+T) - 2*G) / ((k-1)*(k-2))
                cells[i, response] <- yDot
                tableRC[iR[i], iC[i]] <- yDot
                tableSC[iS[i], iC[i]] <- yDot
            }
            delta <- max(abs(table0 - tableRC) / tableRC, na.rm = TRUE)
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
