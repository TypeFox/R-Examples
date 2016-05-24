.imputeRBD <- function(table,
                       epsilon = 1e-15,
                       maxit = 25,
                       trace = FALSE, version = 2, ...) {
    imputeOne1 <- function(table, row, col) {
        Table <- table
        Table[row, col] <- NA
        Table <- cbind(Table, apply(Table, 1, sum, na.rm = TRUE))
        Table <- rbind(Table, apply(Table, 2, sum, na.rm = TRUE))
        B     <- Table[nrow + 1, col]
        T     <- Table[row, ncol+1]
        G     <- Table[nrow+1, ncol+1]
        ydot  <- (nrow * B + ncol * T - G) / ((nrow-1)*(ncol-1))
        table[row, col] <- ydot
        if (trace)
            print(prettyNum(c(1, row, col, B, T, G, ydot)))
        return(table)
    }
    imputeOne2 <- function(table, row, col) {
        G     <- sum(table) - table[row, col]
        B     <- sum(table[, col]) - table[row, col]
        T     <- sum(table[row, ]) - table[row, col]
        ydot  <- (nrow * B + ncol * T - G) / ((nrow-1)*(ncol-1))
        table[row, col] <- ydot
        if (trace)
            print(prettyNum(c(2, row, col, B, T, G, ydot)))
        return(table)
    }
    imputeTable <- function(table) {
        for (i in 1:length(na.row))
            if (version == 1)
                table <- imputeOne1(table, na.row[i], na.col[i])
            else
                table <- imputeOne2(table, na.row[i], na.col[i])
        return(table)
    }
    iterate <- function(table) {
        j <- 0
        delta <- 1
        while ((epsilon < delta | j < 2) & j < maxit) {
            j <- j + 1
            table0 <- table
            table <- imputeTable(table)
            delta <- max(abs(table0 - table) / table, na.rm = TRUE)
            if (trace)
                print(prettyNum(c(0, j, delta)))
        }
        return(table)
    }
    nrow   <- nrow(table)
    ncol   <- ncol(table)
    na.row <- row(table)[is.na(table)]
    na.col <- col(table)[is.na(table)]
    table[is.na(table)] <- mean(table, na.rm = TRUE)
    table  <- iterate(table)
    return(table)
}
