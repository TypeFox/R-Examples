clean <-
function (w, tol.col = 0.5, tol.row = 0.3, name = "") 
{
    w = as.data.frame(w)
    w = as.matrix(w)
    if (sum(is.na(w)) == 0) 
        stop("This dataset does not require cleaning.\n",call.=FALSE)
    else {
        filename = paste("Clean.rep.", name, sep = "")
        zz <- textConnection(filename, "w")
        rep.title = paste("Cleaning report for the matrix: ", 
            name)
        sink(zz)
        cat("\n", rep.title, "\n\n")
        sink()
        sumcol = which(colSums(is.na(w)) != 0, arr.ind = TRUE)
        if (length(sumcol) != 0) {
            dr = dim(w)[1]
            dc = dim(w)[2]
            if (length(sumcol) == 1) {
                per.miss.col = sum(is.na(w[, sumcol]))/dr
                colmiss = colnames(w)[sumcol]
                table.miss = data.frame(cbind(Variables = colmiss, 
                  Percent.of.missing = (per.miss.col * 100)), 
                  row.names = NULL)
                print(table.miss)
                cat("\n")
                sink(zz)
                print(table.miss)
                cat("\n")
                sink()
                if (per.miss.col > tol.col) {
                  cat("Only one variable eliminated: ", colnames(w)[above.tol], 
                    "\n\n")
                  sink(zz)
                  cat("Only one variable eliminated: ", colnames(w)[above.tol], 
                    "\n\n")
                  sink()
                  w = w[, -sumcol]
                  w = as.matrix(w)
                }
            }
            else {
                per.miss.col = colSums(is.na(w[, sumcol]))/dr
                above.tol = sumcol[which(per.miss.col > tol.col, 
                  arr.ind = TRUE)]
                colmiss = colnames(w)[sumcol]
                table.miss = data.frame(cbind(Variables = colmiss, 
                  Percent.of.missing = (per.miss.col * 100)), 
                  row.names = NULL)
                print(table.miss)
                cat("\n")
                sink(zz)
                print(table.miss)
                cat("\n")
                sink()
                if (length(above.tol) == dim(w)[2]) {
                  cat("All variables have missing values above tolerance level.\n\n")
                  sink(zz)
                  cat("All variables have missing values above tolerance level.\n\n")
                  sink()
                }
                else if (length(above.tol) != 0) {
                  col.above.tol = matrix(colnames(w)[above.tol], 
                    length(above.tol), 1)
                  colnames(col.above.tol) = "Variables eliminated"
                  rownames(col.above.tol) = c(1:length(above.tol))
                  print(col.above.tol)
                  cat("\n\n")
                  sink(zz)
                  print(col.above.tol)
                  cat("\n\n")
                  sink()
                  w = w[, -above.tol]
                  w = as.matrix(w)
                }
            }
            dr = dim(w)[1]
            dc = dim(w)[2]
        }
        sumrow = which(rowSums(is.na(w)) != 0, arr.ind = TRUE)
        if (length(sumrow) != 0) {
            if (length(sumrow) == 1) {
                per.miss.row = sum(is.na(w[sumrow, ]))/dc
                if (per.miss.row > tol.row) {
                  cat("Number of instances eliminated: 1\n")
                  cat("Instance eliminated              :", sumrow, 
                    "\n\n")
                  sink(zz)
                  cat("Number of instances eliminated: 1\n")
                  cat("Instance eliminated              :", sumrow, 
                    "\n\n")
                  sink()
                  w = w[-sumrow, ]
                }
            }
            else {
                per.miss.row = rowSums(is.na(w[sumrow, ]))/dc
                above.tol = sumrow[which(per.miss.row > tol.row, 
                  arr.ind = TRUE)]
                if (length(above.tol) == dr) 
                  cat("All instances have missing values above tolerance level.\n")
                else if (length(above.tol) != 0) {
                  cat("Number of instances eliminated:", length(above.tol), 
                    "\n")
                  cat("Instance eliminated           :", as.numeric(above.tol), 
                    "\n\n")
                  sink(zz)
                  cat("Number of instances eliminated:", length(above.tol), 
                    "\n")
                  cat("Instance eliminated           :", as.numeric(above.tol), 
                    "\n\n")
                  sink()
                  w = w[-(above.tol), ]
                }
            }
        }
    }
    w = as.matrix(w)
    cat("Maximum number of values to be imputed: ", sum(is.na(w)), 
        "\n")
    sink(zz)
    cat("Maximum number of values to be imputed: ", sum(is.na(w)), 
        "\n")
    sink()
    close(zz)
    return(w)
}
