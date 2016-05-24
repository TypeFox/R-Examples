VAR <-
function (y, p = 1,  exogen = NULL )
{
  y <- as.matrix(y)
    if (any(is.na(y)))
        stop("\nNAs in y.\n")
    if (ncol(y) < 2)
        stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
    if (is.null(colnames(y))) {
        colnames(y) <- paste("y", 1:ncol(y), sep = "")
            }
    colnames(y) <- make.names(colnames(y))

    obs <- dim(y)[1]
    K <- dim(y)[2]

    sample <- obs - p
    ylags <- embed(y, dimension = p + 1)[, -(1:K)]


    colnames(ylags) <- paste(colnames(y), ".l", rep(1:p, each = ncol(y)), sep = "")
    yend <- y[-c(1:p), ]
        rhs <- ylags

    if (!(is.null(exogen))) {
        exogen <- as.matrix(exogen)
        if (!identical(nrow(exogen), nrow(y))) {
            stop("\nDifferent row size of y and exogen.\n")
        }
        if (is.null(colnames(exogen))) {
            colnames(exogen) <- paste("exo", 1:ncol(exogen),
                sep = "")}

        colnames(exogen) <- make.names(colnames(exogen))
        rhs <- cbind(rhs, exogen[-c(1:p), ])
           }
    datamat <- as.data.frame(rhs)

    y <- yend[,1:K]
    equation <- lm(y~ -1+., data = datamat)
   return(t(equation$coefficients))
}
