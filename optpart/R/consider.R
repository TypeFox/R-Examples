consider <- function (part) 
{
    if (!inherits(part,'partana')) stop('Must pass argument of class partana')

    ctc <- part$ctc
    vals <- part$ctc[row(part$ctc)>col(part$ctc)]
    mus <- as.numeric(table(part$clustering))
    n <- length(part$clustering)

    row <- rep(0,length(vals))
    col <- rep(0,length(vals))

    pnt <- 0
    for (i in 1:(ncol(ctc)-1)) {
        for (j in (i+1):nrow(ctc)) {
            pnt <- pnt + 1
            row[pnt] <- j
            col[pnt] <- i
        }
    }
    row <- row[rev(order(vals))]
    col <- col[rev(order(vals))]
    vals <- rev(sort(vals))
    out <- data.frame(row,col,vals)
    out
}
