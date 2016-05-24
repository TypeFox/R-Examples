## 2014-04-17, 27, 2014-10-30, 2015-11-26
validate <- function (x, test, validrange = c(0, Inf), targets = test) {
    onescenario <- function (scen, minx, maxx, scennum) {
        criterion <- (scen[,test] < minx) | (scen[,test] > maxx)
        criterion[is.na(criterion)] <- TRUE     ## 2014-10-30
        cat('Dropping ', sum(criterion) , 'value(s) for scenario ', scennum, '\n')
        scen[criterion,targets] <- NA
        scen[is.na(scen[,test]),targets] <- NA
        scen
    }
    if (!inherits(x, 'selectedstatistics'))
        stop ("requires 'selectedstatistics'")
    if (tolower(targets)[1] == "all")
        targets <- colnames(x$output[[1]])
    if (!all(c(test, targets) %in% colnames(x$output[[1]])))
        stop ("test or targets not found in input")
    nscen <- length(x$output)
    if (!is.matrix(validrange))
        validrange <- matrix (validrange, byrow = TRUE, nrow = nscen, ncol = 2)
    if (nrow(validrange) != nscen)
        stop ("invalid `validrange'")
    minx <- validrange[,1]
    maxx <- validrange[,2]
    x$output <- mapply(onescenario, x$output, minx, maxx, x$scenarios$scenario, SIMPLIFY = FALSE)
    x
}


