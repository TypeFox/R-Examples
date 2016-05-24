assayTable2frame <-
function (table,
          dr = 2,
          Z = log((1/dr)^(max(Dilution) - Dilution)),
          byOrder = TRUE,
          echoData = TRUE, ...) 
{
    myLine <- function(char = "=", length = 80, start = "", end = "\n") {
        cat(start)
        cat(paste(rep(char, 80), collapse = ""))
        cat(end)
    }
    if (echoData) {
        myLine()
        cat("Data values:\n")
        myLine("-")
        print(rbind(table,
                    Mean =      apply(table, 2, mean),
                    SD   = sqrt(apply(table, 2,  var)),
                    CV   = sqrt(apply(table, 2,  var))/
                                apply(table, 2, mean) * 100), digits = 2)
        myLine()
    }
    Response   <- c(table)
    D          <- dim(table)
    nms        <- rep(dimnames(table)[[2]], rep(D[1], D[2]))
    splitNms   <- strsplit(nms, ":")
    Sample     <- unlist(lapply(splitNms,
                                FUN = function (x)
                                ifelse(length(x) > 1, x[1],
                                       substr(x, 1, 1)
                                       )))
    DoseCh     <- unlist(lapply(splitNms,
                                FUN = function (x)
                                ifelse(length(x) > 1,
                                       ifelse(!byOrder & (length(x) > 2),
                                              x[3], x[2]),
                                       substr(x, 2, nchar(x))
                                       )))
    try(Dose   <- as.numeric(as.character(DoseCh)), silent = TRUE)
    if (any(is.na(Dose)))
        Dose   <- as.numeric(DoseCh)
    Dilution   <- as.numeric(DoseCh)
    if (!byOrder & !is.integer(Dose))
        Z      <- Dose
    SampleStep <- paste0(Sample, ":", Dilution)
    Replicate  <- rep(1:D[1], D[2])
    ## X <- data.frame(Replicate, Response = Response)
    ## m <- mean(unlist(X["Response"]))
    ## AdjustmentStep  <- m - fitted(lm(Response ~ factor(Dilution),  data = X))
    ## AdjustmentPlate <- m - fitted(lm(Response ~ factor(Replicate), data = X))
    data <- data.frame(Replicate,
                       SampleStep,
                       Sample,
                       Dilution,
                       Dose,
                       Z = Z,
                       Response = Response,
                       I = 1:length(Sample))
    return(data)
}
