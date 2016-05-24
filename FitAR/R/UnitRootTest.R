UnitRootTest <- function(z, method=c("MLE", "ExactMLE", "LS", "All"), statistic=c("Z","T"), NumBoot=1000, PValueMethod=c("DH", "ET")) {
if (length(z)<10)
    stop("Minimum length of series recommended is 10!")
#which method
mthd <- match.arg(method)
IsValidmethodQ <- mthd %in% c("MLE", "ExactMLE", "LS", "All")
if (!IsValidmethodQ)
    stop("method = ", mthd, " not known.")
MeanMLEQ <- FALSE
if (mthd == "ExactMLE")
    MeanMLEQ <- TRUE
if (mthd=="MLE" || MeanMLEQ)
    MLEQ <- TRUE
else
    MLEQ <- FALSE
if (mthd == "All")
    ALLQ <- TRUE
else
    ALLQ <- FALSE
#
stat <- match.arg(statistic)
IsValidstatisticQ <- stat %in% c("Z","T")
if (!IsValidstatisticQ)
    stop("statistic = ", stat, " not known.")
if (stat=="Z")
    UseZQ <- TRUE
else
    UseZQ <- FALSE
#
#method for p-value computation
pvm <- match.arg(PValueMethod)
IsValidpvmQ <- pvm %in% c("DH", "ET")
if (!IsValidpvmQ)
    stop("PValueMethod = ", pvm, " not known.")
if (pvm=="DH")
    UseDHQ <- TRUE
else
    UseDHQ <- FALSE
#
if (ALLQ){
    #Observed value of test statistics
        ansObs1 <- FitARz(z, p=1, MeanMLEQ=FALSE)
        ansObs2 <- FitARz(z, p=1, MeanMLEQ=TRUE)
        ansObs3 <- FitARp(z, p=1, MLEQ=FALSE)
        ZObs1 <- getRho(ansObs1)
        ZObs2 <- getRho(ansObs2)
        ZObs3 <- getRho(ansObs3)
        TObs1 <- getT(ansObs1)
        TObs2 <- getT(ansObs2)
        TObs3 <- getT(ansObs3)
        XObs <- c(ZObs1, ZObs2, ZObs3, TObs1, TObs2, TObs3)
#Bootstrap iterations
    K <-numeric(6)
    for (i in 1:NumBoot){
        zB <- cumsum(rnorm(length(z)))
        ansB1 <- FitARz(zB, p=1, MeanMLEQ=FALSE)
        ansB2 <- FitARz(zB, p=1, MeanMLEQ=TRUE)
        ansB3 <- FitARp(zB, p=1, MLEQ=FALSE)
        ZB1 <- getRho(ansB1)
        ZB2 <- getRho(ansB2)
        ZB3 <- getRho(ansB3)
        TB1 <- getT(ansB1)
        TB2 <- getT(ansB2)
        TB3 <- getT(ansB3)
        XB <- c(ZB1, ZB2, ZB3, TB1, TB2, TB3)
        K <- K + as.numeric(XB<=XObs)
        }
    if (UseDHQ)
        pvalue <- (K+1)/(NumBoot+1)
    else
        pvalue <- K/NumBoot
}
else { 
    #Observed value of test statistic
    if (MLEQ)
        ansObs <- FitARz(z, p=1, MeanMLEQ=MeanMLEQ)
    else
        ansObs <- FitARp(z, p=1, MLEQ=FALSE)
    if (stat=="Z")
        XObs <- getRho(ansObs)
    else
        XObs <- getT(ansObs)
#Bootstrap iterations
    K <-0
    for (i in 1:NumBoot){
        zB <- cumsum(rnorm(length(z)))
        if (MLEQ)
            ansB <- FitARz(zB, p=1, MeanMLEQ=MeanMLEQ)
        else
            ansB <- FitARp(zB, p=1, MLEQ=FALSE)
        if (UseZQ)
            XB <- getRho(ansB)
        else
            XB <- getT(ansB)
        if (XB <= XObs) 
            K <- K+1
        }
    if (UseDHQ)
        pvalue <- (K+1)/(NumBoot+1)
    else
        pvalue <- K/NumBoot
    }
pvalue
}
