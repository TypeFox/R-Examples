Cpredict <- function (focal, Krigobject) {
    CKrigidx <- Krigobject$CKrigidx
    if (is.null(CKrigidx)) {
        stop.redef("(!) NULL 'CKrigidx' value in Cpredict(...)")
    }
    if (CKrigidx < 0) {
        stop.redef(paste("(!) Invalid 'CKrigidx' value", CKrigidx,
            "in Cpredict(...)"))
    }
    if (CKrigidx > 100) {
        message.redef(paste("(!) Suspicious 'CKrigidx' value",
            CKrigidx, "in Cpredict(...)"))
        message.redef(paste("(!) Might be correct, but then still suggests a poor organization of code."))
    }
    covfocal <- CcovFocal(focal,CKrigidx)
    return(Krigobject$d + sum(covfocal * Krigobject$c))
}
