`ULPrint` <- 
function(aData, nRows=10, nCols=10, Title="") {
    if (Title != "") cat(Title,"\n")
    RRng <- RCRng(nRows)
    CRng <- RCRng(nCols) 
    print(aData[
        max(RRng[1],1):min(nrow(aData),RRng[2]),
        max(CRng[1],1):min(ncol(aData),CRng[2])])
    if (nRows < nrow(aData)) cat("...")
    cat("\n")
}
