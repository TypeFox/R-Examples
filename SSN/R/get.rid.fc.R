get.rid.fc <- function(binIDs, referenceBinID) {
    ind.match <- .Call("test_fc", binIDs, referenceBinID)
    data.frame(
        fc = ind.match<0,
        binaryID = substr(binIDs, 1, abs(ind.match)),
        stringsAsFactors = FALSE)
}
