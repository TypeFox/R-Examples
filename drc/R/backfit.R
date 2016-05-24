backfit <- function(drcObject)
{
    DL <- drcObject$dataList
    DLdose <- DL$dose
    meansVec <- tapply(DL$origResp, DLdose, mean, na.rm = TRUE) 
    # arranged according to ascending dose values
    # therefore unique doses are sorted below

    backfitValues <- ED(drcObject, meansVec, type = "absolute", display = FALSE, multcomp = TRUE)[["EDdisplay"]][, 1]

    return(cbind(dose = sort(unique(DLdose)), backfit = backfitValues))
}