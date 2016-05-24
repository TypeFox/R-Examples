checkCTechSTu <- function(CTechST, MID){
#checkCTechSTu <- function(CTechST, MID), output: list(CTechST, Status)
#checking the input (list) CTechST for errors
#MID        ... 1 (resp. 2) for the CTechST related to compContourM1u (resp. compContourM2u)
#Status = 0 if the check was successful
#Status = 1 if CTechST was empty or not a list and had to be replaced
#Status > 1 if an important field of CTechST was wrong and had to be replaced
#Note: missing or faulty fields are usually substituted with the default ones. Only
# the faulty fields result in Status > 1
Status <- 0

if (MID == 1){DefCTechST <- getCTechSTM1u()} else {DefCTechST <- getCTechSTM2u()}

FieldNameCA <- names(DefCTechST)

if (!is.list(CTechST)){
    CTechST <- DefCTechST
    Status <- 1
    return(list(CTechST, Status))
} #if

for (IndFNCA in seq_along(FieldNameCA)){
    FieldName <- FieldNameCA[IndFNCA]
    if (is.null(CTechST[[FieldName]])){
        CTechST[[FieldName]] <- DefCTechST[[FieldName]]
        next
    }
    if (length(grep("I$", FieldName))){
        OutList <- checkArray(CTechST[[FieldName]], 1, 0, c(0, 1), c(1, 1), c(1, 1), 1)
        if (OutList[[2]] > 0){
            CTechST[[FieldName]] <- DefCTechST[[FieldName]]
            Status <- IndFNCA + 1
        } #if
    } #if
    if (length(grep("S$", FieldName))){
        StatusF <- !is.character(CTechST[[FieldName]])
        if (StatusF > 0){
            CTechST[[FieldName]] <- DefCTechST[[FieldName]]
            Status <- IndFNCA + 1
        } #if
    } #if
} #for InfFNCA
return(list(CTechST, Status))
}