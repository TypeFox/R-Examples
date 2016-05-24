IntermediateOO <-
function(plist, OOCorrMat){
validate.plist(plist,ncol(OOCorrMat))
return( ordcont(plist, OOCorrMat)$SigmaC )

}
