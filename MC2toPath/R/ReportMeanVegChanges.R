ReportMeanVegChanges <-
function(baseCalib, VTs, vtFracs, changeFracs, vtXpvt=data.frame(NULL)) {

vtNames = VTnames(baseCalib)
nVT = length(VTs)
nYrs = dim(changeFracs)[3]

for (kSrc in 1:nVT) {
   meanFracOfAllCells = mean(vtFracs[kSrc, ])
   if (meanFracOfAllCells>0) {
      # cat(c("\nAveraged over ", nYrs, " years, VTYPE ", VTs[kSrc], " occupies ", meanFracOfAllCells, " fraction of all cells.\n"))
      if (length(vtXpvt)>0) {
         cat(c("\n", "mean transition probabilities for transitions out of", levels(vtXpvt$PVT)[vtXpvt$PVT[kSrc]], "...\n"))
         for (kDest in 1:nVT) if (kSrc!=kDest) {
            meanTransitionProbability = mean(changeFracs[kSrc, kDest,]) 
            if (meanTransitionProbability>0) {
#               cat(c("from VTYPE ", VTs[kSrc], " to VTYPE ", VTs[kDest], "the mean transition probability is ", meanTransitionProbability, "\n"))
#               cat(c("from PVT ", levels(vtXpvt$PVT)[vtXpvt$PVT[kSrc]], " to PVT ", levels(vtXpvt$PVT)[vtXpvt$PVT[kDest]], "the mean transition probability is ", meanTransitionProbability, "\n"))
               cat(c(levels(vtXpvt$PVT)[vtXpvt$PVT[kSrc]], "2", levels(vtXpvt$PVT)[vtXpvt$PVT[kDest]], meanTransitionProbability, "\n"))
            }
         } # end of if (length(vtXpvt>0) ...
      } # end of for (kDest in 1:nVT) if (kSrc!=kDest) ...
   } # end of if (meanFracOfAllCells>0) ...
} # end of loop on kSrc

cat("\n")
cat(c("Mean values over ", dim(vtFracs)[2], " years:\n"))
cat("frac of all cells, VTYPE\n")
for (kSrc in 1:nVT) {
   meanFracOfAllCells = mean(vtFracs[kSrc, ])
   if (meanFracOfAllCells>0) {
      vtName = vtNames[[VTs[kSrc]]]
      cat(c(meanFracOfAllCells, VTs[kSrc], vtName, "\n"))
      }
} # end of loop on kSrc

} # end of ReportMeanVegChanges()
