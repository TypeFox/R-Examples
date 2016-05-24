SaveVegChangeProbabilityMultipliers <-
function(vegChanges, project, climateChangeTransitionTypes, vt2pvtlut) {
# vegChanges is created by this statement at the end of VegTypeChanges.R
# "return(list(tgtFile, years, vts2keep, vtFracsReduced, changeFracsReduced))"
#
# Here we write out one .txt file containing the transition probability multipliers for all
# the climate change transition types. 
# The file is named "vegChangeProbabilityMultipliers.txt"
# The file consists of a header line followed by one line for each timestep for each transition type.
# The header line is just the incoming "tgtfile" string.
# Each succeeding line looks like:
# <timestep><tab>"temporal"<tab><transition type><tab><multiplier value>
# e.g.
# 1<tab>temporal<tab>fdd2fsi<tab>0
# 2<tab>temporal<tab>fdd2fsi<tab>0
# ...
# 14<tab>temporal<tab>fdd2fsi<tab>10.6597
# ...
# The mean over all the time steps of the transition probability for each transition type is written
# to the console.

srcDataFile = vegChanges[[1]]
years = vegChanges[[2]]
VTs = vegChanges[[3]]
nVT = length(VTs)
changeFracs = vegChanges[[5]]
nYrs = dim(changeFracs)[3]
stopifnot((nYrs + 1)==length(years))
nTransitionTypes = length(climateChangeTransitionTypes)
stopifnot(nTransitionTypes>=1)

multiplierFile = "vegChangeProbabilityMultipliers.txt"

#pvts = VTpvts(project)

pvts<-vt2pvtlut$PVT

# header line
cat(srcDataFile, file=multiplierFile, append=FALSE)
cat("\n", file=multiplierFile, append=TRUE)

# all succeeding lines, one for each VTsrc->VTdest pair
for (kSrc in 1:nVT) {
   for (kDest in 1:nVT) if (kSrc!=kDest) {
      meanTransitionProbability = mean(changeFracs[kSrc, kDest,])

      transitionType = paste(c(levels(pvts)[pvts[kSrc]]), "2", levels(pvts)[pvts[kDest]], sep="")
      iType = 0
      found = FALSE
      while (!found && iType<nTransitionTypes) {
         iType = iType + 1
         found = climateChangeTransitionTypes[iType]==transitionType
         }
      if (!found) next
      
      cat(c(transitionType, meanTransitionProbability, "\n"))

      # Now <year><tab>"temporal"<tab><transition type><tab><multiplier value> for each year
      for (yr in 1:nYrs) {
         if (meanTransitionProbability>0) transitionProbabilityMultiplier = 
               changeFracs[kSrc, kDest, yr]/meanTransitionProbability
         else transitionProbabilityMultiplier = 0
         
         cat(yr, file=multiplierFile, append=TRUE)
         cat("\tTemporal\t", file=multiplierFile, append=TRUE)
         cat(transitionType, file=multiplierFile, append=TRUE)
         cat("\t", file=multiplierFile, append=TRUE)
         cat(transitionProbabilityMultiplier, file=multiplierFile, append=TRUE)
         cat("\n", file=multiplierFile, append=TRUE)
         } # end of loop on yr


      } # end of for (kDest in 1:nVT) if (kSrc!=kDest) ...
   } # end of loop on kSrc

}
