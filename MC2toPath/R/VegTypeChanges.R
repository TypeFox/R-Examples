VegTypeChanges <-
function(tgtFile, baseCalibration, vtXpvt=data.frame(NULL)) {
tgtP = open.nc(tgtFile)

tgtLonInfo = dim.inq.nc(tgtP, "lon")
tgtLonDimID = tgtLonInfo$id
tgtLatInfo = dim.inq.nc(tgtP, "lat")
tgtLatDimID = tgtLatInfo$id
tgtYrInfo = dim.inq.nc(tgtP, "year")
tgtYrDimID = tgtYrInfo$id
years = var.get.nc(tgtP, "year")
nVTall = length(VTnames(baseCalibration))	
if (length(vtXpvt)>0) dontskipVTs = vtXpvt$VT
else dontskipVTs = c(1:nVTall)

tgtVarName = "VTYPE"
tgtVarInfo = var.inq.nc(tgtP, tgtVarName)
tgtVarDimIds = tgtVarInfo$dimids
stopifnot(length(tgtVarDimIds)==3)
stopifnot(tgtLonDimID==tgtVarDimIds[1])
stopifnot(tgtLatDimID==tgtVarDimIds[2])
stopifnot(tgtYrDimID==tgtVarDimIds[3])

tgtVar = var.get.nc(tgtP, tgtVarName)
tgtDim = dim(tgtVar)
stopifnot(length(tgtDim)==3)
nCols = tgtDim[1]
nRows = tgtDim[2]
nYrs = tgtDim[3]
stopifnot(nYrs>=2)

nCells = nCols*nRows
dim(tgtVar) = c(nCells, nYrs)

# Count the cells in each vegtype in each year
cat("Count the cells in each vegetation type in each year...\n")
vtCounts = array(0, c(nVTall, nYrs))
for (yr in 1:nYrs) vtCounts[, yr] = tabulate(tgtVar[,yr], nbins=nVTall)   

# Figure out which vegtypes change to which vegtypes
cat("Figure out which vegtypes change to which vegtypes...\n")
changePairs = array(FALSE, c(nVTall, nVTall))
for (yr in 2:nYrs) {
   cat(c("year = ", years[yr], "\n"))
   for (cell in 1:nCells) {
      vtPrev = tgtVar[cell, yr - 1] 
      if (!is.na(vtPrev) && vtPrev==0) vtPrev = NA;
      if (!((1<=vtPrev && vtPrev<=nVTall) || is.na(vtPrev))) cat(c("vtPrev = ", vtPrev, "\n"))
      stopifnot((1<=vtPrev && vtPrev<=nVTall) || is.na(vtPrev))
      vtCurr = tgtVar[cell, yr]
      if (!is.na(vtCurr) && vtCurr==0) vtCurr = NA;
      if (!((1<=vtCurr && vtCurr<=nVTall) || is.na(vtCurr))) cat(c("vtCurr = ", vtCurr, "\n"))
      stopifnot((1<=vtCurr && vtCurr<=nVTall) || is.na(vtCurr))
      if (!is.na(vtPrev) && !is.na(vtCurr))
            changePairs[vtPrev, vtCurr] = changePairs[vtPrev, vtCurr] || (vtPrev!=vtCurr)
      } # end of loop on cell
   } # end of loop on yr

# Figure out which vegtypes can be omitted
cat("Figure out which vegtypes can be omitted...\n")
vts2omit = rep(TRUE, times=nVTall)
if (length(vtXpvt)>0) {
   for (k in 1:length(dontskipVTs)) vts2omit[dontskipVTs[k]] = FALSE
   } # end of if (length(vtXpvt)>0)...
else {
   for (vtPrev in 1:nVTall) {
      for (vtCurr in 1:nVTall) {
         if (changePairs[vtPrev, vtCurr]) vts2omit[vtPrev] = FALSE         
         } # end of loop on vtCurr
      } # end of loop on vtPrev
   } # end of if (length(vtXpvt)>0)... else ...
rm(changePairs)
if (length(vtXpvt)>0) stopifnot(length(vtXpvt$VT)==(nVTall - sum(vts2omit)))

# For each year, first count the number of cells in each transition pair, 
# and then convert the count to a fraction by dividing by the number of cells
# having the original veg type in the first year of the pair of years.
cat("Count the number of cells in each transition pair, etc...\n")
changeCounts = array(0, c(nVTall, nVTall, nYrs - 1))
changeFracs = array(0, c(nVTall, nVTall, nYrs - 1))
for (yr in 2:nYrs) {
   for (cell in 1:nCells) {
      vtPrev = tgtVar[cell, yr - 1]
      vtCurr = tgtVar[cell, yr]
      if (!is.na(vtPrev) && 1<=vtPrev && vtPrev<=nVTall && !is.na(vtCurr) && 1<=vtCurr && vtCurr<=nVTall) 
            changeCounts[vtPrev, vtCurr, yr - 1] = changeCounts[vtPrev, vtCurr, yr - 1] + 1
      } # end of loop on cell
   for (vtPrev in 1:nVTall) {
      for (vtCurr in 1:nVTall) {
         count = changeCounts[vtPrev, vtCurr, yr - 1]
         vtPrevTot = vtCounts[vtPrev, yr - 1]
         if (count>0) {
            stopifnot(vtPrevTot>=1)
            changeFracs[vtPrev, vtCurr, yr - 1] = count/vtPrevTot
            } # end of if (count>0)
         else if (vtPrevTot==0 && yr>2) {
            # There aren't any gridcells left in vtPrev veg type.
            # Replicate the transition multiplier from the previous year.
            # This is to handle the case in Path and Envision where there are actually
            # gridcells remaining in the given type.
            changeFracs[vtPrev, vtCurr, yr - 1] = changeFracs[vtPrev, vtCurr, yr - 2]
            }
         } # end of loop on vtCurr
      } # end of loop on vtPrev
   } # end of loop on yr

# Now reduce the size of the arrays by eliminating vegtypes that never occur.
cat("Reduce the size of the arrays by eliminating vegtypes that never occur...\n")
nVTreduced = nVTall - sum(vts2omit)
vtCountsReduced = array(0, c(nVTreduced, nYrs))
k = 0
for (vt in 1:nVTall) {
   if (!vts2omit[vt]) {
      k = k + 1
      vtCountsReduced[k, ] = vtCounts[vt, ]
      } # end of if (!vts2omit...
   } # end of loop on vt
vtFracsReduced = array(0, c(nVTreduced, nYrs))
maxRoundOffErr = 0
for (yr in 1:nYrs) {
   totCounts = sum(vtCountsReduced[ , yr])
   vtFracsReduced[ , yr] = vtCountsReduced[ , yr]/totCounts
   if (abs(sum(vtFracsReduced[,yr]) - 1)>abs(maxRoundOffErr)) maxRoundOffErr = sum(vtFracsReduced[,yr]) - 1
   } # end of loop on yr
cat(c("maxRoundOffErr = ", maxRoundOffErr, "\n"))
changeFracsReducedByRows = array(0, c(nVTall, nVTreduced, nYrs - 1))
vts2keep = rep(0, times=nVTreduced)
for (yr in 1:(nYrs - 1)) {
   k = 0
   for (vt in 1:nVTall) {
      if (!vts2omit[vt]) {
         k = k + 1
         vts2keep[k] = vt
         changeFracsReducedByRows[, k, yr] = changeFracs[, vt, yr]
         } # end of if (!vts2omit...
      } # end of loop on vt
   } # end of loop on yr
changeFracsReduced = array(0, c(nVTreduced, nVTreduced, nYrs - 1))
k = 0
for (vt in 1:nVTall) {
   if (!vts2omit[vt]) {
      k = k + 1
      changeFracsReduced[k, ,] = changeFracsReducedByRows[vt, ,]
      } # end of if (!vts2omit...
   } # end of loop on vt

ReportMeanVegChanges(baseCalibration, vts2keep, vtFracsReduced, changeFracsReduced, vtXpvt)
return(list(tgtFile = tgtFile,years =  years,vts2keep =  vts2keep,vtFracsReduced= vtFracsReduced, changeFracsReduced = changeFracsReduced))

}
