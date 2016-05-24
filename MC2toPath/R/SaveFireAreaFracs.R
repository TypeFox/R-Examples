SaveFireAreaFracs <-
function(vegChangesLocal, fireFracsLocal) { # , vt2pvtlutLocal) {

srcDataFile = vegChangesLocal[[1]]
years = vegChangesLocal[[2]]
nYr = length(years)
VTs = vegChangesLocal[[3]]
nVT = length(VTs)

vegFracs = vegChangesLocal[[4]] 
# Note:
# vegFracs[nVT, nYr]
# fireFracsLocal[nYr, nVT]

outFile = "fireAreaFracs.csv"
fireAreaFracs = matrix(0, nrow=nYr, ncol=nVT)

#pvts<-vt2pvtlut$PVT

# header line
cat(srcDataFile, file=outFile, append=FALSE)
cat("\n", file=outFile, append=TRUE)

# column headings
cat("year", file=outFile, append=TRUE)
for (vtNdx in 1:nVT) cat(c(", ", VTs[vtNdx]), file=outFile, append=TRUE)
cat("\n", file=outFile, append=TRUE)

# all succeeding lines, one for each year
for (yrNdx in 1:nYr) {
   cat(years[yrNdx], file=outFile, append=TRUE)
   for (vtNdx in 1:nVT) {
      fireAreaFracs[yrNdx, vtNdx] = vegFracs[vtNdx, yrNdx]*fireFracsLocal[yrNdx, vtNdx]
      cat(c(", ", fireAreaFracs[yrNdx, vtNdx]), file=outFile, append=TRUE)
   }
   cat("\n", file=outFile, append=TRUE)
}

return(fireAreaFracs)
} # end of SaveFireAreaFracs()
