library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-4, timestamp=TRUE)

# Avoid being masked by affy::plotDensity()
plotDensity <- aroma.light::plotDensity

dataSet <- "GSE8605"
chipType <- "Mapping10K_Xba142"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType)
print(csR)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv2(csR, drop=FALSE, verbose=verbose)
print(res)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- res$plm
print(plm)

ces <- getChipEffectSet(plm)
print(ces)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quality scores via ChipEffectSet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Boxplots of log2(theta), RLE, and NUSE
toPNG(getFullName(ces), tags=c("QA", "theta", "RLE", "NUSE"), {
  layout(matrix(1:4, ncol=2, byrow=TRUE))
  plotBoxplot(ces, type="theta", transform=log2)
  plotBoxplot(ces, type="RLE", arrays=c(2,4:6))
  plotBoxplot(ces, type="NUSE")
})

# Calculating statistics without plotting
theta <- boxplotStats(ces, type="theta", transform=log2)
nuse <- boxplotStats(ces, type="NUSE")
rle <- boxplotStats(ces, type="RLE")

# Subset of arrays (avoids calculating stats for all arrays)
rleB <- boxplotStats(ces, type="RLE", arrays=c(2,4:6))
for (name in names(rleB)) {
  stopifnot(identical(rleB[name], rle[name]))
}

# Plotting the above statistics
toPNG(getFullName(ces), tags=c("QA", "theta", "RLE", "NUSE", "manual"), {
  layout(matrix(1:4, ncol=2, byrow=TRUE))
  plotBoxplotStats(theta, main="theta")
  plotBoxplotStats(rle[c(2,4:6)], main="RLE")
  plotBoxplotStats(nuse, main="NUSE")
})


# Calculates unit-specific RLE and NUSE scores
units <- 1000+1:5000
theta <- extractMatrix(ces, units=units)
rle <- extractMatrix(ces, units=units, field="RLE")
nuse <- extractMatrix(ces, units=units, field="NUSE")

# Plotting the above statistics
unitsTag <- sprintf("units%s", seqToHumanReadable(units))
toPNG(getFullName(ces), tags=c(unitsTag, "QA", "theta", "RLE", "NUSE"), {
  layout(matrix(1:4, ncol=2, byrow=TRUE))
  plotDensity(log2(theta), main="theta")
  plotDensity(rle[,c(2,4:6)], main="RLE")
  plotDensity(nuse, main="NUSE")
})

# ...same, but basic unit annotation data added
theta <- extractDataFrame(ces, units=units, addNames=TRUE)
print(head(theta))
rle <- extractDataFrame(ces, units=units, field="RLE", addNames=TRUE)
print(head(rle))
nuse <- extractDataFrame(ces, units=units, field="NUSE", addNames=TRUE)
print(head(nuse))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quality scores via QualityAssessmentModel
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
qa <- QualityAssessmentModel(plm)
print(qa)

toPNG(getFullName(qa), tags=c("RLE", "NUSE"), {
  subplots(4, ncol=2, byrow=FALSE)
  par(mar=c(6,5,3,1))
  plotRle(qa)
  plotNuse(qa)
  # Reordered along x-axis
  arrays <- nbrOfArrays(qa):1
  plotRle(qa, arrays=arrays)
  plotNuse(qa, arrays=arrays)
})
