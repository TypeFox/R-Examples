##Lazy Load instead
## data(DAT.simData, package="DATforDCEMRI")

myccarray <- (DAT.simData$mapCC)
mytimevector <- (DAT.simData$vectorTimes)
myroiarray <- (DAT.simData$maskROI)
myaifvector <- (DAT.simData$vectorAIF)
DAT.checkData(file.name="mydcemridata", vector.times=mytimevector, map.CC=myccarray,
mask.ROI=myroiarray, vector.AIF=myaifvector, slice.stop=1)
DAT(file="mydcemridata_s1-s1.RData", slice=1, range.map=1.05, cutoff.map=0.95,
batch.mode=TRUE)

DAT_file <- list.files(pattern="DAT_mydcemridata_s1-s1_s1")[1]

DAT(DAT_file)

