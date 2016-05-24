library("aroma.affymetrix")

## Empty CDF
cdf <- AffymetrixCdfFile()
print(cdf)

## Missing CDF
cdf <- AffymetrixCdfFile(NA_character_, mustExist=FALSE)
print(cdf)


if (setupExampleData(aroma.affymetrix, dataset="FusionSDK_Test3", dirs="annotationData", mustWork=FALSE)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- AffymetrixCdfFile$byChipType("Test3")
  print(cdf)

  # Unit names, types, ...
  names <- getUnitNames(cdf)
  str(names)

  units <- indexOf(cdf, names=names[42:40])
  stopifnot(all.equal(units, 42:40))

  types <- getUnitTypes(cdf)
  print(table(types))

  data <- readUnits(cdf, units=1:10)
  str(data)

  data <- readDataFrame(cdf, units=1:10)
  str(data)

  md5 <- getChecksumFile(cdf)
  print(md5)

  ## Get an existing or create a new onocell CDF
  cdfM <- getMonocellCdf(cdf, verbose=-10)
  print(cdfM)
  t0 <- lastModified(getPathname(cdfM))
  print(t0)

  ## A monocell CDF can be re-created
  cdfM <- getMonocellCdf(cdf, force=TRUE)
  print(cdfM)
  t1 <- lastModified(getPathname(cdfM))
  print(t1)
  stopifnot(t1 >= t0)
} # if (... "AffymetrixDataTestFiles")
