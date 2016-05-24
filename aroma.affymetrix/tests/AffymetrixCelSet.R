library("aroma.affymetrix")

## Empty CEL set
ds <- AffymetrixCelSet()
print(ds)

## Extract non-existing subset on empty set
ds <- AffymetrixCelSet()
dsT <- extract(ds, "foo", onMissing="NA")
print(dsT)

## CEL set with non-existing CEL file
files <- list(AffymetrixCelFile())
ds <- AffymetrixCelSet(files)
print(ds)


if (setupExampleData(aroma.affymetrix, mustWork=FALSE)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- AffymetrixCdfFile$byChipType("HG-Focus")
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CEL set / CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cels <- AffymetrixCelSet$byName("FusionSDK_HG-Focus", cdf=cdf)
  print(cels)

  # A single file
  cel <- cels[[1]]
  print(cel)

  # File checksums
  md5s <- getChecksumFileSet(cels)
  print(md5s)

  # File checksums
  md5 <- getChecksumFile(cel)
  print(md5)

  # Average CEL file
  celR <- getAverageFile(cels)
  print(celR)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AromaPlatformInterface
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  platform <- getAromaPlatform(cdf)
  print(platform)
  stopifnot(equals(getAromaPlatform(cels), platform))
  stopifnot(equals(getAromaPlatform(cel), platform))

  # Unit names file
  unf <- getUnitNamesFile(cdf)
  print(unf)
  stopifnot(equals(getUnitNamesFile(cels), unf))
  stopifnot(equals(getUnitNamesFile(cel), unf))
} # if (require("AffymetrixDataTestFiles"))
