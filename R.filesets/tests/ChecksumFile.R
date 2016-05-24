library("R.filesets")

message("*** ChecksumFile / ChecksumFileSet")

## Empty / missing
dfZ <- ChecksumFile()
print(dfZ)

dfZ <- ChecksumFile(NA_character_)
print(dfZ)


## Example files
path <- system.file("exData", "dataSetA,original", package="R.filesets")
print(path)

## Setting up a file set
ds <- GenericDataFileSet$byPath(path)
print(ds)

## Create copy (so that we can write checksum files)
pathT <- tempdir()
dsC <- copyTo(ds, path=pathT, overwrite=TRUE)
print(dsC)

## Checksum set
dsCZ <- getChecksumFileSet(dsC)
print(dsCZ)

validate(dsCZ, verbose=TRUE)

print(readChecksums(dsCZ))


message("*** ChecksumFile / ChecksumFileSet ... DONE")
