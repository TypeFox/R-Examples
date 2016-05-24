library("R.filesets")

message("*** GenericDataFile ...")

## Example files
path <- system.file("exData", "dataSetA,original", package="R.filesets")
print(path)

## Setting up a file set
ds <- GenericDataFileSet$byPath(path)
print(ds)

df <- ds[[1L]]
print(df)

cat(sprintf("Created on: %s\n", getCreatedOn(df)))
cat(sprintf("Accessed on: %s\n", getLastAccessedOn(df)))

## Missingness
print(is.na(df))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Comparisons
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathT <- file.path(tempdir(), "copy")
ds <- copyTo(ds, path=pathT, overwrite=TRUE)

df <- ds[[1]]
res <- equals(df, df)
print(res)
stopifnot(res)

res <- equals(df, NULL)
print(res)
stopifnot(!res)

res <- equals(df, NA)
print(res)
stopifnot(!res)

res <- equals(df, ds)
print(res)
stopifnot(!res)

print(df)
print(getChecksum(df, verbose=TRUE))
print(readChecksum(df))

print(df)
print(getChecksum(df, force=TRUE, verbose=TRUE))

print(df)
print(getChecksum(df, write=TRUE, force=TRUE, verbose=TRUE))

print(df)
print(getChecksum(df, write=FALSE, force=TRUE, verbose=TRUE))


pathT <- file.path(tempdir(), "foo")
dfC <- copyTo(df, path=pathT, overwrite=TRUE)
print(dfC)
res <- equals(df, dfC)
print(res)
stopifnot(res)
mod <- hasBeenModified(dfC)
stopifnot(!mod)
print(getChecksum(dfC))

mod1 <- hasBeenModified(dfC, update=FALSE)
print(mod1)
mod2 <- hasBeenModified(dfC, update=FALSE)
print(mod2)
stopifnot(identical(mod2, mod1))


## Make sure we can detect differences in timestamps
Sys.sleep(1.5)
raw <- raw(length=getFileSize(dfC))
writeBin(raw, con=getPathname(dfC))
print(dfC)
print(getChecksum(dfC))

res <- equals(df, dfC)
print(res)
stopifnot(!res)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Checksum files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print(hasChecksumFile(df))

print(getChecksum(df))
print(writeChecksum(df))

dfZ <- getChecksumFile(df, verbose=TRUE)
print(dfZ)
validate(dfZ)

dfZ <- getChecksumFile(df, force=TRUE)
print(dfZ)
validate(dfZ)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# copyTo() / renameTo()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** copyTo() / renameTo() on GenericDataFile")

# Copy it to a temporary directory
path <- tempdir()
dfC <- copyTo(df, path=path, overwrite=TRUE)

# Sanity check
stopifnot(getPathname(dfC) != getPathname(df))
stopifnot(equals(dfC, df))

# Try to copy it without overwriting nor skipping
fail <- tryCatch({
  copyTo(df, path=path, overwrite=FALSE, skip=FALSE)
  FALSE
}, error = function(ex) { TRUE })
stopifnot(fail)

# Copy it again by overwriting exiting output
dfC <- copyTo(df, path=path, overwrite=TRUE)
print(dfC)
# Sanity checks
stopifnot(getChecksum(dfC) == getChecksum(df))
stopifnot(getPathname(dfC) != getPathname(df))

## renameTo() may fail on some test systems.
## See R.utils Issue #42.
if (FALSE) {
dfR <- renameTo(dfC, getPathname(dfC), verbose=TRUE)
stopifnot(equals(dfR, df))

filenameC <- getFilename(dfC)
filenameR <- sprintf("%s.foo", filenameC)
dfR <- renameTo(dfC, filenameR, verbose=TRUE)
print(dfR)
stopifnot(equals(dfR, df))

dfC2 <- renameTo(dfR, filenameC, verbose=TRUE)
print(dfC2)
stopifnot(getFilename(dfC2) == filenameC)
stopifnot(equals(dfR, df))
}


# Cleanup
file.remove(getPathname(dfC))

# Sanity checks
stopifnot(!isFile(dfC))
stopifnot(isFile(df))


message("*** copyTo() / renameTo() on GenericDataFile ... DONE")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# linkTo()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** linkTo() on GenericDataFile")

# Link to it in a temporary directory
path <- tempdir()

# On Windows, necessary privileges are required.  If not
# available, generate a warning and not an error.
if (.Platform$OS.type == "windows") {
  dfL <- NULL
  tryCatch({
    dfL <- linkTo(df, path=path)
  }, error = function(ex) {
    print(ex)
    cat("The above exception was caught but ignored for this package system test\n")
  })

  if (!is.null(dfL)) {
    print(dfL)

    # Sanity checks
    stopifnot(getChecksum(dfL) == getChecksum(df))

    # Copy file (via link) - unless a Windows Shortcut link
    isWindowsShortcut <- (getPathname(dfL) == getPathname(df))
    if (!isWindowsShortcut) {
      dfLC <- copyTo(dfL, path=file.path(path, "foo"), overwrite=TRUE, validate=FALSE)
      # Sanity checks
      stopifnot(getChecksum(dfLC) == getChecksum(df))
      stopifnot(getPathname(dfLC) != getPathname(df))
      # Cleanup
      file.remove(getPathname(dfLC))
      # Sanity checks
      stopifnot(!isFile(dfLC))
      stopifnot(isFile(dfL))
      stopifnot(isFile(df))
      # Cleanup
      file.remove(getPathname(dfL))
    } else {
      # Done with the Windows Shortcut link
      dfL <- NULL
    }
  }
} else {
  dfL <- linkTo(df, path=path)
  print(dfL)

  # Sanity checks
  stopifnot(getChecksum(dfL) == getChecksum(df))
  stopifnot(getPathname(dfL) != getPathname(df))

  # Copy file (via link)
  dfLC <- copyTo(dfL, path=file.path(path, "foo"), overwrite=TRUE)
  # Sanity checks
  stopifnot(getChecksum(dfLC) == getChecksum(df))
  stopifnot(getPathname(dfLC) != getPathname(df))
  # Cleanup
  file.remove(getPathname(dfLC))
  # Sanity checks
  stopifnot(!isFile(dfLC))
  stopifnot(isFile(dfL))
  stopifnot(isFile(df))

  # Cleanup
  file.remove(getPathname(dfL))
} # if (.Platform$OS.type == "windows")

# Sanity checks
stopifnot(!isFile(dfL))
stopifnot(isFile(df))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Unknown arguments
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
df <- GenericDataFile(mustExist=FALSE, foobar=42L, .onUnknownArgs="ignore")
df <- GenericDataFile(mustExist=FALSE, foobar=42L, .onUnknownArgs="warning")
res <- try(df <- GenericDataFile(mustExist=FALSE, foobar=42L, .onUnknownArgs="error"), silent=TRUE)
stopifnot(inherits(res, "try-error"))

checksum <- getChecksum(df)
print(checksum)
stopifnot(is.na(checksum))


message("*** GenericDataFile ... DONE")
