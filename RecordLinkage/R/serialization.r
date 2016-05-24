# functions and methods to clone Record Linkage objects and save them
# to disk (applying to RL* object which store data in a database)

# Clone object, including database.
setGeneric(
  name = "clone",
  def = function(object, ...) standardGeneric("clone")
)

setMethod(
  f = "clone",
  signature = "RLBigData",
  definition = function(object, ...)
  {
    result <- object
    result@pairs <- clone.ffdf(object@pairs)
    result@Wdata <- clone.ff(object@Wdata)
    result@WdataInd <- clone.ff(object@WdataInd)
    result
  }
)

setMethod(
  f = "clone",
  signature = "RLResult",
  definition = function(object, ...)
  {
  }
)


# Save object to disk. File format is a zip file of Rdata and ff files
setGeneric(
  name = "saveRLObject",
  def = function(object, file, ...) standardGeneric("saveRLObject")
)

setMethod(
  f = "saveRLObject",
  signature = "RLBigData",
  definition = function(object, file)
  {
    objectfile <- tempfile(pattern = "object")
    ff_file <- tempfile(pattern = "ff_data")
    
    # save RL object as such
    save(object, file=objectfile)
    # save ff objects
    # have to exist as variables outside of S4 object to be saved
    patterns <- object@pairs
    Wdata <- object@Wdata
    WdataInd <- object@WdataInd
    M <- object@M
    U <- object@U
    ffsave(patterns, Wdata, WdataInd, M, U, file=ff_file)
    zip(file, files = c(objectfile, paste(ff_file, c(".ffData", ".RData"), sep="")))
  }
)

setMethod(
  f = "saveRLObject",
  signature = "RLResult",
  definition = function(object, file)
  {
  }
)

# Load object from disk.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("patterns", "Wdata", 
  "WdataInd", "M", "U"))
loadRLObject <- function(file)
{
  # get file contents and unzip file into temporary directory
  tmpdir <- tempdir()
  tryCatch(
    l <- unzip(file, list = TRUE),
    error = function(e) stop(sprintf("An error occured while opening file %s: %s",
      file, print(e))))
  unzip(file, exdir = tmpdir, junkpaths = TRUE)
  filenames <-  basename(as.character(l$Name))
  # find filename of S4 object among extracted files
  objectFileInd <- grep("object", filenames)
  if(length(objectFileInd)==0) stop("File with object data not found in zip file!")
  if(length(objectFileInd) > 1) stop("Multiple files with object data found in zip file!")
  # load S4 object
  load(paste(tmpdir, filenames[objectFileInd], sep="/"))

  # find filenames of ff objects among extracted files
  # first one suffices, remove extension
  ffFileInd <- grep(".ffData", filenames)
  ffFilename <- paste(tmpdir, sub(".ffData", "", filenames[ffFileInd]), sep="/")
  tryCatch(ffload(ffFilename), error = function(e)
    stop(sprintf("Error while loading ff data: %s", print(e))))

  # assign ff components to S4 object
  object@pairs <- patterns
  object@Wdata <- Wdata
  object@WdataInd <- WdataInd
  object@M <- M
  object@U <- U
  object
}

