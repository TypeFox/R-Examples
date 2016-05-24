##' testthat tests for validation of bundles
##'
context("testing validate.XXX.bundle functions")

dbName = "ae"

path2orig = file.path(tempdir(), "emuR_demoData", paste0(dbName, emuDB.suffix))
path2testData = file.path(tempdir(), "emuR_testthat")
path2db = file.path(path2testData, paste0(dbName, emuDB.suffix))

# extract internalVars from environment .emuR_pkgEnv
internalVars = get("internalVars", envir = .emuR_pkgEnv)

#################################
test_that("unaltered bundle (sqlTableRep) validates successfully", {
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = internalVars$testingVars$inMemoryCache, verbose = F)
  
  res = validate_bundleDBI(ae, session = "0000", bundle = "msajc003")
  expect_equal(res$type, 'SUCCESS')
  expect_equal(res$message, '')
})


