##' testthat tests for emuRtrackdata
##'
context("testing emuRtrackdata functions")

dbName = "ae"

path2orig = file.path(tempdir(), "emuR_demoData", paste0(dbName, emuDB.suffix))
path2testData = file.path(tempdir(), "emuR_testthat")
path2db = file.path(path2testData, paste0(dbName, emuDB.suffix))

# extract internalVars from environment .emuR_pkgEnv
internalVars = get("internalVars", envir = .emuR_pkgEnv)

# delete, copy and load
unlink(path2db, recursive = T)
file.copy(path2orig, path2testData, recursive = T)
ae = load_emuDB(path2db, inMemoryCache = internalVars$testingVars$inMemoryCache, verbose = F)


##############################
test_that("correct classes are returned", {
  
  sl = query(ae, "Phonetic=@|i:")
  td = get_trackdata(ae, 
                     seglist = sl, 
                     ssffTrackName = 'fm', verbose = F)
  
  newTd = create_emuRtrackdata(sl, td)
  
  expect_true(inherits(newTd, "emuRtrackdata"))
  
})

##############################
# test_that("cut works correctly", {
#   
#   sl = query(dbName, "Phonetic=@|i:")
#   td = get_trackdata(dbName, 
#                      seglist = sl, 
#                      ssffTrackName = 'fm')
#   
#   newTd = create_emuRtrackdata(sl, td)
#   
#   propRes = cut_td(newTd, 0.5, prop=T)
#   print(propRes)
#   })