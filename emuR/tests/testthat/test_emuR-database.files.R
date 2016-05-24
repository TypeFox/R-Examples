##' testthat tests for database.files
##'
context("testing database.files functions")

# extract internalVars from environment .emuR_pkgEnv
internalVars = get("internalVars", envir = .emuR_pkgEnv)

#######################################
test_that("file operations work", {
  
  dbName = 'ae'
  
  path2orig = file.path(tempdir(), "emuR_demoData", paste0(dbName, emuDB.suffix))
  path2testData = file.path(tempdir(), "emuR_testthat")
  path2db = file.path(path2testData, paste0(dbName, emuDB.suffix))
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  
  ae = load_emuDB(path2db, inMemoryCache = internalVars$testingVars$inMemoryCache, verbose = F)
  
  # extract internalVars from environment .emuR_pkgEnv
  internalVars = get("internalVars", envir = .emuR_pkgEnv)
  
  test_that("import_mediaFiles works", {
    wavPath = system.file('extdata', package='wrassp')
    import_mediaFiles(ae, dir = wavPath, targetSessionName = 'newSes', verbose = F)
    expect_true(file.exists(file.path(path2db, 'newSes_ses')))
    paths = list.files(file.path(path2db, 'newSes_ses'), recursive = T, full.names = T, pattern = 'wav$')
    expect_equal(length(paths), 9)
    paths = list.files(file.path(path2db, 'newSes_ses'), recursive = T, full.names = T, pattern = '_annot.json$')
    expect_equal(length(paths), 9)
  })
  
  test_that("CRUD operations for files work", {
    
    test_that("add = (C)RUD", {
      wrasspExtdataPath = system.file('extdata', package='wrassp')
      wavFilePaths = list.files(wrasspExtdataPath, pattern = "wav$", full.names = T, recursive = T)
      
      outDirPath = file.path(path2testData, 'zcranaVals')
      dir.create(outDirPath)
      zcrana(wavFilePaths, outputDirectory = outDirPath)
      
      add_files(ae, dir = outDirPath, fileExtension = 'zcr', targetSessionName = 'newSes')
      zcrPaths = list.files(path2db, pattern = 'zcr$', recursive = T)
      expect_equal(length(zcrPaths), 9)
      
      # cleanup
      unlink(outDirPath, recursive = T)
      
      
    })
    
    test_that("list = C(R)UD", {
      df = list_files(ae)
      expect_equal(dim(df),c(55, 3))
    })
    
    
  })

})
