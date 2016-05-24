##' testthat tests for create_filePairList
##'
context("testing create_filePairList function")


path2demoData = file.path(tempdir(), "emuR_demoData")
path2testData = file.path(tempdir(), "emuR_testthat")
path2tgCol = file.path(path2demoData, "TextGrid_collection")


ext1 = 'wav'
ext2 = 'TextGrid'

wavPaths = list.files(path2tgCol, pattern=paste(ext1, "$", sep = ""), recursive=T, full.names=T)
tgPaths = list.files(path2tgCol, pattern=paste(ext2, "$", sep = ""), recursive=T, full.names=T)

testDirName = 'test_createFilePairList'

path2testDir = file.path(path2testData, testDirName)

##############################
test_that("bad calls cause errors", {
  
  expect_error(create_filePairList('asdf', '', '', ''), 'ext1Path2rootDir does not exist!')
  expect_error(create_filePairList(path2tgCol, 'asdf', '', ''), 'ext2Path2rootDir does not exist!')
  
})

##############################
test_that("error is generated when nr of ext1 files > ext2 files", {
  # create testdir
  dir.create(path2testDir)
  
  # copy files 
  file.copy(wavPaths, path2testDir)
  file.copy(tgPaths[-length(tgPaths)], path2testDir)
  
  expect_error(create_filePairList(path2testDir, path2testDir, 'wav', 'TextGrid'))
  
  # clean up
  unlink(path2testDir, recursive = T)
  
})


##############################
test_that("correct filePairList is generated when nr of ext1 files < ext2 files", {
  # create testdir
  dir.create(path2testDir)
  
  # copy files 
  file.copy(wavPaths[-length(wavPaths)], path2testDir)
  file.copy(tgPaths, path2testDir)
  
  fpl = create_filePairList(path2testDir, path2testDir, 'wav', 'TextGrid')
  
  expect_equal(dim(fpl)[1], 6)
  expect_equal(dim(fpl)[2], 2)
  
  # clean up
  unlink(path2testDir, recursive = T)
  
})

##############################
test_that("error is thrown if dirs are empty", {
  
  # create testdir
  dir.create(path2testDir)
  
  expect_error(create_filePairList(path2testDir, path2testDir, 'wav', 'TextGrid'))
    
  # clean up
  unlink(path2testDir, recursive = T)
  
})

##############################
test_that("error is thrown if one ext2 does not have same base name", {
  
  # create testdir
  dir.create(path2testDir)
  
  # copy files 
  file.copy(wavPaths, path2testDir)
  file.copy(tgPaths, path2testDir)
  
  #rename file
  file.rename(file.path(path2testDir, basename(tgPaths[3])), file.path(path2testDir, 'asdf.TextGrid'))
  
  expect_error(create_filePairList(path2testDir, path2testDir, 'wav', 'TextGrid'))
  
  # clean up
  unlink(path2testDir, recursive = T)
  
})

