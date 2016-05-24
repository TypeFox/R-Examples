##' testthat tests for convert_TextGridCollection
##'

requireNamespace("RSQLite", quietly = T)

context("testing convert_TextGridCollection function")

path2demoData = file.path(tempdir(), "emuR_demoData")
path2testData = file.path(tempdir(), "emuR_testthat")
path2tgCol = file.path(path2demoData, "TextGrid_collection")

emuDBname = 'convert-TextGridCollection-testDB'

path2newDb = file.path(path2testData, paste0(emuDBname, emuDB.suffix))


# clean up
unlink(path2newDb, recursive = T)

##############################
test_that("bad calls cause errors", {
  
  # create dir
  dir.create(path2newDb)
  
  # existing targetDir causes errors
  expect_error(convert_TextGridCollection(dir = path2tgCol, 
                                          dbName = emuDBname,
                                          targetDir = path2testData, 
                                          verbose=F), regexp = "already exists!", ignore.case = T)
  # clean up
  unlink(path2newDb, recursive = T)
  
})

##############################
test_that("correct emuDB is created", {
  
  convert_TextGridCollection(dir = path2tgCol, 
                             dbName = emuDBname,
                             path2testData, verbose=F)
  
  test_that("emuDB has correct file format on disc", {
    # 2 files in top level
    tmp = list.files(path2newDb)
    expect_equal(length(tmp), 2)
    
    # 14 files in 0000_ses
    tmp = list.files(file.path(path2newDb,'0000_ses'), recursive = T)
    expect_equal(length(tmp), 14)
  })
  
  test_that("emuDB _DBconfig.json is correct", {
    # read config
    dbCfgJSONLns = readLines(file.path(path2newDb, paste0(emuDBname, '_DBconfig.json')), warn = FALSE)
    dbCfgJSON = paste(dbCfgJSONLns,collapse='')
    dbCfgPersisted = jsonlite::fromJSON(dbCfgJSON, simplifyVector=FALSE)
    
    # correct name
    expect_equal(dbCfgPersisted$name, emuDBname)
    # no ssffTrackDefs
    expect_equal(length(dbCfgPersisted$ssffTrackDefinitions), 0)
    # no linkDefs
    expect_equal(length(dbCfgPersisted$linkDefinitions), 0)
    # levelDef stuff
    expect_equal(length(dbCfgPersisted$levelDefinitions), 11)
    expect_equal(dbCfgPersisted$levelDefinitions[[9]]$name, 'Phonetic')
    
    # EMUwebAppConfig stuff
    expect_equal(length(dbCfgPersisted$EMUwebAppConfig$perspectives), 1)
    expect_equal(dbCfgPersisted$EMUwebAppConfig$perspectives[[1]]$signalCanvases$order[[1]], 'OSCI')
    expect_equal(length(dbCfgPersisted$EMUwebAppConfig$perspectives[[1]]$levelCanvases$order), 11)
    
  })
  
  test_that("emuDB _annot.json is correct", {
    # read annot
    annotJSONLns = readLines(file.path(path2newDb, '0000_ses/msajc003_bndl/msajc003_annot.json'), warn = FALSE)
    annotJSON = paste(annotJSONLns,collapse='')
    annotPersisted = jsonlite::fromJSON(annotJSON,simplifyVector=FALSE)
    # general stuff
    expect_equal(annotPersisted$name, 'msajc003')
    expect_equal(annotPersisted$annotates, 'msajc003.wav')
    expect_equal(length(annotPersisted$links), 0)
    expect_equal(length(annotPersisted$levels), 11)
    expect_equal(annotPersisted$levels[[9]]$name, 'Phonetic')
    
    # test a couple of items
    
    # second segment
    expect_that(annotPersisted$levels[[9]]$items[[2]]$sampleStart, equals(3749))
    expect_that(annotPersisted$levels[[9]]$items[[2]]$sampleDur, equals(1389))
    expect_that(annotPersisted$levels[[9]]$items[[2]]$labels[[1]]$value, equals('V'))
    
    # 18th segment
    expect_that(annotPersisted$levels[[9]]$items[[18]]$sampleStart, equals(30124))
    expect_that(annotPersisted$levels[[9]]$items[[18]]$sampleDur, equals(844))
    expect_that(annotPersisted$levels[[9]]$items[[18]]$labels[[1]]$value, equals('@'))
    
    # 35th segment
    # item[33] = {id: XYZ, labels: [{name: ‘lab', value: ‘l'}], sampleStart: 50126, sampleDur: 1962}
    expect_that(annotPersisted$levels[[9]]$items[[35]]$sampleStart, equals(50126))
    expect_that(annotPersisted$levels[[9]]$items[[35]]$sampleDur, equals(1962))
    expect_that(annotPersisted$levels[[9]]$items[[35]]$labels[[1]]$value, equals('l'))
    
  })
  
  # clean up
  unlink(path2newDb, recursive = T)
  
})

##############################
test_that("only specified tiers are converted when tierNames is set", {
  
  convert_TextGridCollection(dir = path2tgCol, 
                             dbName = emuDBname,
                             path2testData, tierNames=c("Phonetic", "Tone"), verbose=F)
  
  test_that("emuDB has correct file format on disc", {
    # 2 files in top level
    tmp = list.files(path2newDb)
    expect_equal(length(tmp), 2)
    
    # 14 files in 0000_ses
    tmp = list.files(file.path(path2newDb,'0000_ses'), recursive = T)
    expect_equal(length(tmp), 14)
  })
  
  test_that("emuDB _DBconfig.json is correct", {
    # read config
    dbCfgJSONLns=readLines(file.path(path2newDb, paste0(emuDBname, '_DBconfig.json')),warn=FALSE)
    dbCfgJSON=paste(dbCfgJSONLns,collapse='')
    dbCfgPersisted=jsonlite::fromJSON(dbCfgJSON,simplifyVector=FALSE)
    
    # correct name
    expect_equal(dbCfgPersisted$name, emuDBname)
    # no ssffTrackDefs
    expect_equal(length(dbCfgPersisted$ssffTrackDefinitions), 0)
    # no linkDefs
    expect_equal(length(dbCfgPersisted$linkDefinitions), 0)
    # levelDef stuff
    expect_equal(length(dbCfgPersisted$levelDefinitions), 2)
    expect_equal(dbCfgPersisted$levelDefinitions[[1]]$name, 'Phonetic')
    expect_equal(dbCfgPersisted$levelDefinitions[[1]]$type, 'SEGMENT')
    expect_equal(dbCfgPersisted$levelDefinitions[[2]]$name, 'Tone')
    expect_equal(dbCfgPersisted$levelDefinitions[[2]]$type, 'EVENT')
    
    # EMUwebAppConfig stuff
    expect_equal(length(dbCfgPersisted$EMUwebAppConfig$perspectives), 1)
    expect_equal(dbCfgPersisted$EMUwebAppConfig$perspectives[[1]]$signalCanvases$order[[1]], 'OSCI')
    expect_equal(length(dbCfgPersisted$EMUwebAppConfig$perspectives[[1]]$levelCanvases$order), 2)
    
  })
  
  test_that("emuDB _annot.json is correct", {
    # read annot
    annotJSONLns=readLines(file.path(path2newDb, '0000_ses/msajc003_bndl/msajc003_annot.json'),warn=FALSE)
    annotJSON=paste(annotJSONLns,collapse='')
    annotPersisted=jsonlite::fromJSON(annotJSON,simplifyVector=FALSE)
    # general stuff
    expect_equal(annotPersisted$name, 'msajc003')
    expect_equal(annotPersisted$annotates, 'msajc003.wav')
    expect_equal(length(annotPersisted$links), 0)
    expect_equal(length(annotPersisted$levels), 2)
    expect_equal(annotPersisted$levels[[1]]$name, 'Phonetic')
    expect_equal(annotPersisted$levels[[1]]$type, 'SEGMENT')
    expect_equal(annotPersisted$levels[[2]]$name, 'Tone')
    expect_equal(annotPersisted$levels[[2]]$type, 'EVENT')
    # test a couple of items
    
    # second segment
    expect_that(annotPersisted$levels[[1]]$items[[2]]$sampleStart, equals(3749))
    expect_that(annotPersisted$levels[[1]]$items[[2]]$sampleDur, equals(1389))
    expect_that(annotPersisted$levels[[1]]$items[[2]]$labels[[1]]$value, equals('V'))
    
    # 18th segment
    expect_that(annotPersisted$levels[[1]]$items[[18]]$sampleStart, equals(30124))
    expect_that(annotPersisted$levels[[1]]$items[[18]]$sampleDur, equals(844))
    expect_that(annotPersisted$levels[[1]]$items[[18]]$labels[[1]]$value, equals('@'))
    
    # 35th segment
    # item[33] = {id: XYZ, labels: [{name: ‘lab', value: ‘l'}], sampleStart: 50126, sampleDur: 1962}
    expect_that(annotPersisted$levels[[1]]$items[[35]]$sampleStart, equals(50126))
    expect_that(annotPersisted$levels[[1]]$items[[35]]$sampleDur, equals(1962))
    expect_that(annotPersisted$levels[[1]]$items[[35]]$labels[[1]]$value, equals('l'))
    
  })
  
  
  # clean up
  unlink(path2newDb, recursive = T)
  
})
