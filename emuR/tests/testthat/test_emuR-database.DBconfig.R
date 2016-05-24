##' testthat tests for database.DBconfig
##'
context("testing database.DBconfig functions")

dbName = 'ae'
useInMemoryCache = F

path2orig = file.path(tempdir(), "emuR_demoData", paste0(dbName, emuDB.suffix))
path2testData = file.path(tempdir(), "emuR_testthat")
path2db = file.path(path2testData, paste0(dbName, emuDB.suffix))


##############################
test_that("get_levelDefinition returns correct levelDef", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  #########################
  # get dbObj
  dbConfig = load_DBconfig(ae)
  
  ld = get_levelDefinition(ae, 'Phonetic')
  expect_equal(ld$name, 'Phonetic')
  expect_equal(ld$type, 'SEGMENT')
  expect_equal(ld$attributeDefinitions[[1]]$name, 'Phonetic')
  expect_equal(ld$attributeDefinitions[[1]]$type, 'STRING')
})

##############################
test_that("CRUD operations work for ssffTrackDefinitions", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  test_that("add = (C)RUD", {
    expect_error(add_ssffTrackDefinition(ae, 'fm'))
    expect_error(add_ssffTrackDefinition(ae, 'fm', 'bla'))
    expect_error(add_ssffTrackDefinition(ae, 'newTrackName', 'badColName', 'pit', 
                                         onTheFlyFunctionName = 'mhsF0', interactive = T))
    
    add_ssffTrackDefinition(ae, 'newTrackName', 'pitch', 'pit', 
                            onTheFlyFunctionName = 'mhsF0', interactive = F)
    
    pitFilePaths = list.files(path2db, pattern = 'pit$', recursive = T)
    expect_equal(length(pitFilePaths), 7)
    
  })
  
  test_that("list = C(R)UD", {
    df = list_ssffTrackDefinitions(ae)
    expect_equal(df$name, c('dft','fm', 'newTrackName'))
    expect_equal(df$columnName, c('dft','fm', 'pitch'))
    expect_equal(df$fileExtension, c('dft','fms', 'pit'))
  })
  
  
  test_that("remove = CRU(D)", {
    # bad name
    expect_error(remove_ssffTrackDefinition(ae, name="asdf"))
    remove_ssffTrackDefinition(ae, name="newTrackName", deleteFiles = T)
    # check that _DBconfig entry is deleted
    dbConfig = load_DBconfig(ae)
    expect_equal(dbConfig$ssffTrackDefinitions[[1]]$name, "dft")
    expect_equal(dbConfig$ssffTrackDefinitions[[2]]$name, "fm")
    
    # check that files have been deleted
    filePaths = list_files(ae, "pit")
    expect_equal(nrow(filePaths), 0)
    
  })
  
})

##############################
test_that("CRUD operations work for levelDefinitions", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  
  test_that("add = (C)RUD", {
    expect_error(add_levelDefinition(ae, 'Phonetic', 'SEGM')) # bad type
    expect_error(add_levelDefinition(ae, 'Phonetic', 'SEGMENT')) # already exists
    
    add_levelDefinition(ae, 'Phonetic2', 'SEGMENT')
    
    dbConfig = load_DBconfig(ae)
    expect_equal(length(dbConfig$levelDefinitions), 10)
    
  })
  
  test_that("list = C(R)UD", {
    df = list_levelDefinitions(ae)
    expect_equal(as.vector(df$name[8:10]), c('Tone','Foot', 'Phonetic2'))
    expect_equal(as.vector(df$type[8:10]), c('EVENT','ITEM', 'SEGMENT'))
    expect_equal(as.vector(df$nrOfAttrDefs[1:4]), c(1, 1, 1, 3))
  })
  
  test_that("remove = CRU(D)", {
    
    expect_error(remove_levelDefinition(ae, name="asdf")) # bad name
    expect_error(remove_levelDefinition(ae, name="Phonetic")) # linkDef present

    DBI::dbGetQuery(ae$connection, paste0("INSERT INTO session VALUES ('", ae$UUID,
                                          "', '0001')")) # add item

    DBI::dbGetQuery(ae$connection, paste0("INSERT INTO bundle VALUES ('", ae$UUID,
                                          "', '0001', 'fakeBundle', 'fakeBundle.wav', 20000, '785c7cdb6d4bd5e8b5cd7c56a5946ddf')")) # add item
    
    DBI::dbGetQuery(ae$connection, paste0("INSERT INTO items VALUES ('", ae$UUID,
                                     "', '0001', 'fakeBundle', 1, 'Phonetic2', 'ITEM', 20000, 1, NULL, NULL, NULL)")) # add item
    
    expect_error(remove_levelDefinition(ae, name="Phonetic2")) # item present
    
    DBI::dbGetQuery(ae$connection, paste0("DELETE FROM items WHERE db_uuid='", 
                                     ae$UUID,"'")) # items present
    
    remove_levelDefinition(ae, name="Phonetic2")
    dbConfig = load_DBconfig(ae)
    expect_equal(length(dbConfig$levelDefinition), 9)
    
  })
  
})  

##############################
test_that("CRUD operations work for attributeDefinitions", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  
  test_that("add = (C)RUD", {
    expect_error(add_attributeDefinition(ae, 'Word', 'Word')) # present attrDef
    
    add_attributeDefinition(ae, 'Word', 'testAttrDef')
    df = list_attributeDefinitions(ae, 'Word')
    expect_true('testAttrDef' %in% df$name)
  })
  
  test_that("list = C(R)UD", {
    df = list_attributeDefinitions(ae, 'Word')
    expect_equal(df$name, c('Word', 'Accent', 'Text', 'testAttrDef'))
    expect_equal(df$type, c('STRING', 'STRING', 'STRING', 'STRING'))
    expect_equal(df$hasLabelGroups, c(F, F, F, F))
    expect_equal(df$hasLegalLabels, c(F, F, F, F))
  })
  

  test_that("remove = CRU(D)", {
    expect_error(remove_attributeDefinition(ae, 'Word', 'Word'))
    expect_error(remove_attributeDefinition(ae, 'Word', 'Accent'))
    remove_attributeDefinition(ae, 'Word', 'testAttrDef')
    df = list_attributeDefinitions(ae, 'Word')
    expect_equal(nrow(df), 3)
  })
  
})  

##############################
test_that("CRUD operations work for legalLabels", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  test_that("set = (C)RUD", {
    # non character vector causes error:
    expect_error(set_legalLabels(ae, 
                    levelName = 'Word', 
                    attributeDefinitionName = 'Word',
                    legalLabels=c(1:3)))
    
    
    set_legalLabels(ae,
                    levelName = 'Word',
                    attributeDefinitionName = 'Word',
                    legalLabels=c('A', 'B', 'C'))
  })
  
  test_that("get = C(R)UD", {
    ll = get_legalLabels(ae, 
                         levelName = 'Word', 
                         attributeDefinitionName = 'Word')
    
    expect_equal(ll, c('A', 'B', 'C'))
  })
  

  test_that("remove = CRU(D)", {
    remove_legalLabels(ae, 
                       levelName = 'Word', 
                       attributeDefinitionName = 'Word')
    
    ll = get_legalLabels(ae, 
                         levelName = 'Word', 
                         attributeDefinitionName = 'Word')
    
    expect_true(is.na(ll))
  })
  
})  

##############################
test_that("CRUD operations work for labelGroups", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  test_that("add = (C)RUD", {
    # bad call already def. labelGroup
    expect_error(add_attrDefLabelGroup(ae,
                                       levelName = 'Phoneme', 
                                       attributeDefinitionName = 'Phoneme',
                                       labelGroupName = 'vowel',
                                       labelGroupValues = c('sdf', 'f')))
    
    add_attrDefLabelGroup(ae,
                          levelName = 'Word', 
                          attributeDefinitionName = 'Word',
                          labelGroupName = 'newGroup',
                          labelGroupValues = c('sdf', 'f'))
    
  })
  
  test_that("list = C(R)UD", {
    df = list_attrDefLabelGroups(ae,
                                 levelName = 'Utterance', 
                                 attributeDefinitionName = 'Utterance')
    expect_equal(nrow(df), 0)
    
    df = list_attrDefLabelGroups(ae,
                                 levelName = 'Phoneme', 
                                 attributeDefinitionName = 'Phoneme')
    expect_equal(nrow(df), 6)
    expect_true(df[6,]$values == "H")
    
    df = list_attrDefLabelGroups(ae,
                                 levelName = 'Word', 
                                 attributeDefinitionName = 'Word')
    expect_true(df[1,]$name == "newGroup")
    expect_true(df[1,]$values == "sdf; f")
  })
  
  test_that("remove = CRU(D)", {
    expect_error(remove_attrDefLabelGroup(ae,
                                          levelName = 'Word', 
                                          attributeDefinitionName = 'Word',
                                          labelGroupName = 'notThere'))
    
    remove_attrDefLabelGroup(ae,
                             levelName = 'Word', 
                             attributeDefinitionName = 'Word',
                             labelGroupName = 'newGroup')
    
    df = list_attrDefLabelGroups(ae,
                                 levelName = 'Word', 
                                 attributeDefinitionName = 'Word')
    expect_equal(nrow(df), 0)
  })
  
})  

##############################
test_that("CRUD operations work for linkDefinitions", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  test_that("add = (C)RUD", {
    # bad call (bad type)
    expect_error(add_linkDefinition(ae, "ONE_TO_TWO"))
    # bad call (link exists)
    expect_error(add_linkDefinition(ae, "ONE_TO_ONE", 
                                    superlevelName ="Syllable", 
                                    sublevelName = "Tone"))
    # bad call undefined superlevelName 
    expect_error(add_linkDefinition(ae, "ONE_TO_MANY", 
                                    superlevelName ="undefinedLevel", 
                                    sublevelName = "Tone"))
    
    
    add_linkDefinition(ae, "ONE_TO_MANY", 
                       superlevelName ="Phoneme", 
                       sublevelName = "Tone")
    
  })
  
  test_that("list = C(R)UD", {
    df = list_linkDefinitions(ae)
    expect_equal(ncol(df), 3)
    expect_equal(nrow(df), 10)
    expect_true(df$type[10] == "ONE_TO_MANY")
    expect_true(df$superlevelName[10] == "Phoneme")
    expect_true(df$sublevelName[10] == "Tone")
  })
  
  test_that("remove = CRU(D)", {
    # bad call -> bad superlevelName
    expect_error(remove_linkDefinition(ae, 
                                       superlevelName ="badName", 
                                       sublevelName = "Tone"))
    # bad call -> bad sublevelName
    expect_error(remove_linkDefinition(ae, 
                                       superlevelName ="Word", 
                                       sublevelName = "badName"))
    # bad call -> links present
    expect_error(remove_linkDefinition(ae, 
                                       superlevelName ="Syllable", 
                                       sublevelName = "Tone"))
    
    remove_linkDefinition(ae, 
                          superlevelName ="Phoneme", 
                          sublevelName = "Tone")
    
    df = list_linkDefinitions(ae)
    expect_equal(ncol(df), 3)
    expect_equal(nrow(df), 9)
    
  })
  
})  


##############################
test_that("CRUD operations work for labelGroups", {
  
  # delete, copy and load
  unlink(path2db, recursive = T)
  file.copy(path2orig, path2testData, recursive = T)
  ae = load_emuDB(path2db, inMemoryCache = useInMemoryCache, verbose = F)
  
  
  test_that("add = (C)RUD", {
    add_labelGroup(ae, 
                   name = 'testLG',
                   values = c('a', 'b', 'c'))  
  })
  
  test_that("list = C(R)UD", {
    df = list_labelGroups(ae)
    expect_true(df$name == 'testLG')
    expect_true(df$values =='a; b; c')
  })
  
  test_that("remove = CRU(D)", {
    # bad call -> bad name
    expect_error(remove_labelGroup(ae, 
                                   name = 'badName'))
    
    remove_labelGroup(ae, 
                      name = 'testLG')
    df = list_labelGroups(ae)
    expect_equal(nrow(df), 0)
  })
})  

# 
test_that("delete", {
  unlink(path2db, recursive = T)
})


