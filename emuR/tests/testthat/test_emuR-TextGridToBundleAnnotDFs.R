##' testthat tests for TextGridToBundleAnnotDFs
##'
context("testing TextGridToBundleAnnotDFs function")

path2demoData = file.path(tempdir(),"emuR_demoData")

path2tg = file.path(path2demoData, "TextGrid_collection", "msajc003.TextGrid")

tgBundleAnnotDFs = TextGridToBundleAnnotDFs(path2tg, sampleRate = 20000, annotates = "msajc003.wav", name = "msajc003")

##############################
test_that("correct SEGMENT values are parsed and calculated in data.frame tables", {  
  
  # get Phonetic table
  phoneticTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Phonetic"))
  
  expect_that(phoneticTbl[1,]$type, equals('SEGMENT'))
  
  # first segment of Phonetic
  # item[0] = {id: XYZ, labels: [{name: ‘lab', value: ‘V'}], sampleStart: 3749, sampleDur: 1389}
  expect_that(phoneticTbl[1,]$sample_start, equals(0))
  
  # second segment
  expect_that(phoneticTbl[2,]$sample_start, equals(3749))
  expect_that(phoneticTbl[2,]$sample_dur, equals(1389))
  tmpItemID = phoneticTbl[2,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('V'))
  
  # 18th segment
  # item[16] = {id: XYZ, labels: [{name: ‘lab', value: ‘@'}], sampleStart: 30124, sampleDur: 844}
  expect_that(phoneticTbl[18,]$sample_start, equals(30124))
  expect_that(phoneticTbl[18,]$sample_dur, equals(844))
  tmpItemID = phoneticTbl[18,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('@'))
  
  # 35th segment
  # item[33] = {id: XYZ, labels: [{name: ‘lab', value: ‘l'}], sampleStart: 50126, sampleDur: 1962}
  expect_that(phoneticTbl[35,]$sample_start, equals(50126))
  expect_that(phoneticTbl[35,]$sample_dur, equals(1962))
  tmpItemID = phoneticTbl[35,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('l'))
  
})

##############################
test_that("correct EVENT values are parsed and calculated in SQLite items table", {
  
  # get Tone table
  toneTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Tone"))
  
  # first event
  # item[0] = {id: XYZ, labels: [{name: ’tone', value: ‘H*'}], samplePoint: 8381}
  expect_that(toneTbl[1,]$sample_point, equals(8381))
  tmpItemID = toneTbl[1,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('H*'))
  
  # 4th event
  # item[3] = {id: XYZ, labels: [{name: ’tone', value: ‘H*'}], samplePoint: 38255}
  expect_that(toneTbl[4,]$sample_point, equals(38255))
  tmpItemID = toneTbl[4,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('H*'))
  
  # 7th event
  # item[6] = {id: XYZ, labels: [{name: ’tone', value: ‘L%'}], samplePoint: 51552}
  expect_that(toneTbl[7,]$sample_point, equals(51552))
  tmpItemID = toneTbl[7,]$item_id
  labelsRow = dplyr::filter_(tgBundleAnnotDFs$labels, ~(item_id==tmpItemID))
  
  expect_that(labelsRow$label, equals('L%'))
  
})  

##############################
test_that("SEGMENTs & EVENTs have correct itemIDs in data.frame tables", {
  
  # get Tone table
  toneTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Tone"))
  
  # get Phonetic table
  phoneticTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Phonetic"))

  # increment IDs for EVENTs
  expect_equal(toneTbl[2,]$item_id, toneTbl[1,]$item_id + 1)
  expect_equal(toneTbl[3,]$item_id, toneTbl[2,]$item_id + 1)
  
  # increment ids for SEGMENTs
  expect_equal(phoneticTbl[2,]$item_id, phoneticTbl[1,]$item_id + 1)
  expect_equal(phoneticTbl[3,]$item_id, phoneticTbl[2,]$item_id + 1)
  
})


##############################
test_that("data.frame labels table has correct values", {
  # get Tone table
  toneTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Tone"))
  
  # get Phonetic table
  phoneticTbl = dplyr::filter_(tgBundleAnnotDFs$items, ~(level=="Phonetic"))
  
  # check phoneticsTable are ok
  expect_equal(phoneticTbl[1,]$item_id, 86)
  expect_equal(sum(phoneticTbl[1,]$labelIdx), 0)
  expect_equal(phoneticTbl[1,]$level, 'Phonetic')
  
  # check toneTbl are ok
  expect_equal(toneTbl[1,]$item_id, 122)
  expect_equal(sum(toneTbl[1,]$labelIdx), 0)
  expect_equal(toneTbl[1,]$level, 'Tone')
  
  # check labelTbl
  labelsTbl = dplyr::filter_(tgBundleAnnotDFs$labels, ~(name=="Phonetic"))
  expect_equal(paste0(labelsTbl$label, collapse = ''), 'VmVNstH@:frEnzSi:w@zkH@nsId@dbju:dH@f@l')
  labelsTbl = dplyr::filter_(tgBundleAnnotDFs$labels, ~(name=="Tone"))
  expect_equal(paste0(labelsTbl$label, collapse = ''), 'H*H*L-H*H*L-L%')
  
})

