##' testthat tests for get_trackdata
##'
context("testing get_trackdata function")

dbName = "ae"

path2orig = file.path(tempdir(), "emuR_demoData", paste0(dbName, emuDB.suffix))

# extract internalVars from environment .emuR_pkgEnv
internalVars = get("internalVars", envir = .emuR_pkgEnv)

ae = load_emuDB(path2orig, inMemoryCache = internalVars$testingVars$inMemoryCache, verbose = F)


# get segmentlist of type segment
path2segl <- list.files(system.file("extdata", package = "emuR"), pattern = glob2rx("ae-n.seg"), full.names = TRUE)
n = read.emusegs(path2segl)

# get segmentlist of type event
path2segl <- list.files(system.file("extdata", package = "emuR"), pattern = glob2rx("ae-hStar.seg"), full.names = TRUE)
hStar = read.emusegs(path2segl)

##############################
test_that("correct classes are returned", {
  
  td = get_trackdata(ae, n, 'fm', verbose=F)
  expect_that(class(td), equals('trackdata'))
  
  td = get_trackdata(ae, n, 'fm', cut=.5, verbose=F)
  expect_that(class(td), equals('data.frame'))
  
  td = get_trackdata(ae, n, 'fm', cut=.5, npoints=3, verbose=F)
  expect_that(class(td), equals('trackdata'))
  
  td = get_trackdata(ae, n, 'fm', cut=.5, npoints=1, verbose=F)
  expect_that(class(td), equals('data.frame'))
  
  td = get_trackdata(ae, hStar, 'fm', verbose=F)
  expect_that(class(td), equals('data.frame'))
  
  td = get_trackdata(ae, hStar, 'fm', npoints=3, verbose=F)
  expect_that(class(td), equals('trackdata'))
  
  sl = query(ae, "Phonetic=@|i:")
  td = get_trackdata(ae, sl, "fm", resultType = "emuRtrackdata", verbose=F)
  expect_true(inherits(td, 'emuRtrackdata'))
})

#############################
test_that("bad calls", {
  expect_error(get_trackdata(ae, verbose=F)) # seglist missing error
  expect_error(get_trackdata(seglist=n, verbose=F)) # S3 error
  expect_error(get_trackdata(ae, n, verbose=F)) # ssffTrackName missing error
  expect_error(get_trackdata(ae, n, cut=2, verbose=F)) # cut > 1 error
  expect_error(get_trackdata(ae, n, cut=-1, verbose=F)) # cut < 0 error
  expect_error(get_trackdata(ae, n, 'fm', npoints=3, verbose=F)) # npoint with no cut argument error
  expect_error(get_trackdata(ae, n, 'fm', onTheFlyParams = formals(open), verbose=F)) # no onTheFlyFunctionName error
  expect_error(get_trackdata(aeDB, n, 'fm', onTheFlyOptLogFilePath = '/path/to/bla/', verbose=F)) # onTheFlyOptLogFilePath error
  expect_error(get_trackdata(ae, n, resultType = "emuRtrackdata", verbose = F)) # bad resultType for seglist of type 'emusegs'
})

##############################
test_that("returned trackdata$data field has correct length", {
  td = get_trackdata(ae, n, 'fm',cut=.5, npoints=3, verbose=F)
  expect_that(dim(td$data)[1], equals(length(n$utts)*3))
  
  td = get_trackdata(ae, n, 'fm',cut=.5, npoints=5, verbose=F)
  expect_that(dim(td$data)[1], equals(length(n$utts)*5))
  
})

##############################
test_that("all sorts of cut values work", {
  cutVals = seq(0, 1, 0.04)
  for(cutV in cutVals){
    td = get_trackdata(ae, n, 'fm',cut=cutV, verbose=F)
    expect_that(class(td), equals('data.frame'))
  }
})

###############################
test_that("n points greater than boundaries word", {
  
    td = get_trackdata(ae, n, 'fm', cut=0.5, npoints = 20, verbose=F)
    expect_that(class(td), equals('trackdata'))
})


#############################
test_that("on-the-fly calculations work", {
  for(wrasspFun in names(wrasspOutputInfos)){
    if(length(wrasspOutputInfos[[wrasspFun]]$tracks) > 0 && wrasspOutputInfos[[wrasspFun]]$outputType == "SSFF"){
      td = get_trackdata(ae, n, wrasspOutputInfos[[wrasspFun]]$tracks[1], onTheFlyFunctionName = wrasspFun, verbose=F)
    }
    expect_that(class(td), equals('trackdata'))
  }
})

##############################
test_that("on-the-fly calculations work if ssffTrackName is not set", {
  td = get_trackdata(ae, n, onTheFlyFunctionName = "ksvF0", verbose=F)
  expect_equal(dim(td$index)[1], 12)
})

##############################
test_that("data fields are the same as hardcoded values (taken from original emu.track(n, 'fm') command)", {
  # note that values have slightly changed due to the recalulation with wrassp
  td = get_trackdata(ae, n, 'fm', verbose=F)
  expect_that(td$data[10,1], equals(256))
  expect_that(td$data[10,2], equals(1521))
  expect_that(td$data[10,3], equals(2382))
  expect_that(td$data[10,4], equals(3573))
  # on-the-fly values should be the same
  td = get_trackdata(ae, n, onTheFlyFunctionName = "forest", verbose=F)
  expect_that(td$data[10,1], equals(256))
  expect_that(td$data[10,2], equals(1521))
  expect_that(td$data[10,3], equals(2382))
  expect_that(td$data[10,4], equals(3573))
  
  
})


