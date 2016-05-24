# Test reading of Packed Experiment (pxp) files
# 
# Author: jefferis
###############################################################################

context("Test handling of Igor pxp files")

pxp<-read.pxp(system.file("igor/WedJul407c2_001.pxp", package = 'IgorR'))

test_that("Read Igor packed experiment file", {
      
      expected_names<-c("vars", "WavSelect", "ChanSelect", "ChanWaveList", "yLabel", 
          "Group", "Set1", "Set2", "SetX", "CT_TimeStamp", "CT_TimeIntvl", 
          "RecordA0", "RecordB0", "RecordA1", "RecordB1", "RecordA2", "RecordB2", 
          "RecordA3", "RecordB3", "RecordA4", "RecordB4", "RecordA5", "RecordB5", 
          "RecordA6", "RecordB6", "RecordA7", "RecordB7", "RecordA8", "RecordB8", 
          "RecordA9", "RecordB9", "RecordA10", "RecordB10", "RecordA11", 
          "RecordB11", "RecordA12", "RecordB12", "RecordA13", "RecordB13", 
          "RecordA14", "RecordB14", "RecordA15", "RecordB15", "RecordA16", 
          "RecordB16", "RecordA17", "RecordB17", "RecordA18", "RecordB18", 
          "RecordA19", "RecordB19", "RecordA20", "RecordB20", "RecordA21", 
          "RecordB21", "RecordA22", "RecordB22", "RecordA23", "RecordB23", 
          "RecordA24", "RecordB24", "RecordA25", "RecordB25", "RecordA26", 
          "RecordB26", "RecordA27", "RecordB27", "RecordA28", "RecordB28", 
          "RecordA29", "RecordB29", "RecordA30", "RecordB30", "RecordA31", 
          "RecordB31", "RecordA32", "RecordB32", "RecordA33", "RecordB33", 
          "RecordA34", "RecordB34", "RecordA35", "RecordB35", "RecordA36", 
          "RecordB36", "RecordA37", "RecordB37", "RecordA38", "RecordB38", 
          "RecordA39", "RecordB39", "RecordA40", "RecordB40", "RecordA41", 
          "RecordB41", "RecordA42", "RecordB42", "RecordA43", "RecordB43", 
          "RecordA44", "RecordB44", "RecordA45", "RecordB45", "RecordA46", 
          "RecordB46", "RecordA47", "RecordB47", "RecordA48", "RecordB48", 
          "RecordA49", "RecordB49", "RecordA50", "RecordB50", "RecordA51", 
          "RecordB51", "RecordA52", "RecordB52", "RecordA53", "RecordB53", 
          "RecordA54", "RecordB54", "RecordA55", "RecordB55", "RecordA56", 
          "RecordB56", "RecordA57", "RecordB57", "RecordA58", "RecordB58", 
          "RecordA59", "RecordB59", "RecordA60", "RecordB60", "RecordA61", 
          "RecordB61", "RecordA62", "RecordB62", "RecordA63", "RecordB63", 
          "RecordA64", "RecordB64", "RecordA65", "RecordB65", "RecordA66", 
          "RecordB66", "RecordA67", "RecordB67", "RecordA68", "RecordB68", 
          "RecordA69", "RecordB69", "RecordA70", "RecordB70", "RecordA71", 
          "RecordB71", "RecordA72", "RecordB72", "RecordA73", "RecordB73", 
          "RecordA74", "RecordB74", "RecordA75", "RecordB75", "RecordA76", 
          "RecordB76", "RecordA77", "RecordB77", "RecordA78", "RecordB78", 
          "RecordA79", "RecordB79", "RecordA80", "RecordB80", "RecordA81", 
          "RecordB81", "RecordA82", "RecordB82", "RecordA83", "RecordB83", 
          "RecordA84", "RecordB84", "RecordA85", "RecordB85", "RecordA86", 
          "RecordB86", "RecordA87", "RecordB87", "RecordA88", "RecordB88", 
          "RecordA89", "RecordB89", "RecordA90", "RecordB90", "RecordA91", 
          "RecordB91", "RecordA92", "RecordB92", "RecordA93", "RecordB93", 
          "RecordA94", "RecordB94", "RecordA95", "RecordB95", "RecordA96", 
          "RecordB96", "RecordA97", "RecordB97", "RecordA98", "RecordB98", 
          "RecordA99", "RecordB99", "MulticlampVCForRAccess", "ChanA", 
          "ChanB", "Notes")
      
      expect_that(names(pxp),
          equals(expected_names),'pxp contents identified')
      first5vars<-structure(list(sysVars = c(-2.2097864151001, -71.0904388427734, 
                  0.0737364292144775, -17.2362651824951, 0.504302203655243, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), FileFormat = 1.91, 
              FileDateTime = 3266426893, NumWaves = 100, TotalNumWaves = 200), .Names = c("sysVars", 
              "FileFormat", "FileDateTime", "NumWaves", "TotalNumWaves"))
      expect_that(pxp$vars[1:5],
          equals(first5vars),'read in some Igor variables')
      
      expect_that(range(pxp[['RecordA0']]),
          equals(c(-207, 159.4375)),'read in some numeric wave data')
    })

test_that("Igor Wave to R time series", {
      RecordA0<-pxp[['RecordA0']]
      wA0=WaveToTimeSeries(RecordA0)
      expect_that(wA0,is_a('ts'))
      
      expect_that(tsp.igorwave(RecordA0),
          is_equivalent_to(c(0, 19.975, 40)),'check correct tsp')
      
      expect_that(tsp.igorwave(RecordA0),
          is_equivalent_to(tsp(wA0)),'check tsp attributes from wave and time series match')
    })

test_that("Convert Igor Wave to R time series when reading pxp file", {
      expect_warning(pxp<-read.pxp(system.file("igor/WedJul407c2_001.pxp", package = 'IgorR'),
                    regex='Record',ReturnTimeSeries=TRUE), 
                    regexp = 'object must have one or more observations')
      record_names = paste("Record",
          rep(c("A","B"),99), rep(0:99,rep(2,100)), sep="")
      RecordA0<-pxp[['RecordA0']]
      expect_that(RecordA0,is_a('ts'))
      # restrict to waves only
      pxp=pxp[!sapply(pxp,inherits,'list')]
      expect_that(names(pxp), equals(record_names))
    })

test_that("Read pxp file loading only waves matching regex", {
      pxp<-read.pxp(system.file("igor/WedJul407c2_001.pxp", package = 'IgorR'),
                    regex='Record')
      record_names = paste("Record",
          rep(c("A","B"),99), rep(0:99,rep(2,100)), sep="")
      # restrict to waves only
      pxp=pxp[!sapply(pxp,inherits,'list')]
      expect_that(names(pxp), equals(record_names))
    })

test_that("Read pxp file structure only", {
      pxp<-read.pxp(system.file("igor/WedJul407c2_001.pxp", package = 'IgorR'),
                    StructureOnly = TRUE)
      record_names = paste("Record",
          rep(c("A","B"),99), rep(0:99,rep(2,100)), sep="")
      expect_true( all(sapply(pxp[record_names],is.na)) )
      
      # stored structure
      pxp.str.baseline<-readRDS('igor/WedJul407c2_001_str.rds')
      expect_equivalent(pxp, pxp.str.baseline)
    })

test_that("Read pxp files containing variables with higher characters", {
      pxp<-read.pxp("igor/ExperimentWithHigherChars.pxp")
      expect_that(pxp$vars$myscalar,equals(1))
      expect_that(pxp$vars$mystring,equals("Hello!"))
      # should be the same on windows 1252, latin-1 and macintosh 
      expect_that(pxp$vars$micron,equals(iconv("\u00B5",from='UTF-8',to="UTF-8")))
      # different for windows 1252, latin-1 and macintosh 
      expect_that(pxp$vars$pix6,
                  equals(iconv("\u03C0\u03C0\u03C0\u03C0\u03C0\u03C0",
                               from='UTF-8',to='UTF-8')))
    })

test_that("Read pxp history", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_false(is.null(pxp$history))
      expect_true(nchar(pxp$history) > 0)
    })

test_that("Read pxp recreation macro", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_false(is.null(pxp$recmacro))
      expect_true(nchar(pxp$recmacro) > 0)
    })

test_that("Read pxp main procedure", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_false(is.null(pxp$mainproc))
      expect_equal(pxp$mainproc, "#pragma rtGlobals=3\t\t// Use modern global access method and strict wave access.\r\n")
    })

test_that("Read pxp embedded procedure", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_false(is.null(pxp$Proc0))
      expect_equal(pxp$Proc0, "#pragma rtGlobals=3\t\t// Use modern global access method and strict wave access.\r\n")
    })

test_that("Read pxp plain notebook", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_false(is.null(pxp$Notebook0))
      expect_equal(pxp$Notebook0, "plain notebook!")
    })

test_that("Read pxp ignores formatted notebook", {
      pxp<-read.pxp("igor/ExperimentWithProcHistoryAndNBs.pxp", ExtractText=T)
      expect_true(is.null(pxp$Notebook1))
    })

test_that("Read pxp can cope with invalid folder names", {
  pxp<-read.pxp("igor/bad_foldername4.pxp")
  expect_true("_default_" %in% names(pxp))
})
