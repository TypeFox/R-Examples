##' testthat test to see if generated logfile stays the same
##'
##' @author Raphael Winkelmann

context("test logging capabilities")

test_that("logging file content is the same as hard coded string", {
  
  str = "beginTime : 0 \ncenterTime : FALSE \nendTime : 0 \nwindowShift : 5 \nwindowSize : 20 \neffectiveLength : TRUE \nwindow : BLACKMAN \nanalysisOrder : 0 \nenergyNormalization : FALSE \nlengthNormalization : FALSE \ntoFile : FALSE \nexplicitExt :  \noutputDirectory :  \nforceToLog : TRUE \n => on files:"
  
  altDir = tempdir()
  path2log = paste(altDir, "/wrasspTESTTHATlog.txt", sep="")
  
  if(file.exists(path2log)){
    unlink(path2log)
  }

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  for (func in names(wrasspOutputInfos)){
    for(wavFile in wavFiles[1]){
      funcFormals = formals(func)
      funcFormals$listOfFiles = wavFile
      funcFormals$toFile = FALSE
      funcFormals$forceToLog = TRUE
      funcFormals$optLogFilePath = path2log
      res = do.call(func,as.list(funcFormals))
    }
  }
  lines = suppressWarnings(readLines(path2log))
  logFileStr = paste(lines[6:20],collapse="\n")

  expect_that(logFileStr, equals(str))

  expect_that(grep("######## zcrana performed ########", lines), equals(232))

  blackGrep = grep("window : BLACKMAN ", lines)== c(12,  61,  78,  96, 160, 180, 201)

  expect_that(sum(blackGrep), equals(7))

  unlink(path2log)

})
