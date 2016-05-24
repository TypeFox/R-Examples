##' testthat test to see if writing to then reading from disc
##' changes anything 
##'
##' @author Raphael Winkelmann
context("test fileIO")

test_that("read things that are written to disc are the same as origs", {
  
  altDir = tempdir()
  
  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)
  
  for (func in names(wrasspOutputInfos)){
    for(wavFile in wavFiles){
      funcFormals = formals(func)
      funcFormals$listOfFiles = wavFile
      funcFormals$outputDirectory = altDir
      funcFormals$explicitExt = "testthat"
      res = do.call(func,as.list(funcFormals))
      path2new = paste(altDir, basename(wavFile), sep="/")
      sp=unlist(strsplit(path2new, ".", fixed = T))
      sp[length(sp)] = funcFormals$explicitExt
      fromFile = read.AsspDataObj(paste(sp, collapse="."))
      funcFormals$toFile = FALSE
      inMem = do.call(func,as.list(funcFormals))
      # test attributes if they are the same 
      for (at in names(attributes(fromFile))){
        if(at != "filePath"){
          expect_that(length(attr(fromFile, at)), equals(sum(attr(fromFile, at)==attr(inMem, at))))
        }
      }
      # test if data is the same
      expect_that(sum(unlist(inMem[attr(inMem,"names")]) == unlist(fromFile[attr(fromFile,"names")]))
                  , equals(length(unlist(inMem[attr(inMem,"names")]))))
    }
  }
  
  # clean up
  files <- list.files(altDir, "testthat$", full.names=T)
  
  for (file in files){
    unlink(file)
  }
  
})
