traitRateParams <- function(mainType, n, fn, outDir, outFile,
                            charModelParam1, charModelParam2,
                            gammaParam, seqModelParam1){
  
  ## assemble control file
  ## ---------------------
  ctrl <- c(paste("_mainType", mainType),
            paste("_treeFile", fn[1]),
            paste("_characterFile", fn[3]),
            paste("_seqFile", fn[2]),
            paste("_outDir", outDir),
            paste("_outFile", outFile),
            "_logFile log.txt",
            "_scaledTreeFile scaled.tree",
            paste("_charModelParam1", charModelParam1),
            paste("_charModelParam2", charModelParam2),
            paste("_gammaParam", gammaParam),
            paste("_seqModelParam1", seqModelParam1),
            "_relRate 1",
            "_seqModelType HKY",
            "_logValue 3",
            "_bScaleTree 1",
            "_stochasicMappingIterations 100",
            "_treeLength 1.0")
  
  ## number of iterations
  if ( mainType == "runTraitBootstrap" ){
    if ( missing(n) ) n <- 200
    ctrl <- c(ctrl,
              paste("_", c("start", "end"),
                    "SimulationsIter ", c(0, n - 1), sep = ""))
  }
  write(ctrl, file = "params.txt")
}