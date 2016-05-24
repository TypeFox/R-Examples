context("test makeProblem()")

sp <- searchpaths()
fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
microData <- get(load(fn))
dim.region <- data.frame(levels=c('@','@@','@@','@@','@@'), codes=c('Total', 'A','B','C','D'), stringsAsFactors=FALSE)
dim.gender <- data.frame(levels=c('@','@@','@@'), codes=c('Total', 'male','female'), stringsAsFactors=FALSE)
dimList <- list(region=dim.region, gender=dim.gender)
dimVarInd <- c(1,2)
freqVarInd <- numVarInd <- weightInd <- sampWeightInd <- NULL

problem <- makeProblem(
  data=microData,
  dimList=dimList,
  dimVarInd=dimVarInd,
  freqVarInd=freqVarInd,
  numVarInd=numVarInd,
  weightInd=weightInd,
  sampWeightInd=sampWeightInd)

expect_that(problem, is_a("sdcProblem"))
expect_that(get.problemInstance(problem@problemInstance, "nrVars"), equals(15))

