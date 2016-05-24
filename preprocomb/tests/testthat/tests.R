
modifiediris <- droplevels(iris[-c(1:60),])

## Test DataClass
testdataobject <- initializedataclassobject(modifiediris)
expect_is(testdataobject, "DataClass")
expect_true(all(testdataobject@variance==TRUE, testdataobject@finite==TRUE, testdataobject@completeobs==TRUE, testdataobject@classbalance==TRUE, testdataobject@ntopratiotwoplus==TRUE, testdataobject@mindimensions==TRUE))

## PhaseClass

expect_error(testphase_error <- setphase("testphase_error", "something", preimpute=TRUE))
testphase_success <- setphase("testphase_success", "naomit", preimpute=TRUE)
expect_is(testphase_success, "PhaseClass")

## PreprocessorClass

numpreproccessors <- paste("List of", length(getpreprocessors()))
expect_output(str(testpreprocessors(data=modifiediris)), numpreproccessors)

## GridClass
testgrid <- setgrid(phases=c("outliers"), data=modifiediris)
testgrid2 <- setgrid(phases=c("outliers","scaling"), data=modifiediris)
expect_is(testgrid, "GridClass")

## Preprocombclass

testcomb <- preprocomb(gridclassobject=testgrid, models="knn", nholdout=3, predict=TRUE, cluster=TRUE, outlier=TRUE, search="exhaustive")
testcomb2 <- preprocomb(gridclassobject=testgrid2, models="knn", nholdout=2, predict=FALSE, cluster=TRUE, outlier=TRUE, search="exhaustive")
expect_is(testcomb, "PreProCombClass")
expect_true(!all(testcomb@rawall[,2:7]==0))


