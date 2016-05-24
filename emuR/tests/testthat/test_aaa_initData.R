context("init data for testing")
# AAA in the name so it is run as first test
path2demoData = file.path(tempdir(),"emuR_demoData")
path2testhatFolder = file.path(tempdir(),"emuR_testthat")

unlink(path2demoData, recursive = T)
unlink(path2testhatFolder, recursive = T)

create_emuRdemoData(precache = T)
create_BPFcollectionManipulated(path2demoData)

dir.create(path2testhatFolder)
