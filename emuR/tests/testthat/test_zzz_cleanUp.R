context("clean up test data")

unlink(file.path(tempdir(),"emuR_demoData"), recursive = T)
unlink(file.path(tempdir(),"emuR_testthat"), recursive = T)