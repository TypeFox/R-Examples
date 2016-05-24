context("SR research test files")

test_that("test_files_load",{
files <- c('bino1000.asc.gz','bino250.asc.gz','bino500.asc.gz','binoRemote250.asc.gz','binoRemote500.asc.gz','mono1000.asc.gz','mono2000.asc.gz','mono250.asc.gz','mono500.asc.gz','monoRemote250.asc.gz','monoRemote500.asc.gz')
for  (f in files)
{
    fpath <- system.file(paste0("extdata/",f),package="eyelinker")
    tst <- read.asc(fpath)
    expect_equal(sort(names(tst)),c("blinks", "fix","info","msg","raw", "sacc"  ))
}
})

