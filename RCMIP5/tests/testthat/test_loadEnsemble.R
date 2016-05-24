# Testing code for the RCMIP5 'loadEnsemble.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code:
#   source("loadEnsemble.R")
#   library(testthat)
#   test_file("tests/testthat/test_loadEnsemble.R")

context("loadEnsemble")

test_that("loadEnsemble handles bad input", {
    expect_error(loadEnsemble("","","","","",path="does_not_exist"))  # path does not exist
    expect_error(loadEnsemble("","","","","",path=c("path1","path2")))       # multi-value path
    expect_error(loadEnsemble("","","","","",path=1))                  # non-character path
    expect_error(loadEnsemble(variable=1,"","","",""))                          # non-character
    expect_error(loadEnsemble("",model=1,"","",""))                          # non-character
    expect_error(loadEnsemble("","",experiment=1,"",""))                          # non-character
    expect_error(loadEnsemble("","","",ensemble=1,""))                          # non-character
    expect_error(loadEnsemble("","","","",domain=1))                # non-character
    expect_error(loadEnsemble(variable=c("",""),"","","",""))                   # multi-value
    expect_error(loadEnsemble("",model=c("",""),"","",""))                   # multi-value
    expect_error(loadEnsemble("","",experiment=c("",""),"",""))                   # multi-value
    expect_error(loadEnsemble("","","",ensemble=c("",""),""))                   # multi-value
    expect_error(loadEnsemble("","","","",domain=c("","")))         # multi-value
    expect_error(loadEnsemble("","","","",verbose=1))               # non-logical verbose
    expect_error(loadEnsemble("","","","",recursive=1))             # non-logical recursive
})

test_that("loadEnsemble handles no files found", {            # no NetCDF files found
    w <- getOption('warn')
    options(warn=-1)
    expect_warning(loadEnsemble("","","","","", path=("testdata_none")))
    expect_is(loadEnsemble("","","","","", path=("testdata_none")), "NULL")
    options(warn=w)
})

test_that("loadEnsemble loads monthly data", {
    
    skip_on_cran()
    
    path <- "../../sampledata/monthly"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- loadEnsemble('prc','GFDL-CM3', 'rcp85', 'r1i1p1', '[^_]+', path=path, verbose=F)
    expect_is(d, "cmip5data")
    d <- loadEnsemble('prc','GFDL-CM3','rcp85','r1i1p1','[^_]+', path=path, verbose=F)     
    expect_is(d, "cmip5data")
    expect_equal(length(d$files), 1)                                 # should be one file
})

test_that("loadEnsemble loads annual data", {
    
    skip_on_cran()
    
    path <- "../../sampledata/annual"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- loadEnsemble('co3', 'HadGEM2-ES', 'rcp85', 'r1i1p1', '[^_]+', path=path, verbose=F)
    expect_is(d, "cmip5data")
})

test_that("loadEnsemble loads 4D data", {
    
    skip_on_cran()
    
    path <- "../../sampledata/annual"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- loadEnsemble('ph','MPI-ESM-LR','historical','r1i1p1', '[^_]+', path=path, verbose=F)     # test data set
    expect_is(d, "cmip5data")
    expect_is(d$Z, "array")
    expect_is(d$val, "array")
    
    d <- loadEnsemble('co3','HadGEM2-ES','rcp85','r1i1p1', '[^_]+', 
                      path=path, verbose=F)     # test data set
    expect_is(d,"cmip5data")
    expect_is(d$Z, "array")
    
    path <- "../../sampledata/monthly"
    d <- loadEnsemble('tsl','GFDL-CM3','historicalGHG','r1i1p1', '[^_]+', 
                      path=path, verbose=F)     # test data set
    expect_is(d,"cmip5data")
    expect_is(d$Z, "array")
})

test_that("loadEnsemble checks unique domain", {
    expect_error(loadEnsemble("co3","fakemodel1-ES","rcp85","r1i1p1", '[^_]+',
                              path='testdata_twodomains/', verbose=F))
})

test_that("loadEnsemble assigns ancillary data", {
    
    skip_on_cran()
    
    path <- "../../sampledata/annual"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- loadEnsemble('co3','HadGEM2-ES','rcp85','r1i1p1', '[^_]+', path=path,verbose=F)
    expect_is(d, "cmip5data")
    expect_true(!is.null(d$provenance))
})

test_that("loadEnsemble handles 2D lon and lat", {
    
    skip_on_cran()
    
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- loadEnsemble('tos','GFDL-ESM2G', 'historical', 'r1i1p1', '[^_]+', path=path, verbose=F)
    expect_is(d, "cmip5data")
    expect_is(d$lon, "numeric")
    expect_is(d$lat, "numeric")
    expect_null(dim(d$lon))
    expect_null(dim(d$lat))    
})

test_that("loadEnsemble handles data with time length=1", {

    skip_on_cran()
    
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    # This is a real CMIP5 file with one single month
    # loadEnsemble should add an extra dimension (of length 1) to the data
    d <- loadEnsemble("spco2", "HadGEM2-ES", "rcp85", domain="Omon",
                      ensemble="r1i1p1", path=path, verbose=F)    
    expect_is(d, "cmip5data")
    # TODO: check dimensions
})


test_that("loadEnsemble handles time-only data", {
    
    skip_on_cran()
    
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
        
    # This is a real CMIP5 file with no lon or lat, just time
    # loadEnsemble should read OK, and add extra lon/lat dimensions of length 1
    d <- loadEnsemble("co2mass", "GFDL-ESM2M", "historical", domain="Amon", ensemble="r1i1p1",
                      path=path, verbose=F)
    
    expect_is(d, "cmip5data")
    # TODO: check dimensions
})
