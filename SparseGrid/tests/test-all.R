library('SparseGrid')

if(require("testthat", quietly = TRUE)) {
    test_package("SparseGrid")
} else {
    print( "package 'testthat' not available, cannot run unit tests" )
}
