###
### $Id: tictoc.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.toc <- function(delay, expected) {
    Sys.sleep(delay)
    output <- matlab::toc(FALSE)
    all.equal(output,
              expected,
              tolerance = 0.10)
}

matlab::tic()
try(test.toc(4, 4))
try(test.toc(2, 4+2))

