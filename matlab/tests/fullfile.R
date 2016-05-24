###
### $Id: fullfile.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.fullfile <- function(input, expected) {
    output <- do.call(getFromNamespace("fullfile", "matlab"), input)
    identical(output, expected)
}

fullfile.expected <- file.path(path.expand("~"), "somedir", "foo.txt")

test.fullfile(list(dir    = path.expand("~"),
                   subdir = "somedir",
                   file   = "foo.txt"), fullfile.expected)

