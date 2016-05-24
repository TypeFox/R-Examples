##### Helper Functions ##############################################
getSpecs <- function() {
    specs <- list()
    specs[['sp1']] <- data.frame(widths = c(1, 2, 3),
                                  col.names = c('Col1', 'Col2', 'Col3'))
    specs[['sp2']] <- data.frame(widths = c(3, 2, 1),
                                  col.names = c('C1', 'C2', 'C3'))
    specs
}

selector <- function(line, specs) {
    s <- substr(line, 1, 1)
    spec_name <- ''
    if (s == '1')
        spec_name <- 'sp1'
    else if (s == '2')
        spec_name <- 'sp2'

    spec_name
}

emptySelector <- function(line, specs) {
    ''
}
###################################################################

test_that("empty spec returns an empty list", {
    ff <- tempfile()
    cat(file = ff, '123456', '287654', '987654', '198765', sep = '\n')

    myspecs <- list()

    data <- read.multi.fwf(ff, multi.specs = myspecs, select = selector)
    unlink(ff)

    expect_equal(length(data), 0)
})

test_that("empty selector returns an empty list", {
    ff <- tempfile()
    cat(file = ff, '123456', '287654', '987654', '198765', sep = '\n')

    myspecs <- getSpecs()

    data <- read.multi.fwf(ff, multi.specs = myspecs, select = emptySelector)
    unlink(ff)

    expect_equal(length(data), 2)
})

test_that("data and spec lists have the same length", {
    ff <- tempfile()
    cat(file = ff, '123456', '287654', '987654', '198765', sep = '\n')

    myspecs <- getSpecs()

    data <- read.multi.fwf(ff, multi.specs = myspecs, select = selector)
    unlink(ff)

    expect_equal(length(data), length(myspecs))
})

test_that("a file with no matches gives a data list with null elements", {
    ff <- tempfile()
    cat(file = ff, 'lalala', 'dododo', 'hehehe', 'hohoho', sep = '\n')

    myspecs <- getSpecs()

    data <- read.multi.fwf(ff, multi.specs = myspecs, select = selector)
    unlink(ff)

    for (i in data) {
        expect_null(i)
    }
})

test_that("2 'sp1' and 1 'sp2'", {
    ff <- tempfile()
    cat(file = ff, '123456', '287654', '987654', '198765', sep = '\n')

    myspecs <- getSpecs()

    data <- read.multi.fwf(ff, multi.specs = myspecs, select = selector)
    unlink(ff)

    expect_equal(nrow(data$sp1), 2)
    expect_equal(nrow(data$sp2), 1)
})
