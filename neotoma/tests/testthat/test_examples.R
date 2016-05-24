## load packages
library("testthat")
library("neotoma")

context("Run Neotoma examples only when not on CRAN")

test_that("Examples run without error", {
    ## we don't want this to run on CRAN
    skip_on_cran()

    ## List of example topics we want to check
    egs <- c('compile_downloads',
             'compile_taxa',
             'counts',
             'get_chroncontrol',
             'get_contact',
             'get_dataset',
             'get_download',
             #'get_geochron',
             'get_publication',
             'get_site',
             'get_table',
             'get_taxa',
             'write_agefile')

    refnames <- paste0("example-ref-", egs, ".rds")

    for (i in seq_along(egs)) {
        egout <- try(example(topic = egs[i], package = "neotoma", ask = FALSE,
                             character.only = TRUE, run.dontrun = TRUE,
                             echo = TRUE))
        expect_that(inherits(egout, "try-error"), is_false(),
                    label = paste("Error raised in example:", egs[i]))
        ## expect_that(egout, equals_reference(refnames[i]))
    }
})
