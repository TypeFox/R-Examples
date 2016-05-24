context("Logging and caching")

if (run.integration.tests) {
    with(test.authentication, {
        with(temp.option(httpcache.log=""), {
            msg <- capture.output(z <- crGET(getOption("crunch.api")))
            msg <- strsplit(msg, " ")
            expect_identical(length(msg), 1L)
            expect_identical(msg[[1]][2], "CACHE")
            expect_identical(msg[[1]][3], "HIT")
        })
    })
}
