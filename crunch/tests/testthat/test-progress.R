context("Polling progress")

with_mock_HTTP({
    test_that("If progress polling gives up, it tells you what to do", {
        with(temp.option(crunch.timeout=0.5), {
            expect_error(pollProgress("/api/progress/1.json", wait=0.25),
            paste('Your process is still running on the server. It is',
                'currently 22.5% complete. Check',
                '`uncached(crGET("/api/progress/1.json"))` until it reports',
                '100% complete'),
            fixed=TRUE)
        })
    })

    counter <- 1
    with_mock(
        ## GET something slightly different each time through so we can
        ## approximate polling a changing resource
        `httr::GET`=function (url, ...) {
            if (is.null(url)) {
                stop("No URL found", call.=FALSE)
            }
            url <- paste0(url, counter, ".json") ## Add counter
            counter <<- counter + 1 ## Increment
            url <- sub("^\\/", "", url) ## relative to cwd
            out <- handleShoji(fromJSON(url, simplifyVector=FALSE))
            return(list(
                status_code=200,
                times=structure(nchar(url), .Names="total"),
                request=list(method="GET", url=url),
                response=out
            ))
        },
        test_that("Progress polling goes until 100", {
            expect_identical(pollProgress("/api/progress/", wait=.05), 100)
        })
    )
})
