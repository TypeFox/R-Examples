context("UNFv6: UNF Common Tests")

test_that("Numerics", {
    expect_equal(unf6(0)$unf, 
                 unf6("+0.e+")$unf, 
                 "YUvj33xEHnzirIHQyZaHow==")
    expect_equal(unf6(1)$unf, 
                 unf6("+1.e+")$unf, 
                 "tv3XYCv524AfmlFyVOhuZg==")
    expect_equal(unf6(-300)$unf, 
                 unf6("-3.e+2")$unf, 
                 "ZTXyg54FoMfRDWZl6oWmFQ==")
    expect_equal(unf6(3.1415)$unf, 
                 unf6("+3.1415e+")$unf, 
                 "vOSZmXXXpKfQcqZ0Cuu5/w==")
    expect_equal(unf6(0.00073)$unf, 
                 unf6("+7.3e-4")$unf, 
                 "qhw3qzg3fEK0NNfoVxk4jQ==")
    expect_equal(unf6(0.0123456789012345)$unf, 
                 unf6("+1.234568e-2")$unf, 
                 "tbCwhb/PKbElzXByo6hLAw==")
})

test_that("Nonfinites", {
    expect_equal(unf6(NaN)$unf, 
                 unf6("+nan")$unf, 
                 "GNcR8/UCnImaPpw47gdPNg==")
    expect_equal(unf6(Inf)$unf, 
                 unf6("+inf")$unf, 
                 "MdAI70WZdDHnu6qmkpqUQg==")
    expect_equal(unf6(-Inf)$unf, 
                 unf6("-inf")$unf, 
                 "A7orv3pgAhljFnGjQVLCog==")
})

test_that("Character strings", {
    expect_equal(unf6("A character String")$unf, "FYqU7uBl885eHMbpco1ooA==", label="Short character strings")
    expect_equal(unf6("A quite long character string, so long that the number of characters in it happens to be more than the default cutoff limit of 128.")$unf,
                 unf6("A quite long character string, so long that the number of characters in it happens to be more than the default cutoff limit of 1")$unf,
                 "/BoSlfcIlsmQ+GHu5gxwEw==", label = "Long character string")
    #x1 <- "på Færøerne"
    #expect_equal(unf6(x1)$unf, "KHM6bKVaVaxWDDsmyerfDA==", label = "Non-ASCII characters in string 1")
    #x2 <- "p\303\245 F\303\246r\303\270erne"
    #expect_equal(unf6(x2)$unf, "KHM6bKVaVaxWDDsmyerfDA==", label = "Non-ASCII characters in string 2")
    #x3 <- `Encoding<-`(x2, "UTF-8")
    #expect_equal(unf6(x3)$unf, "KHM6bKVaVaxWDDsmyerfDA==", label = "Non-ASCII characters in string 3")
    expect_equal(unf6("")$unf, "ECtRuXZaVqPomffPDuOOUg==", label = "Empty character string")
    expect_equal(unf6(NA)$unf, "cJ6AyISHokEeHuTfufIqhg==", label = "Missing value")
})

test_that("Datetimes", {
    expect_equal(unf6("2014-01-13T20:47:18")$unf, "eaMxex5EHi2LunomVc0SDw==")
    expect_equal(unf6(strptime("2014-01-13T20:47:18", format = "%Y-%m-%dT%H:%M:%S"))$unf, "eaMxex5EHi2LunomVc0SDw==")
    expect_equal(unf6("2014-01-14T01:47:18Z")$unf, "1Pku/Z/EIRtmpdEepAb1MA==")
    expect_equal(unf6(strptime("2014-01-14T01:47:18Z", format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"), timezone = "UTC")$unf, "1Pku/Z/EIRtmpdEepAb1MA==")
})

test_that("Boolean/logical values", {
    expect_equal(unf6(TRUE)$unf, unf6("+1.e+")$unf, "tv3XYCv524AfmlFyVOhuZg==")
    expect_equal(unf6(FALSE)$unf, unf6("+0.e+")$unf, "YUvj33xEHnzirIHQyZaHow==")
})
