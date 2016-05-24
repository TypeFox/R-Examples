context("UNFv6: Character Encoding")

test_that("Characters with known encoding coerced to UTF-8", {
    d <- data.frame(x = c("\u00E6","\u00E5","\u00F8"), stringsAsFactors = FALSE)
    d_latin1 <- data.frame(x = sapply(d$x, "Encoding<-", "latin1"), stringsAsFactors = FALSE)
    d_utf8 <- data.frame(x = sapply(d$x, "Encoding<-", "UTF-8"), stringsAsFactors = FALSE)
    expect_equal(unf(d)$unf, unf(d_utf8)$unf)
    expect_equal(unf(d)$unf, unf(d_latin1)$unf)
    rm("d")
    rm("d_latin1")
    rm("d_utf8")
})

test_that("ASCII characters unaffected by encoding during file I/O", {
    d <- data.frame(x = c("a","b","c"), stringsAsFactors = FALSE)
    write.csv(d, file="temp.csv", row.names=FALSE, fileEncoding="UTF-8")
    d2 <- read.csv("temp.csv", fileEncoding="latin1")
    expect_equal(unf(d)$unf, unf(d2)$unf)
    unlink("temp.csv")
    rm("d")
    rm("d2")
})

test_that("non-ASCII characters unaffected by correct UTF-8 encoding during file I/O", {
    d <- data.frame(x = c("\u00E6","\u00E5","\u00F8"), stringsAsFactors = FALSE)
    write.csv(d, file="temp.csv", row.names=FALSE, fileEncoding = "UTF-8")
    d2 <- read.csv("temp.csv", encoding = "UTF-8", stringsAsFactors = FALSE)
    expect_equal(unf(d)$unf, unf(d2)$unf)
    unlink("temp.csv")
    rm("d2")
})

test_that("non-ASCII characters unaffected by incorrect encoding during file I/O", {
    d <- data.frame(x = c("\u00E6","\u00E5","\u00F8"), stringsAsFactors = FALSE)
    write.csv(d, file="temp.csv", row.names=FALSE, fileEncoding = "UTF-8")
    d2 <- read.csv("temp.csv", encoding = "latin1", stringsAsFactors = FALSE)
    #expect_equal(unf(d)$unf, unf(d2)$unf)
    unlink("temp.csv")
    rm("d")
    rm("d2")
})
