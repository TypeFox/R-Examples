# Simple Run-Time Tests for the read.properties function
# 
# Author: Tobias Verbeke
###############################################################################

require(properties)

# test 0: simple test of three properties
propFilePath <- system.file("testFiles", "test0.properties", package = "properties")
(props <- read.properties(file = propFilePath))
expectedProps <- list(key1 = "value1", key2 = "value2", key3 = "value3")
stopifnot(all.equal(props, expectedProps))

(props <- read.properties(file = propFilePath, fields = c("key2", "key3")))
stopifnot(all.equal(props, expectedProps[-1]))

# test 2: values with spaces, underscores etc.
propFilePath <- system.file("testFiles", "test2.properties", package = "properties")
(props <- read.properties(file = propFilePath))
expectedProps <- list(reportFile = "lncap_2010-10-22.pdf", experimentName = "LnCap",
        sweaveFile = "platereaderTemplate.Rnw", rScript = "platereaderScript.R",
        reportAuthor = "George Orwell (Big Brother Inc.)")
stopifnot(all.equal(props, expectedProps))
    
(props <- read.properties(file = propFilePath, fields = c("reportAuthor", "rScript")))
stopifnot(all.equal(props, expectedProps[c(5, 4)]))

# test 3: regression test for more than 5 properties
propFilePath <- system.file("testFiles", "test3.properties", package = "properties")
(props <- read.properties(file = propFilePath))
expectedProps <- list(key1 = "value1", key2 = "value2", key3 = "value3", key4 = "value4", 
    key5 = "value5", key6 = "value6", key7 = "value7")
stopifnot(all.equal(props, expectedProps))

(props <- read.properties(file = propFilePath, fields = c("key2", "key7")))
stopifnot(all.equal(props, expectedProps[c(2, 7)]))

# test 4: empty property values
propFilePath <- system.file("testFiles", "test4.properties", package = "properties")
(props <- read.properties(file = propFilePath))
expectedProps <- list(key1 = "value1", key2 = "", key3 = "value3", key4 = "")
stopifnot(all.equal(props, expectedProps))

(props <- read.properties(file = propFilePath, fields = c("key2", "key3")))
stopifnot(all.equal(props, expectedProps[2:3]))

