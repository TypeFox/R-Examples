context("Input")

non_unique <- data.frame(V1 = c("a", "a"), V2 = c("b", "b"), value = 1:2)
txt <- "{\"dataset\":{\"dimension\":{\"V1\":{\"category\":{\"index\":[\"a\"]}},\"id\":[\"V1\"],\"size\":[1]},\"value\":[1]}}"

test_that("wrong input fails", {
    expect_error(fromJSONstat(1), "is not a character vector")
    expect_error(fromJSONstat(character(0)), "not greater than 0")
    expect_error(fromJSONstat(txt, 1), "is not a string")
    expect_error(fromJSONstat(txt, letters), "is not a string")
    expect_error(fromJSONstat(txt, "a"), "naming must be \"label\" or \"id\"")
    expect_error(fromJSONstat(txt, use_factors = "a"), "is not a flag")
    expect_error(fromJSONstat(txt, use_factors = NA),
                 "contains 1 missing values")
    expect_error(toJSONstat(1), "(?:.*is not a data frame)(?:.* is not a list)")
    expect_error(toJSONstat(fromJSONstat(txt), letters), "is not a string")
    expect_error(toJSONstat(list()), "not greater than 0")
    expect_error(toJSONstat(list(1)), "is not a data frame")
    expect_error(toJSONstat(data.frame()), "not greater than 0")
    expect_error(toJSONstat(data.frame(value = 1)), "not greater than 1")
    expect_error(toJSONstat(data.frame(value = 1, V1 = NA)), "missing values")
    expect_error(toJSONstat(fromJSONstat(txt), "a"),
                 "is not a column in dataset")
    expect_error(toJSONstat(data.frame(value = 1, id = 1)),
                 "not allowed column names")
    expect_error(toJSONstat(data.frame(value = 1, size = 1)),
                 "not allowed column names")
    expect_error(toJSONstat(data.frame(value = 1, role = 1)),
                 "not allowed column names")
    expect_error(toJSONstat(non_unique),
                 "non-value columns must constitute a unique ID")
})

test_that("name of value column works", {
    df1 <- data.frame(V1 = "a", value = 1)
    expect_match(toJSONstat(df1), "\"value\":\\[1\\]")
    expect_match(toJSONstat(df1, value = ""), "\"value\":\\[1\\]")
    df2 <- data.frame(V1 = "a", v = 1)
    expect_match(toJSONstat(df2, value = "v"), "\"value\":\\[1\\]")
})

test_that("round-trip works", {
    df1 <- fromJSONstat(txt)
    df2 <- fromJSONstat(toJSONstat(df1))
    expect_equal(df1, df2)
})
