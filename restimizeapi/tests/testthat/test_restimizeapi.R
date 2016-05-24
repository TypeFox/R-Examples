test_that (
    ".CheckParams year is null and function stops",
    {
        expect_error (.CheckParams (year=NULL))
    }
)

test_that (
    ".CheckParams year is character string and function stops",
    {
        expect_error (.CheckParams (year="hello"))
    }
)

test_that (
    ".CheckParams quarter is character string and function stops",
    {
        expect_error (.CheckParams (year=1999, quarter="hello"))
    }
)

test_that (
    ".CheckParams year is valid and quarter is NOT between one and four",
    {
        expect_error (.CheckParams (year=1999, quarter=0))
        expect_error (.CheckParams (year=1999, quarter=5))
    }
)

test_that (
    ".CheckParams year is valid and quarter is between one and four.",
    {   
        .CheckParams (year=1999, quarter=1)
        .CheckParams (year=1999, quarter=4)
    }
)

test_that (
    "GetCompanyReleases rejects a call where a ticker and only a quarter is provided.",
    {   
        expect_error (GetCompanyReleases ("MSFT", quarter=3))
    }
)

test_that (
    "GetCompanyEstimates rejects a call where a ticker and only a quarter is provided.",
    {   
        expect_error (GetCompanyEstimates ("MSFT", quarter=3))
    }
)

test_that (
    "GetEstimates rejects null startDate.",
    {   
    expect_error (GetEstimates (startDate = NULL, endDate = "2015-05-20"))
    }
)

test_that (
    "GetEstimates rejects null startDate.",
    {   
        expect_error (GetEstimates (startDate = "2015-05-20", endDate = NULL))
    }
)

test_that (
    "GetReleaseConsensus rejects null ids.",
    {
        expect_error (GetGetReleaseConsensus (id = NULL))
    }
)