test_that("empty lists return empty df", {
    cases <- data.frame(id=c(), icd9list=c())
    expect_that(melt_icd9list(cases, "id", "icd9list"), equals(data.frame()))
})

test_that("numeric input to icd9list", {
    cases <- data.frame(id=c(1,2), icd9list=c(100.1,200.1))
    expect_that(melt_icd9list(cases, "id", "icd9list"), 
                equals(data.frame(id=c(1,2), icd9cm=c("1001","2001"))))
})

test_that("real icd9 codes", {
    cases <- data.frame(id=c(1,2), icd9list=c('162.4,070.30,155.0,401.9','996.52,E878.8,V45.86'))
    expect_that(melt_icd9list(cases, "id", "icd9list"), 
                equals(data.frame(id=c(1, 1, 1, 1, 2, 2, 2), 
                                  icd9cm=c("1624", "07030", "1550", "4019", "99652", "E8788", "V4586"))))
})
